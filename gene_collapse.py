"""
gene_collapse.py
----------------
Score and rank genes for de novo variant enrichment using weighted mutation rate models.

For each protein-coding gene, this script computes a p-value reflecting whether the
observed de novo variants in a cohort carry unexpectedly high pathogenicity scores,
given the gene's background mutation rate. Three complementary tests are run:

    pvalue             -- score enrichment tested within the protein only
    pvalue_proteome    -- score enrichment tested both within the protein and
                         against the proteome-wide score distribution
    pvalue_only_proteome -- score enrichment tested against the proteome-wide
                         score distribution only

Scores are drawn from a user-specified variant effect model (e.g. popEVE) and
normalized relative to the 90th percentile of the proteome-wide distribution before
use as weights in a Poisson simulation.

Usage
-----
    python gene_collapse.py \\
        --genes       genes_info.txt \\
        --all-rates   all_rates.txt \\
        --variants    ddd_variants.txt \\
        --outdir      results/ \\
        --study       DDD \\
        --model       popEVE

Example with all optional arguments:
    python gene_collapse.py \\
        --genes       genes_info.txt \\
        --all-rates   all_rates.txt \\
        --variants    ddd_variants.txt \\
        --outdir      results/ \\
        --study       DDD \\
        --norm        75 \\
        --model       popEVE \\
        --pos-path \\
        --rounding    4 \\
        --n-samples   10000

Input file formats (all tab-separated)
---------------------------------------
genes_info.txt   : One row per protein. Required columns:
                     protein        -- protein identifier (must match variants file)
                     gene           -- HGNC gene symbol
                     chrom          -- chromosome (use 'X' for X-linked genes)
                     rates_filepath -- path to the per-protein rates .csv file

all_rates.txt    : Proteome-wide mutation rates. Required columns:
                     prob           -- per-site mutation probability
                     <model>        -- pathogenicity score for the chosen model

ddd_variants.txt : Observed de novo variants. Required columns:
                     study          -- cohort label (e.g. 'DDD', 'ASD_Affected')
                     protein        -- protein identifier
                     mutant         -- mutation identifier (must match rates files)
                     variant_id     -- unique variant identifier with subject id (used for deduplication)
                     <model>        -- pathogenicity score for the chosen model

Output
------
A DataFrame saved to:
    <outdir>/<study>_<model>_<norm>_<rounding>_<n_samples>.csv

Columns: protein, s_obs, n, n_exp, s_exp, pvalue, pvalue_proteome,
         pvalue_only_proteome, gene
"""

import argparse
import warnings
import multiprocessing as mp

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute gene-level enrichment p-values for de novo variant scores "
            "using a Poisson-weighted simulation framework."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required inputs
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--genes", required=True, metavar="FILE",
        help="Tab-separated gene/protein info file (e.g. genes_info.txt). "
             "Must contain columns: protein, gene, chrom, rates_filepath.",
    )
    required.add_argument(
        "--all-rates", required=True, metavar="FILE",
        help="Tab-separated proteome-wide mutation rates file (e.g. all_rates.txt). "
             "Must contain columns: prob, <model>.",
    )
    required.add_argument(
        "--variants", required=True, metavar="FILE",
        help="Tab-separated observed de novo variants file (e.g. ddd_variants.txt). "
             "Must contain columns: study, protein, mutant, variant_id, <model>.",
    )
    required.add_argument(
        "--outdir", required=True, metavar="DIR",
        help="Directory where the output pickle file will be written.",
    )
    required.add_argument(
        "--study", required=True,
        choices=["ASD_Affected", "ASD_Unaffected", "DDD"],
        help="Cohort to analyse. Controls which cohort size constants are used.",
    )

    # Optional / tuning arguments
    parser.add_argument(
        "--norm", type=float, default=75.0, metavar="FLOAT",
        help="Percentile used for score normalization (scores below this percentile "
             "are set to zero).",
    )
    parser.add_argument(
        "--model", default="popEVE", metavar="STR",
        help="Column name of the variant effect model scores to use.",
    )
    parser.add_argument(
        "--pos-path", action="store_true", default=False,
        help="If set, restrict analysis to variants with a positive pathogenicity score.",
    )
    parser.add_argument(
        "--rounding", type=int, default=4, metavar="INT",
        help="Number of decimal places to round normalized scores to. "
             "Finer rounding increases precision but slows simulation.",
    )
    parser.add_argument(
        "--n-samples", type=int, default=10000, metavar="INT",
        help="Number of Monte Carlo samples used to estimate each p-value. "
             "Higher values give more precise p-values at greater runtime cost.",
    )
    parser.add_argument(
        "--workers", type=int, default=4, metavar="INT",
        help="Number of parallel worker processes.",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Cohort size constants
# ---------------------------------------------------------------------------

# Total number of individuals per cohort (used to scale per-site mutation rates
# to cohort-level expected counts).
COHORT_SIZES = {
    "ASD_Affected":   20493,
    "ASD_Unaffected":  7743,
    "DDD":            31058,
}

# Sex-stratified cohort sizes (needed for X-chromosome rate adjustment).
COHORT_SEX_SIZES = {
    "ASD_Affected":  {"male": 16076, "female":  3772},
    "ASD_Unaffected": {"male":  3942, "female":  3882},
    "DDD":           {"male": 17422, "female": 13636},
}


# ---------------------------------------------------------------------------
# Rate / chromosome helpers
# ---------------------------------------------------------------------------

def compute_x_factor(male_n: int, female_n: int, alpha: float = 3.4) -> float:
    """Compute the X-chromosome mutation rate scaling factor.

    De novo mutations on the X chromosome arise at different rates in male and
    female transmissions, and the sex ratio within a cohort therefore affects
    the expected mutation count. This function derives a scaling factor relative
    to the autosomal rate so that a single ``cohort_prob`` column can be used
    for both autosomal and X-linked genes.

    Parameters
    ----------
    male_n:
        Number of male probands in the cohort.
    female_n:
        Number of female probands in the cohort.
    alpha:
        Male-to-female mutation rate ratio (default 3.4, from Kong et al. 2012).

    Returns
    -------
    float
        Multiplicative factor applied to X-linked per-site mutation probabilities.
    """
    autosomal = 2 * (male_n + female_n)
    # Mothers transmit to all children; fathers transmit only to daughters
    female_transmit = male_n + female_n
    male_transmit = female_n

    male_k = 2 / (1 + (1 / alpha))
    female_k = 2 / (1 + alpha)
    return ((male_transmit * male_k) + (female_transmit * female_k)) / autosomal


# ---------------------------------------------------------------------------
# Score normalization
# ---------------------------------------------------------------------------

def normalize_scores(weights: np.ndarray, all_scores: np.ndarray, norm: float,
                     rounding: int, pos_path: bool) -> np.ndarray:
    """Shift and truncate variant effect scores relative to a proteome percentile.

    Scores below the ``norm``-th percentile of the proteome-wide distribution
    are set to zero, effectively treating them as neutral. The remaining scores
    are rounded to ``rounding`` decimal places so that the simulation can group
    variants with identical weights, reducing memory use.

    Parameters
    ----------
    weights:
        Array of raw model scores for the variants/sites being normalized.
    all_scores:
        Proteome-wide array of raw model scores (used to compute the percentile).
    norm:
        Percentile threshold (0–100). Scores below this percentile become 0.
    rounding:
        Number of decimal places to round normalized scores to.
    pos_path:
        If True, higher scores indicate greater pathogenicity and are kept as-is.
        If False, scores are multiplied by -1 before normalization so that the
        most negative raw scores (most pathogenic) become the largest weights.

    Returns
    -------
    np.ndarray
        Normalized, non-negative, rounded scores.
    """
    if not pos_path:
        weights = weights * -1
    weights = weights - np.percentile(all_scores, norm)
    weights[weights < 0] = 0
    return weights.round(rounding)


# ---------------------------------------------------------------------------
# P-value computation (Monte Carlo Poisson simulation)
# ---------------------------------------------------------------------------

def get_prob(s_obs: float, n_exp: float, scores: np.ndarray,
             probs: np.ndarray, n_samples: int) -> float:
    """P-value using only the within-protein score distribution.

    Estimates P(S >= s_obs) by sampling mutation counts from a Poisson(n_exp)
    distribution and, for each count n, drawing n scores from the protein's
    own score distribution.

    Parameters
    ----------
    s_obs:
        Observed sum of normalized scores for variants seen in this protein.
    n_exp:
        Expected number of de novo mutations under the null model.
    scores:
        Unique normalized score values for sites within this protein.
    probs:
        Probability weights for each score value (sums to 1).
    n_samples:
        Number of Monte Carlo draws per Poisson count.

    Returns
    -------
    float
        Estimated one-sided p-value.
    """
    p0 = stats.poisson.pmf(0, n_exp) if s_obs == 0 else 0.0
    p1 = stats.poisson.pmf(1, n_exp)

    p1_weight_prob = (np.random.choice(scores, (1, n_samples), p=probs).sum(axis=0) >= s_obs).sum() / n_samples
    tot_prob = p0 + p1 * p1_weight_prob

    for n in range(2, 250):
        n_prob = stats.poisson.pmf(n, n_exp)
        weight_prob = (np.random.choice(scores, (n, n_samples), p=probs).sum(axis=0) >= s_obs).sum() / n_samples
        tot_prob += weight_prob * n_prob

        # Stop once the Poisson tail probability is negligible
        if 1 - stats.poisson.cdf(n, n_exp) < 1e-12 and n > n_exp:
            break

    return tot_prob


def get_all_prob(s_obs: float, n_exp: float, scores: np.ndarray, probs: np.ndarray,
                 all_probs: np.ndarray, all_scores: np.ndarray, n_samples: int) -> float:
    """P-value combining within-protein and proteome-wide score distributions.

    Both the within-protein score distribution and the proteome-wide distribution
    must independently produce a sum >= s_obs. This joint test is more conservative
    and guards against artefacts from the within-protein distribution alone.

    Parameters
    ----------
    s_obs, n_exp, scores, probs:
        As in :func:`get_prob`.
    all_probs:
        Proteome-wide probability weights for each score in ``all_scores``.
    all_scores:
        Proteome-wide unique normalized score values.
    n_samples:
        Number of Monte Carlo draws per Poisson count.

    Returns
    -------
    float
        Estimated one-sided p-value.
    """
    p0 = stats.poisson.pmf(0, n_exp) if s_obs == 0 else 0.0
    p1 = stats.poisson.pmf(1, n_exp)

    p1_weight_prob     = (np.random.choice(scores,     (1, n_samples), p=probs).sum(axis=0)     >= s_obs).sum() / n_samples
    p1_all_weight_prob = (np.random.choice(all_scores, (1, n_samples), p=all_probs).sum(axis=0) >= s_obs).sum() / n_samples
    tot_prob = p0 + p1 * p1_weight_prob * p1_all_weight_prob

    for n in range(2, 250):
        n_prob        = stats.poisson.pmf(n, n_exp)
        weight_prob     = (np.random.choice(scores,     (n, n_samples), p=probs).sum(axis=0)     >= s_obs).sum() / n_samples
        all_weight_prob = (np.random.choice(all_scores, (n, n_samples), p=all_probs).sum(axis=0) >= s_obs).sum() / n_samples
        tot_prob += weight_prob * all_weight_prob * n_prob

        if 1 - stats.poisson.cdf(n, n_exp) < 1e-12 and n > n_exp:
            break

    return tot_prob


def get_only_all_prob(s_obs: float, n_exp: float, scores: np.ndarray, probs: np.ndarray,
                      all_probs: np.ndarray, all_scores: np.ndarray, n_samples: int) -> float:
    """P-value using only the proteome-wide score distribution.

    Ignores the within-protein score distribution entirely; mutations are drawn
    from the proteome-wide distribution. Useful as a complementary test when the
    per-protein score distribution is sparse.

    Parameters
    ----------
    s_obs, n_exp, scores, probs:
        As in :func:`get_prob` (``scores``/``probs`` are unused here but kept
        for a consistent call signature).
    all_probs, all_scores, n_samples:
        As in :func:`get_all_prob`.

    Returns
    -------
    float
        Estimated one-sided p-value.
    """
    p0 = stats.poisson.pmf(0, n_exp) if s_obs == 0 else 0.0
    p1 = stats.poisson.pmf(1, n_exp)

    p1_all_weight_prob = (np.random.choice(all_scores, (1, n_samples), p=all_probs).sum(axis=0) >= s_obs).sum() / n_samples
    tot_prob = p0 + p1 * p1_all_weight_prob

    for n in range(2, 250):
        n_prob          = stats.poisson.pmf(n, n_exp)
        all_weight_prob = (np.random.choice(all_scores, (n, n_samples), p=all_probs).sum(axis=0) >= s_obs).sum() / n_samples
        tot_prob += all_weight_prob * n_prob

        if 1 - stats.poisson.cdf(n, n_exp) < 1e-12 and n > n_exp:
            break

    return tot_prob
  
  
# ---------------------------------------------------------------------------
# Per-protein test worker
# ---------------------------------------------------------------------------
 
# Module-level globals populated by the pool initializer so that the worker
# function is a plain top-level function (required for pickling by multiprocessing).
_worker_state: dict = {}
  
def _init_worker(df_variants, model, autosomal, x_factor,
                 all_probs, norm_all_scores, all_scores,
                 norm, rounding, n_samples, pos_path):
    """Pool initializer: store shared state in a module-level dict.
 
    Each worker process calls this once at startup. Storing data here avoids
    pickling the large DataFrames and arrays on every task, and allows
    """
    _worker_state["df_variants"]     = df_variants
    _worker_state["model"]           = model
    _worker_state["autosomal"]       = autosomal
    _worker_state["x_factor"]        = x_factor
    _worker_state["all_probs"]       = all_probs
    _worker_state["norm_all_scores"] = norm_all_scores
    _worker_state["all_scores"]      = all_scores
    _worker_state["norm"]            = norm
    _worker_state["rounding"]        = rounding
    _worker_state["n_samples"]       = n_samples
    _worker_state["pos_path"]        = pos_path
 
 
def test_protein(params):
    """Test a single protein for de novo variant score enrichment.
 
    This is the top-level worker function dispatched to each pool process.
    All shared data is read from the module-level ``_worker_state`` dict
    populated by :func:`_init_worker`.
 
    Parameters
    ----------
    params:
        A row from ``df_genes.values`` with elements
        [protein, gene, chrom, rates_filepath].
 
    Returns
    -------
    list
        ``[protein, s_obs, N, n_exp, s_exp, pvalue, pvalue_proteome,
        pvalue_only_proteome]``, or an empty list if no score data is available
        for this protein.
    """
    protein, gene, chrom, rates_filepath = params
 
    # Unpack shared state set by the pool initializer
    df_variants     = _worker_state["df_variants"]
    model           = _worker_state["model"]
    autosomal       = _worker_state["autosomal"]
    x_factor        = _worker_state["x_factor"]
    all_probs       = _worker_state["all_probs"]
    norm_all_scores = _worker_state["norm_all_scores"]
    all_scores      = _worker_state["all_scores"]
    norm            = _worker_state["norm"]
    rounding        = _worker_state["rounding"]
    n_samples       = _worker_state["n_samples"]
    pos_path        = _worker_state["pos_path"]
 
    # Observed variants for this protein
    df_mut = df_variants[df_variants.protein == protein]
 
    # Load per-protein mutation rate table
    rates = pd.read_csv(rates_filepath)
    rates.prob = rates.prob.astype(float)
 
    # Scale per-site probabilities to the whole cohort; X-linked genes
    # require an additional sex-ratio correction.
    if chrom == "X":
        rates["cohort_prob"] = rates.prob * autosomal * x_factor
    else:
        rates["cohort_prob"] = rates.prob * autosomal
 
    # Normalize model scores relative to the proteome-wide distribution
    rates["weight"] = normalize_scores(rates[model].to_numpy(), all_scores, norm, rounding, pos_path)
    rates = rates.groupby(["mutant", "weight"]).sum().reset_index()
 
    # Match observed variants to their normalized weights
    df_mut = pd.merge(df_mut[["mutant"]], rates[["mutant", "weight"]])
    df_mut = df_mut.dropna(subset="weight")
 
    # Summary statistics
    n_exp = rates.cohort_prob.sum()
    s_exp = np.sum(rates.cohort_prob * rates.weight)
    s_obs = df_mut.weight.sum()
    N     = df_mut.shape[0]
 
    # Collapse to unique (score, total_prob) pairs for efficient simulation
    probs  = rates.groupby("weight").cohort_prob.sum()
    scores = probs.index.astype(float).to_numpy()
    if len(scores) == 0:
        return []
    probs = (probs / probs.sum()).astype(float).to_numpy()
 
    return [
        protein, s_obs, N, n_exp, s_exp,
        get_prob(s_obs, n_exp, scores, probs, n_samples),
        get_all_prob(s_obs, n_exp, scores, probs, all_probs, norm_all_scores, n_samples),
        get_only_all_prob(s_obs, n_exp, scores, probs, all_probs, norm_all_scores, n_samples),
    ]
 

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
 
def main():
    args = parse_args()
 
    # ------------------------------------------------------------------
    # Load inputs
    # ------------------------------------------------------------------
    df_genes    = pd.read_csv(args.genes,     sep="\t")
    df_all      = pd.read_csv(args.all_rates, sep="\t")
    df_variants = pd.read_csv(args.variants,  sep="\t")
 
    # ------------------------------------------------------------------
    # Cohort-level mutation rate setup
    # ------------------------------------------------------------------
    male_n   = COHORT_SEX_SIZES[args.study]["male"]
    female_n = COHORT_SEX_SIZES[args.study]["female"]
    autosomal = 2 * (male_n + female_n)
    x_factor  = compute_x_factor(male_n, female_n)
 
    # ------------------------------------------------------------------
    # Filter variants to the chosen study and model
    # ------------------------------------------------------------------
    df_variants = df_variants[df_variants.study.isin([args.study])]
    df_variants = df_variants[~df_variants[args.model].isna()]
    df_variants = df_variants.drop_duplicates(subset="variant_id")
 
    # ------------------------------------------------------------------
    # Proteome-wide score normalization (used in all three p-value tests)
    # ------------------------------------------------------------------
    all_scores      = df_all[args.model].astype(float).to_numpy()
    norm_all_scores = normalize_scores(all_scores, all_scores, args.norm, args.rounding, args.pos_path)
    all_probs       = df_all.prob
    all_probs       = (all_probs / all_probs.sum()).astype(float).to_numpy()
 
    # ------------------------------------------------------------------
    # Build per-protein mapping for gene annotation
    # ------------------------------------------------------------------
    protein_to_gene = df_genes.set_index("protein")["gene"].to_dict()
 
    # ------------------------------------------------------------------
    # Run per-protein tests in parallel
    # ------------------------------------------------------------------
    init_args = (
        df_variants, args.model, autosomal, x_factor,
        all_probs, norm_all_scores, all_scores,
        args.norm, args.rounding, args.n_samples, args.pos_path,
    )
 
    with mp.Pool(args.workers, initializer=_init_worker, initargs=init_args) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(test_protein, df_genes.values),
                total=df_genes.shape[0],
                desc="Testing proteins",
            )
        )
 
    # ------------------------------------------------------------------
    # Assemble and save results
    # ------------------------------------------------------------------
    df_results = pd.DataFrame(
        [row for row in results if len(row) > 0],
        columns=["protein", "s_obs", "n", "n_exp", "s_exp",
                 "pvalue", "pvalue_proteome", "pvalue_only_proteome"],
    )
    df_results["gene"] = df_results.protein.map(protein_to_gene)
 
    outdir = args.outdir.rstrip("/") + "/"
    out_path = (
        f"{outdir}{args.study}_{args.model}_"
        f"{args.norm}_{args.rounding}_{args.n_samples}.csv"
    )
    df_results.to_csv(out_path,index=False)
    print(f"Results written to {out_path}")
 
 
if __name__ == "__main__":
    main()