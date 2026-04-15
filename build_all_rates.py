"""
build_all_rates.py
------------------
Aggregate per-protein mutation rate files into a single proteome-wide rates
table (all_rates.txt) for use as input to gene_collapse.py.

For each per-protein rates file, missense variants are extracted and their
mutation probabilities are summed by (rounded) model score. The results are
then concatenated across all proteins and summed again, yielding a single
table of proteome-wide score–probability pairs.

Usage
-----
    python build_all_rates.py --rates-dir /path/to/rates/ --model popEVE

    python build_all_rates.py \\
        --rates-dir /n/groups/marks/projects/Poly/users/rose/data/rates_v4/ \\
        --model     popEVE \\
        --rounding  4 \\
        --workers   40 \\
        --outfile   all_rates.txt

Input
-----
rates-dir: Directory containing per-protein rate files. Each file must be a
           CSV (any separator pandas can infer) with at least these columns:
               Consequence -- variant consequence annotation (e.g. 'missense_variant')
               prob        -- per-site mutation probability
               mutant      -- mutation identifier
               <model>     -- pathogenicity score for the chosen model

Output
------
A tab-separated file with columns:
    <model>   -- rounded pathogenicity score
    prob      -- summed mutation probability across the proteome for that score
"""

import argparse
import multiprocessing as mp
from glob import glob

import pandas as pd
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build the proteome-wide all_rates.txt input for gene_collapse.py.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--rates-dir", required=True, metavar="DIR",
        help="Directory containing per-protein mutation rate files.",
    )
    parser.add_argument(
        "--model", default="popEVE", metavar="STR",
        help="Column name of the variant effect model scores to aggregate.",
    )
    parser.add_argument(
        "--rounding", type=int, default=4, metavar="INT",
        help="Number of decimal places to round model scores to before summing "
             "probabilities. Should match the value used in gene_collapse.py.",
    )
    parser.add_argument(
        "--workers", type=int, default=40, metavar="INT",
        help="Number of parallel worker processes.",
    )
    parser.add_argument(
        "--outfile", default="all_rates.txt", metavar="FILE",
        help="Path for the output tab-separated rates file.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Per-file worker
# ---------------------------------------------------------------------------

# Module-level globals populated by the pool initializer (required for
# multiprocessing picklability).
_model: str = ""
_rounding: int = 4


def _init_worker(model: str, rounding: int) -> None:
    """Pool initializer: store model name and rounding in each worker process."""
    global _model, _rounding
    _model    = model
    _rounding = rounding


def _aggregate_file(filepath: str) -> pd.DataFrame:
    """Load one per-protein rates file and return a score–probability table.

    Filters to missense variants, drops rows with missing model scores, rounds
    scores to the configured precision, then sums probabilities within each
    score bin.

    Parameters
    ----------
    filepath:
        Path to a single per-protein rates CSV file.

    Returns
    -------
    pd.DataFrame
        Two-column DataFrame with columns [<model>, 'prob'], one row per
        unique rounded score value present in this protein.
    """
    df = pd.read_csv(filepath, usecols=["Consequence", _model, "mutant", "prob"])
    df = df[df.Consequence.str.contains("missense")]
    df = df.dropna(subset=[_model])
    df[_model] = df[_model].round(_rounding)
    df = df.groupby(_model).prob.sum().reset_index()
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    filepaths = glob(args.rates_dir.rstrip("/") + "/*")
    if not filepaths:
        raise FileNotFoundError(f"No files found in rates directory: {args.rates_dir}")
    print(f"Found {len(filepaths)} rate files. Aggregating with {args.workers} workers...")

    with mp.Pool(args.workers, initializer=_init_worker, initargs=(args.model, args.rounding)) as pool:
        results = list(tqdm(pool.imap_unordered(_aggregate_file, filepaths), total=len(filepaths)))

    df_all = pd.concat(results, ignore_index=True)
    df_all = df_all.groupby(args.model).prob.sum().reset_index()

    df_all.to_csv(args.outfile, sep="\t", index=False)
    print(f"Proteome-wide rates written to {args.outfile} ({len(df_all)} score bins)")


if __name__ == "__main__":
    main()
