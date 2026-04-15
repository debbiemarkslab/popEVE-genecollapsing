# popEVE Gene Collapsing

Score and rank genes for de novo variant enrichment using weighted mutation rate models.

For each protein-coding gene, `gene_collapse.py` computes a p-value reflecting whether the observed de novo variants in a cohort carry unexpectedly high pathogenicity scores given the gene's background mutation rate. Scores are drawn from a user-specified variant effect model (e.g. popEVE), normalized relative to the 90th percentile of the proteome-wide distribution, and used as weights in a Poisson simulation. Three complementary tests are run:

| Test | Description |
|---|---|
| `pvalue` | Score enrichment tested within the protein only |
| `pvalue_proteome` | Score enrichment tested both within the protein and against the proteome-wide score distribution |
| `pvalue_only_proteome` | Score enrichment tested against the proteome-wide score distribution only |

`build_all_rates.py` is a helper script that generates the proteome-wide rates file required as input to `gene_collapse.py`.

---

## Scripts

### `build_all_rates.py`

Aggregates per-protein mutation rate files into a single proteome-wide rates table. Run this first to produce `all_rates.txt`.

**Usage**

```bash
python build_all_rates.py \
    --rates-dir /path/to/rates/ \
    --model     popEVE \
    --rounding  4 \
    --workers   40 \
    --outfile   all_rates.txt
```

**Arguments**

| Argument | Required | Default | Description |
|---|---|---|---|
| `--rates-dir` | Yes | — | Directory containing per-protein mutation rate files |
| `--model` | No | `popEVE` | Column name of the variant effect model scores to aggregate |
| `--rounding` | No | `4` | Decimal places to round model scores to before summing probabilities. Should match the value used in `gene_collapse.py` |
| `--workers` | No | `40` | Number of parallel worker processes |
| `--outfile` | No | `all_rates.txt` | Path for the output tab-separated rates file |

**Input: per-protein rate files** (`--rates-dir`)

One file per protein, in any separator format readable by pandas. Each file must contain the following columns:

| Column | Description |
|---|---|
| `Consequence` | Variant consequence annotation — only rows containing `"missense"` are retained |
| `prob` | Per-site mutation probability |
| `mutant` | Mutation identifier |
| `<model>` | Pathogenicity score for the chosen model (e.g. `popEVE`) |

**Output: `all_rates.txt`**

Tab-separated file with one row per unique rounded score value across the proteome.

| Column | Description |
|---|---|
| `<model>` | Rounded pathogenicity score (e.g. `popEVE`) |
| `prob` | Summed mutation probability across the proteome for that score |

---

### `gene_collapse.py`

Runs the per-gene enrichment tests. Requires `all_rates.txt` produced by `build_all_rates.py`.

**Usage**

```bash
# Minimal
python gene_collapse.py \
    --genes     genes_info.txt \
    --all-rates all_rates.txt \
    --variants  ddd_variants.txt \
    --outdir    results/ \
    --study     DDD

# All arguments
python gene_collapse.py \
    --genes     genes_info.txt \
    --all-rates all_rates.txt \
    --variants  ddd_variants.txt \
    --outdir    results/ \
    --study     DDD \
    --model     popEVE \
    --norm      75 \
    --pos-path \
    --rounding  4 \
    --n-samples 10000 \
    --workers   4
```

**Arguments**

| Argument | Required | Default | Description |
|---|---|---|---|
| `--genes` | Yes | — | Tab-separated gene/protein info file (e.g. `genes_info.txt`) |
| `--all-rates` | Yes | — | Tab-separated proteome-wide mutation rates file (e.g. `all_rates.txt`) |
| `--variants` | Yes | — | Tab-separated observed de novo variants file (e.g. `ddd_variants.txt`) |
| `--outdir` | Yes | — | Directory where the output pickle file will be written |
| `--study` | Yes | — | Cohort to analyse (`ASD_Affected`, `ASD_Unaffected`, or `DDD`). Controls which cohort size constants are used |
| `--model` | No | `popEVE` | Column name of the variant effect model scores to use |
| `--norm` | No | `75.0` | Percentile used for score normalization — scores below this percentile are set to zero |
| `--pos-path` | No | `False` | If set, scores are kept as-is (higher = more pathogenic). If unset, scores are multiplied by −1 before normalization so that the most negative raw scores become the largest weights |
| `--rounding` | No | `4` | Decimal places to round normalized scores to. Finer rounding increases precision but slows simulation |
| `--n-samples` | No | `10000` | Number of Monte Carlo samples used to estimate each p-value. Higher values give more precise p-values at greater runtime cost |
| `--workers` | No | `4` | Number of parallel worker processes |

**Input: genes file** (`--genes`)

Tab-separated, one row per protein.

| Column | Description |
|---|---|
| `protein` | Protein identifier — must match values in the variants file |
| `gene` | HGNC gene symbol |
| `chrom` | Chromosome — use `X` for X-linked genes |
| `rates_filepath` | Path to the per-protein rates `.pkl` file |

**Input: all_rates file** (`--all-rates`)

Tab-separated proteome-wide mutation rates, as produced by `build_all_rates.py`.

| Column | Description |
|---|---|
| `<model>` | Pathogenicity score for the chosen model (e.g. `popEVE`) |
| `prob` | Per-site mutation probability summed across the proteome for that score |

**Input: variants file** (`--variants`)

Tab-separated observed de novo variants.

| Column | Description |
|---|---|
| `study` | Cohort label (e.g. `DDD`, `ASD_Affected`) |
| `protein` | Protein identifier — must match values in the genes file |
| `mutant` | Mutation identifier — must match values in the per-protein rates files |
| `variant_id` | Unique variant identifier including subject ID, used for deduplication |
| `<model>` | Pathogenicity score for the chosen model (e.g. `popEVE`) |

**Output**

A pickled pandas DataFrame saved to:

```
<outdir>/<study>_<model>_<norm>_<rounding>_<n_samples>.pkl
```

| Column | Description |
|---|---|
| `protein` | Protein identifier |
| `gene` | HGNC gene symbol |
| `s_obs` | Observed sum of normalized pathogenicity scores for variants in this protein |
| `n` | Observed number of de novo variants |
| `n_exp` | Expected number of de novo variants under the null mutation rate model |
| `s_exp` | Expected score sum under the null model |
| `pvalue` | P-value from within-protein score enrichment test |
| `pvalue_proteome` | P-value from joint within-protein and proteome-wide score enrichment test |
| `pvalue_only_proteome` | P-value from proteome-wide score enrichment test only |
