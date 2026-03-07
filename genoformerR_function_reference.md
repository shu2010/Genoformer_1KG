# genoformerR — Function Reference (v0.3.1)

## Overview

genoformerR implements a transformer-based, ancestry-aware polygenic risk scoring
pipeline. Functions are grouped into seven stages that can be run individually or
via the end-to-end `gf_run_pipeline()` wrapper.

**v0.3.1 changes** — Both `gf_ancestry()` and `gf_windowed_pca()` now support
reference-based projection: PCA is fitted on a reference cohort (e.g. 1KGP3) and
study or new-individual samples are projected into that fixed subspace. This is
required for consistent ancestry inference when scoring individuals outside the
training cohort. Three new functions support windowed PCA persistence:
`gf_save_window_pcas()`, `gf_load_window_pcas()`.

---

## Stage 0 — Initialisation

### `gf_init(envname, python, reinstall)`
Creates (or reuses) a conda environment named `"genoformer"`, installs PyTorch +
scikit-learn + numpy, and sources `model.py` / `train.py` into the Python session.
Must be called once per R session before any modelling function.

```r
gf_init()                        # uses "genoformer" conda env
gf_init(reinstall = TRUE)        # force-reinstall Python deps
```

---

## Stage 1 — Data acquisition & QC

### `gf_download_1kg(outdir, chroms, panel_only, ncores)`
Downloads per-chromosome VCFs and the sample panel file from the EBI 1KGP3 FTP.
Parallelised across chromosomes. Returns paths of downloaded files.

The **panel file** (`integrated_call_samples_v3.20130502.ALL.panel`) is a plain
tab-separated file with four columns: `sample`, `pop`, `super_pop`, `gender`.
It contains no genotype data — only population metadata. It is the required input
to `gf_ancestry(panel_file=)`.

### `gf_qc(vcf_dir, pgs_file, outdir, maf, hwe, geno_miss, mind_miss, ld_window, ld_step, ld_r2, plink2, chroms)`
Full PLINK2 genotype QC pipeline:
1. VCF to pgen per chromosome
2. Merge chromosomes
3. Apply MAF / HWE / missingness filters
4. LD prune (for ancestry PCA)
5. Intersect with PGS scoring file variants

Returns a list with `pgen_prefix`, `pruned_prefix`, `dosage_raw`, `n_samples`,
`n_snps_qc`, `n_snps_pgs`.

---

## Stage 2 — Global ancestry inference

### `gf_ancestry(dosage_mat, panel_file, ref_dosage_raw, n_pcs, method, n_trees, seed)`

Infers genome-wide superpopulation ancestry. Operates in one of two modes
determined by whether `ref_dosage_raw` is provided.

**Mode A — reference projection** (`ref_dosage_raw` provided, recommended for
external cohorts): PCA is fitted on the 1KGP3 reference samples in `ref_dosage_raw`.
Study samples in `dosage_mat` are then projected into that fixed PC space using
`stats::predict.prcomp`. The Random Forest classifier is trained on the reference
labels from `panel_file` and applied to the projected study coordinates. Ensures
ancestry proportions are comparable across independent scoring runs.

**Mode B — study-only** (`ref_dosage_raw = NULL`): PCA is fitted on `dosage_mat`
directly. Only valid when `dosage_mat` already contains 1KGP3 samples whose IDs
appear in `panel_file`. Errors with a helpful message if no sample IDs match.

```r
# Mode A — external study cohort (general case)
anc <- gf_ancestry(
  dosage_mat     = study_dosage,
  panel_file     = "data/1kg/integrated_call_samples_v3.20130502.ALL.panel",
  ref_dosage_raw = "data/plink/kg3_pruned.raw"
)

# Mode B — dosage_mat IS 1KGP3
anc <- gf_ancestry(
  dosage_mat = kg3_dosage,
  panel_file = "data/1kg/integrated_call_samples_v3.20130502.ALL.panel"
)
```

Returns a `gf_ancestry` S3 object:

| Field | Description |
|---|---|
| `pcs` | (N x n_pcs) PC matrix for study samples |
| `ref_pcs` | (N_ref x n_pcs) reference PC matrix (Mode A only) |
| `anc_probs` | (N x 5) soft ancestry proportions (AFR, AMR, EAS, EUR, SAS) |
| `labels` | MAP superpopulation per sample |
| `entropy` | Admixture entropy per sample |
| `pca_obj` | Fitted `prcomp` / sklearn PCA — save for projection at scoring time |
| `clf_obj` | Trained `ranger` / sklearn GBM classifier — save alongside `pca_obj` |
| `sample_ids` | Sample IDs from dosage matrix |

Save the full object as `ancestry.rds` after training. At scoring time load it and
call `stats::predict(anc$pca_obj, newdata=new_dosage)` to project new individuals,
then `predict(anc$clf_obj, ...)` to obtain probabilities.

### `print(anc)` / `plot(anc, pcs)`
Print summary table or PCA scatter coloured by inferred superpopulation.

---

## Stage 3 — Local ancestry inference

### `gf_local_ancestry(method, ...)` — unified dispatcher
Single entry point. Returns a **(N x M x K) array** — one K-vector of ancestry
posteriors per sample per SNP. All downstream functions accept this array directly.

| `method` | Tool | Phasing | Reference panel | Consistent across runs |
|---|---|---|---|---|
| `"flare"` | FLARE Java jar | No | Yes — 1KGP3 VCF | Yes — reference anchors K axes |
| `"rfmix2"` | rfmix2 binary | Yes (SHAPEIT4) | Yes | Yes |
| `"windowed_pca"` | pure R | No | Optional | Only with `ref_window_pcas=` |

### `gf_local_ancestry_flare(vcf_file, ref_vcf, sample_map, gmap_dir, outdir, chroms, k_pops, flare_jar, java, min_maf, nthreads, em_its)`

Runs FLARE (Hidden Markov Model) per chromosome. FLARE uses the reference panel
in two ways: (1) emission probabilities from per-population allele frequencies at
each SNP position; (2) transition probabilities from the genetic map recombination
rates. The `ref_vcf` (1KGP3 joint-called VCF) provides allele frequencies.
The `sample_map` defines K populations and maps reference VCF samples to them.

Accepts unphased study VCF. Because FLARE's K dimensions are anchored to the named
reference populations, the coordinate system is automatically consistent between
training and scoring runs — no equivalent of `gf_save_window_pcas` is needed.

Returns named list of output prefixes. Pass to `gf_parse_flare()`.

`sample_map` format (2-column headerless TSV):
```
HG00096    EUR
NA19701    AMR
HG01880    AFR
```

### `gf_parse_flare(flare_prefixes, snp_ids, snp_positions, snp_chroms, sample_ids, global_anc, entropy_threshold, n_cores)`

Parses FLARE `.anc.vcf.gz` files. The `ANS` FORMAT field gives two population
integers per sample (one per haplotype); each contributes 0.5 weight to its
population, producing diploid proportions. Non-admixed samples (entropy below
`entropy_threshold`) are filled with global ancestry as fallback. Population
names are read from the `##ancestry=` VCF header line.
Returns **(N x M x K) array**.

### `gf_windowed_pca(dosage_mat, snp_chroms, ref_dosage_mat, ref_window_pcas, k_pops, window_size, step_size, min_var, scale_snps, n_cores)`

Pure R sliding-window PCA. Within each genomic window, fits PCA and
softmax-normalises scores to produce pseudo-ancestry proportions summing to 1.
Three modes:

**Mode A** (`ref_dosage_mat` provided): PCA fitted on reference per window; study
samples projected using saved rotation. Recommended for external study cohorts.

**Mode B** (neither argument): PCA fitted on `dosage_mat`. Valid for initial 1KGP3
training runs. Issues warning if N < 30.

**Mode C** (`ref_window_pcas` provided): projects `dosage_mat` into saved per-window
subspaces from a prior training run. Required for scoring new individuals.

The returned array carries a `"window_pcas"` attribute containing fitted `prcomp`
objects keyed by `"chr:window_start"`. Always save these after training.

```r
# Training on 1KGP3 — Mode B
la <- gf_windowed_pca(kg3_dosage, snp_chroms, k_pops=5,
                       window_size=500, step_size=250)
gf_save_window_pcas(la, "models/window_pcas.rds")

# Scoring new individual — Mode C
ref_pcas <- gf_load_window_pcas("models/window_pcas.rds")
la_new   <- gf_windowed_pca(new_dosage, snp_chroms, ref_window_pcas=ref_pcas)
```

### `gf_save_window_pcas(local_anc_arr, path)`
Extracts the `"window_pcas"` attribute and writes it to an RDS file. Must be called
after any training windowed PCA run. Without this file, new individuals cannot be
projected into the training subspace.

### `gf_load_window_pcas(path)`
Loads per-window `prcomp` objects from RDS. Pass as `ref_window_pcas=` to
`gf_windowed_pca()` when scoring new individuals.

### `gf_phase(pgen_prefix, outdir, gmap_dir, chroms, shapeit4, scaffold, ncores, pbwt_depth)`
SHAPEIT4 wrapper — required pre-processing for `method = "rfmix2"` only.
Not needed for FLARE or windowed PCA.

### `gf_local_ancestry_rfmix2(vcf_file, ref_vcf, sample_map, gmap_dir, outdir, chroms, k_pops, rfmix2, em_iterations)`
Legacy RFMix2 wrapper. Requires phased VCFs. Pass output to `gf_align_rfmix2()`.

### `gf_align_rfmix2(rfmix_prefixes, snp_ids, snp_positions, snp_chroms, sample_ids, global_anc, entropy_threshold, n_cores)`
Parses RFMix2 `.fb.tsv` forward-backward posteriors. Nearest-SNP interpolation
aligns output to dosage SNP order. Returns **(N x M x K) array**.

### `gf_align_local_anc(raw_prefixes, ...)`
Auto-detects FLARE vs RFMix2 output by file extension and dispatches accordingly.

### `gf_local_anc_summary(local_anc_arr, snp_chroms, pop_labels)`
Diagnostic summary: mean ancestry per labelled population per chromosome.

---

## Stage 4 — Model construction

### `gf_build_model(n_snps, n_ld_blocks, d_model, n_heads, n_layers, anc_dim, k_pops, d_local, use_local_anc, dropout, device)`

Constructs the GenoFormer PyTorch module via reticulate. Takes only integer scalars
derived from the data — no file data enters this function.

| Parameter | Default | Derived from |
|---|---|---|
| `n_snps` | 80000 | `ncol(dosage_mat)` |
| `n_ld_blocks` | 1700 | `max(ld_blocks) + 1` |
| `d_model` | 256 | — |
| `n_heads` | 8 | — |
| `n_layers` | 6 | — |
| `anc_dim` | 5 | `ncol(anc$anc_probs)` |
| `k_pops` | 5 | `dim(local_anc)[3]` |
| `d_local` | 64 | — |
| `use_local_anc` | TRUE | — |

Returns a `gf_model` S3 object with `model_py`, `device`, `n_params`, `config`.

### `print(model)`
Prints device, parameter count, and full architecture config.

---

## Stage 5 — Training

### `gf_train(model, dosage_mat, anc_probs, pgs_weights, ld_blocks, chroms, phenotype, pop_labels, local_anc, epochs, batch_size, lr, weight_decay, anc_weight, calib_weight, val_split, save_path, verbose)`

Multi-task training via AdamW + cosine LR schedule. Three loss terms:
1. PRS MSE on z-scored phenotype
2. Population CE on auxiliary head x `anc_weight`
3. Ancestry calibration loss x `calib_weight` (penalises population-stratified residuals)

Input tensor mapping:

| R argument | Shape | Where used in model |
|---|---|---|
| `dosage_mat` | (B, M) float32 | SNP token — `dosage_proj` |
| `pgs_weights` → `log_beta` | (B, M) float32 | SNP token — `beta_proj` |
| `ld_blocks` | (B, M) int64 | SNP token — `ld_emb` |
| `chroms` | (B, M) int64 | SNP token — `chr_emb` |
| `anc_probs` | (B, 5) float32 | FiLM on every transformer block |
| `local_anc` | (B, M, K) float32 | SNP token fusion layer |
| `phenotype` | (B,) float32 | MSE + calibration loss |
| `pop_labels` | (B,) int64 | CE auxiliary loss |

`local_anc` is optional; if NULL, model trains in global-only mode.
Returns `list(model, history)`.

---

## Stage 6 — Inference

### `gf_predict(model, dosage_mat, anc_probs, pgs_weights, ld_blocks, chroms, effect_alleles, dosage_alleles, local_anc, batch_size)`

Batched transformer inference. Returns a `gf_prs` tibble:

| Column | Description |
|---|---|
| `transformer_prs` | z-scored GenoFormer PRS (primary output) |
| `classical_prs` | z-scored classical weighted-dosage sum |
| `local_anc_available` | TRUE if local_anc was passed |
| `AFR`, `AMR`, `EAS`, `EUR`, `SAS` | global soft ancestry proportions |

If `use_local_anc=TRUE` but `local_anc=NULL`, warns and broadcasts global ancestry
per-SNP as fallback. For admixed individuals, always provide `local_anc`.

Score standardisation note: both columns are z-scored within the current batch.
For percentile reporting use population means/SDs saved from training predictions.

---

## Stage 7 — Evaluation

### `gf_evaluate(prs_tbl, phenotype, pop_labels, populations, ref_pop, n_bootstrap)`

Cross-ancestry portability evaluation. Per population: Pearson R² with bootstrap
CIs, calibration slope, portability ratio (R²_pop / R²_ref_pop).
Returns `gf_eval` with `metrics`, `portability_plot`, `distribution_plot`.

---

## Model persistence

### `gf_save_model(model, path, save_config)`
Saves PyTorch state dict. `save_config=TRUE` writes `_config.json` alongside.

### `gf_load_model(path, config, device)`
Loads weights. If `config=NULL` auto-detects `_config.json` next to the `.pt` file.

---

## End-to-end pipeline

### `gf_run_pipeline(..., local_anc_method)`
Runs all stages in sequence. `local_anc_method`: `"flare"`, `"windowed_pca"`,
`"rfmix2"`, or `"none"`. Returns named list: `ancestry`, `local_anc_arr`, `model`,
`history`, `prs`, `evaluation`, `paths`.

---

## S3 classes

| Class | Created by | Key fields |
|---|---|---|
| `gf_ancestry` | `gf_ancestry()` | `pcs`, `ref_pcs`*, `anc_probs`, `labels`, `entropy`, `pca_obj`, `clf_obj` |
| `gf_model` | `gf_build_model()` | `model_py`, `device`, `n_params`, `config` |
| `gf_prs` | `gf_predict()` | `transformer_prs`, `classical_prs`, `local_anc_available`, ancestry cols |
| `gf_eval` | `gf_evaluate()` | `metrics`, `portability_plot`, `distribution_plot` |

\* `ref_pcs` present only in Mode A.

---

## Files to preserve after training

| File | Written by | Purpose at scoring time |
|---|---|---|
| `genoformer_best.pt` | `gf_save_model()` | Model weights |
| `genoformer_best_config.json` | `gf_save_model(save_config=TRUE)` | Architecture reconstruction |
| `ancestry.rds` | `saveRDS(anc, ...)` | PCA + RF for new-individual global ancestry |
| `window_pcas.rds` | `gf_save_window_pcas()` | Windowed PCA projection (if used) |
| `pgs_snps.txt` | `gf_qc()` | SNP extraction from new VCF |
| `train_prs_scores.rds` | `saveRDS(prs_train, ...)` | Population means/SDs for percentile reporting |

---

## External dependencies

| Tool | Used by | Required for |
|---|---|---|
| PLINK2 | `gf_qc()` | Genotype QC and dosage export |
| conda / Miniconda | `gf_init()` | Python backend |
| PyTorch | `gf_build_model()`, `gf_train()`, `gf_predict()` | Model training & inference |
| FLARE (`flare.jar`, Java >= 11) | `gf_local_ancestry_flare()` | Local ancestry (recommended) |
| SHAPEIT4 | `gf_phase()` | Pre-phasing for RFMix2 only |
| RFMix2 | `gf_local_ancestry_rfmix2()` | Local ancestry (legacy) |
| ranger | `gf_ancestry(method="r")` | Global ancestry R-only path |
