# SMC-ABC Perturbation Kernels on the Simplex — MATLAB Code

MATLAB source code accompanying:

> **Bayesian Parameter Inference on the Simplex: Convenient Transformations and Perturbation Kernels for Sequential Monte Carlo Sampling**
> Pedro A. Saa, Gabriel I. Miranda-Contreras, Christopher Drovandi, Lars K. Nielsen.

Two families are compared:

- **Dirichlet kernel (DK)** — perturbs particles directly on the simplex, tuned through a single concentration hyper-parameter `κ`.
- **Multivariate Gaussian kernel (MGK)** — applied *after* mapping the simplex onto free variables, using four transformations (LT, LGT, ILRT, BLT).

---

## Table of contents

- [Quick start](#quick-start)
- [Figure numbering](#figure-numbering)
- [Directory layout](#directory-layout)
- [Case studies](#case-studies)
- [Kernel names: code ↔ paper](#kernel-names-code--paper)
- [DK tuning methods](#dk-tuning-methods)
- [Common ASMC settings](#common-asmc-settings)
- [Shared helper functions](#shared-helper-functions)
- [Output: `population` structure](#output-population-structure)
- [Results files and naming convention](#results-files-and-naming-convention)
- [Requirements](#requirements)

---

## Quick start

Each case study lives in its own folder and is driven by a `run_*.m` function taking the **kernel index** as its only argument. `cd` into the folder first — the scripts use relative `addpath`, `load` and `save`.

Figure names below are the file names (`Fig4.png`, …); the printed figure numbers differ, see [Figure numbering](#figure-numbering).

```matlab
% Dirichlet kernel tuning study (Fig3.png)
cd DK_tuning
run_example_scripts                       % 4 models x 3 datasets x 3 tuning methods
make_figures_dk                           % metrics table + panels A-B
ABC_solution_comparison                   % panel C
make_figure_S1                            % Supplementary Fig. S1

% Multinomial benchmark with a known posterior (Fig4.png)
cd ../multinomial_posterior
for k = [1 3 4 6 7]; run_multinomial(k); end
make_figures_multinomial                  % metrics table + assembled figure

% Pair-HMM global alignment (Fig5.png)
cd ../hmm_global_alignment
for k = [1 3 4 6 7]; run_global_aligment_hmm(k); end
make_figures_hmm                          % metrics table + assembled figure

% COVID-19 screening Bayes net (Fig6.png)
cd '../COVID-19 model'
for k = [1 3 4 6 7]; run_bayes_net(k); end
make_figures_covid                        % metrics table + assembled figure
classifier_performance                    % Table 1 and Supplementary Table S6
```

Pre-computed results are included, so every analysis and figure script runs without re-running the (expensive) samplers. The drivers create `results/` if it is missing.

Two further scripts are diagnostics whose output is not an item of the paper: `make_network_figure` (in `COVID-19 model/`) draws the DAG of the Bayes net, and `getPosteriorStructure('covid' | 'hmm' | 'all')`, run from the *root* of `scripts/`, writes the posterior dependence structure to `analysis/`.

---

## Figure numbering

Every script, comment and output file in `scripts/` refers to a figure by its **file name** in `submission/figures_final/`. Because the Dirichlet-kernel illustration (`Fig2.png`) is placed in the Methods section, the printed numbers are shifted:

| File | Printed in the paper | Content | Produced by |
|---|---|---|---|
| `Fig1.pdf` | Fig. 1 | perturbation strategies (schematic) | not produced here |
| `Fig3.png` | **Fig. 2** | DK tuning strategies | `DK_tuning/` |
| `Fig4.png` | **Fig. 3** | multinomial benchmark, known posterior | `multinomial_posterior/` |
| `Fig5.png` | **Fig. 4** | pair-HMM global alignment | `hmm_global_alignment/` |
| `Fig6.png` | **Fig. 5** | COVID-19 screening Bayes net | `COVID-19 model/` |
| `Fig2.png` | **Fig. 6** | Dirichlet kernel concentration (Methods) | not produced here |

So "Fig. 4" in a comment means `Fig4.png`, printed as Fig. 3.

---

## Directory layout

```
scripts/
├── getPosteriorStructure.m     # diagnostic: Sigma in ILR coordinates, D_cross, variation matrices
├── analysis/                   # its output (tables + figures), written on demand
│
├── DK_tuning/                  # Fig3.png — multinomial models, DK tuning comparison
│   ├── run_example_scripts.m       <- DRIVER
│   ├── multinomial_example.m       <- ASMC sampler for the multinomial model
│   ├── getKernelParam.m            <- the 3 DK tuning methods
│   ├── ABC_solution_comparison.m   <- ABC vs. analytical posterior (panel C)
│   ├── make_figures_dk.m           <- metrics table + panels A-B
│   ├── make_figure_S1.m            <- exact ABC posterior at eps_T = 0 (Supp. Fig. S1)
│   ├── whistogram.m                <- weighted histogram
│   ├── data/                       settings1-4.mat (3-, 4-, 5-, 6-component models)
│   └── results/                    tunning_results_1-4.mat, metrics_dk.*, figures
│
├── multinomial_posterior/      # Fig4.png — 6-component multinomial, known posterior
│   ├── run_multinomial.m           <- DRIVER (serial: a proposal costs microseconds,
│   │                                  so parfor dispatch would dominate the cost)
│   ├── multinomial_example.m       <- ASMC sampler, all 5 kernels
│   ├── simulateMultinomialSummary.m
│   ├── make_figures_multinomial.m  <- metrics table + assembled figure
│   ├── figStyle.m                  <- shared palette / fonts (copy, see note below)
│   └── results/                    results_multinomial{1,3,4,6,7}.mat, metrics, figure
│
├── hmm_global_alignment/       # Fig5.png — pair-HMM global sequence alignment
│   ├── run_global_aligment_hmm.m       <- DRIVER
│   ├── global_aligment_hmm_example.m   <- ASMC sampler, all 5 kernels
│   ├── simulateHmmAlignments.m         <- ancestral sampler over M alignments
│   ├── getHmmSummary.m                 <- the q = 17 summary statistics
│   ├── simulateHmmSummary.m
│   ├── setupHmmPool.m                  <- parallel pool for the proposal loop
│   ├── test_hmm_predictions.m          <- posterior predictive check (diagnostic only;
│   │                                      no figure of the paper depends on it)
│   ├── make_figures_hmm.m              <- metrics table + assembled figure
│   ├── figStyle.m                      <- shared palette / fonts (copy, see note below)
│   └── results/                        results_hmm_full{1,3,4,6,7}.mat, metrics, figure
│
└── COVID-19 model/             # Fig6.png + Table 1 — COVID-19 screening Bayes net
    ├── run_bayes_net.m             <- DRIVER
    ├── bayes_net_example.m         <- ASMC sampler, all 5 kernels
    ├── build_net.m                 <- defines the 15-node DAG -> Bayes_net.mat
    ├── getCovidSummary.m           <- the q = 31 summary statistics
    ├── simulateCovidSummary.m      <- cohort simulator (reads theta directly)
    ├── getCovidOffsets.m           <- flat-vector offsets of each node's CPD
    ├── predictCovidProb.m          <- closed-form P(COV = + | symptoms, odours)
    ├── classifier_performance.m    <- Table 1 and Supp. Table S6: metrics + 95% credible intervals
    ├── youdenThreshold.m           <- operating threshold shared with make_figures_covid
    ├── make_network_figure.m       <- drawing of the DAG (diagnostic)
    ├── initializePriors_final_sim.m
    ├── setupCovidPool.m
    ├── make_figures_covid.m        <- metrics table + assembled figure
    ├── figStyle.m                  <- shared palette / fonts (copy, see note below)
    ├── xdata.mat, dataset.xlsx     KOR-test screening data
    ├── Bayes_net.mat               DAG structure written by build_net.m
    ├── results/                    results_bayes_net{1,3,4,6,7}.mat, metrics, figures
    └── bnt-master/                 vendored Bayes Net Toolbox; used only to build the
                                    DAG and to cross-check the closed-form predictions
```

---

## Case studies

| Folder | Model | Parameters on the simplex | Paper item |
|---|---|---|---|
| [`DK_tuning/`](DK_tuning/) | Multinomial likelihood, 3–6 categories, 3 dataset sizes | one Dirichlet vector, `D = 3…6` | **Fig3.png**, Fig. S1 |
| [`multinomial_posterior/`](multinomial_posterior/) | 6-component multinomial, `n = 20` trials | one Dirichlet vector, `D = 6` | **Fig4.png** |
| [`hmm_global_alignment/`](hmm_global_alignment/) | Pair-HMM for global genomic alignment | two blocks, `D = 2` and `D = 3` | **Fig5.png** |
| [`COVID-19 model/`](COVID-19%20model/) | 15-node COVID-19 screening Bayes net with a hidden olfactory-status node | 35 blocks, 37 free parameters | **Fig6.png**, **Table 1**, Table S6 |

### Which script produces which item of the paper

| Paper item (printed number) | File | Produced by |
|---|---|---|
| Fig. 2A–B | `Fig3.png` | `DK_tuning/make_figures_dk.m` |
| Fig. 2C | `Fig3.png` | `DK_tuning/ABC_solution_comparison.m` |
| Fig. 3 | `Fig4.png` | `multinomial_posterior/make_figures_multinomial.m` |
| Fig. 4 | `Fig5.png` | `hmm_global_alignment/make_figures_hmm.m` |
| Fig. 5 | `Fig6.png` | `COVID-19 model/make_figures_covid.m` |
| Table 1 | — | `COVID-19 model/classifier_performance.m` (output `T`, Youden threshold) |
| Supplementary Table S6 | — | `COVID-19 model/classifier_performance.m` (`T` and `T50`, the fixed 0.5 threshold) |
| Supplementary Table S3 | — | collated from the three `metrics_<model>.csv` files |
| Supplementary Table S5 | — | the q = 17 statistics of `hmm_global_alignment/getHmmSummary.m`, typeset by hand |
| Supplementary Tables S2, S4 | — | settings reported by hand from `DK_tuning/data/settings*.mat` and the drivers |
| Supplementary Table S1 | — | not reproduced here: only the adopted scaling, `γ = 1/D`, is implemented (`getKernelParam.m`, method 2) |
| Supplementary Fig. S1 | `FigS1_exact.png` | `DK_tuning/make_figure_S1.m` |

The assembled figures leave room for the model schematics, which are added to the final artwork by hand: panel A is blank in `Fig5.png` and `Fig6.png` as written by the scripts. `make_figures_multinomial.m` likewise writes a blank panel A followed by panels B–E, whereas the published `Fig4.png` has no schematic and carries those four panels as A–D; the relettering was done on the artwork.

`multinomial_posterior` is the **convergence benchmark**: the counts are sufficient, so at `tolFinal = 0` the ABC posterior is the *exact* posterior `Dir(prior + y)`. Each kernel is scored against it component by component — posterior mean, median, interquartile range and central 95% interval, all computed from the weighted particles and compared with the analytic marginals `Beta(alpha_i, alpha_0 - alpha_i)`.

---

## Kernel names: code ↔ paper

The `kernel` string selects the perturbation strategy; `getKernelSpec(idx)` maps the index used in the drivers and file names onto it.

| Index | `kernel` string | Paper name | What it does |
|---|---|---|---|
| 1 | `'cw_DK'` | **DK** | Dirichlet kernel applied separately to each Dirichlet block. Proposes inside the simplex by construction. |
| 3 | `'MGK'` | **LT-MGK** | Linear Transformation: keeps the first `D-1` components, `U_LT = [I_{D-1}, 0]`. No Jacobian correction; support is the unit simplex. |
| 4 | `'ilrt_MGK'` | **ILRT-MGK** | Isometric Log-Ratio Transformation, `z = U'_ilrt · ln(x)`. Unbounded support in `R^{D-1}`, metric-preserving. |
| 6 | `'beta_logitMGK'` | **BLT-MGK** | Betancourt's stick-breaking map to `D-1` independent Betas, followed by the logit `z_j = ln(y_j/(1-y_j))`. Unbounded support in `R^{D-1}`. |
| 7 | `'gamma_MGK'` | **LGT-MGK** | Log-Gamma Transformation: `z_j = ln(v_j)` on the underlying Gamma variates. Unbounded, no dimension reduction (`D` free variables). |

These are exactly the five strategies reported in the paper. **Indices 2 and 5 are unused**.

Each transformation is set up in a dedicated `if strcmp(kernel, ...)` block at the top of the sampler, where the forward map, inverse map, transformed log-prior `log_pz` and Jacobian determinant `detJ` are defined as anonymous functions.

---

## DK tuning methods

The `method` argument selects how the DK concentration `κ_t` is adapted at each ASMC iteration. Implemented in `getKernelParam.m` (DK_tuning) and `getDKparam.m` (the other folders).

| `method` | Paper name | Implementation |
|---|---|---|
| `1` | **MLE concentration** | `fmincon` (SQP) on the weighted Dirichlet log-likelihood on the log-scale `α = exp(a)`, with the analytic gradient; then `κ = Σ α_MLE − D`. **Best performer in the paper**, and the one used by every `cw_DK` run. |
| `2` | **Minimum KL divergence** | `fminsearch` on `compute_kld`, the KL between the target and a weighted mixture of DKs; returns the deflated estimate `κ = κ_KL / D` (see Supplementary Table S1). |
| `3` | **Total variance** | Weighted variation matrix `T_{l,m} = Var[ln(θ_l/θ_m)]`; returns `κ = 2D / Σ T`. |

`method` is ignored (pass `[]`) for all MGK kernels, whose hyper-parameter is twice the weighted sample covariance — computed in the `else` branch of `getGeneralKernelParam.m`. The covariance is **shrunk towards its block diagonal** when near-singular.

---

## Common ASMC settings

| Variable | Meaning | Value used |
|---|---|---|
| `N` | number of particles | `1e4` in the three main case studies |
| `alpha` | percentage of "alive" particles kept per iteration | `50` → `α_alive = 0.5` |
| `tolStart` | initial tolerance | `Inf` (set from the prior sample) |
| `tolFinal` | final tolerance `ε_T` | `0` (multinomial), `12` (HMM), `35` (COVID) |
| `pmin_acc` | minimum acceptance rate `p_min` | `1e-5`, an **argument** of the sampler |

`alpha` is a **percentage** in the code (`50`) but a fraction (`0.5`) in the paper.

The tolerance schedule is `ε_t = max{ min{ ρ^(N_α), ε_{t-1} − δ_t }, ε_T }`. The percentile `ρ^(N_α)` normally drives it; the `ε_{t-1} − δ_t` term binds only when the percentile fails to move, which keeps the sequence strictly decreasing and guarantees termination. The decrement is multiplicative for the HMM and COVID models and additive (`δ_t = 2`) for the multinomial ones, where `ρ` is an `l1` norm between count vectors of equal total and is therefore always even.

There is **no cap on the number of proposals**. A run stops when it reaches `tolFinal`, or when an iteration is below `pmin_acc` after 1e6 proposals. In the latter case the partial batch is discarded, the previous complete population is restored, and the result is returned with `converged = false`.

Every driver calls `rng('default')` before sampling. The HMM and COVID samplers parallelise the proposal loop with `parfor` over **chunks**, each drawing from its own `mrg32k3a` substream, so results are reproducible and independent of the worker count.

---

## Shared helper functions

Each case-study folder is **self-contained**: it holds its own copy of every helper it needs, so you can `cd` into one and run it without touching the others. `compute_kld.m`, `getPopulationMetrics.m`, `logsumexp.m` and `probDirichlet.m` are duplicated in all four folders, `figStyle.m`, `getDKparam.m` and `mappingFxn_beta.m` in the three benchmark folders. Those copies are identical — if you change one, change the others. `getKernelSpec.m` and `getGeneralKernelParam.m` also exist in the three benchmark folders but differ on purpose: the former in the result-file names, the latter in how the Dirichlet blocks are located (from `priors` in the multinomial and HMM folders, from the block `indexes` in the Bayes net).

| Function | Purpose |
|---|---|
| `probDirichlet(alpha, xdata, flag, weights)` | Dirichlet log-density (`flag = 1`) or weighted log-likelihood for unconstrained optimisation over `ln α` (`flag = 2`). |
| `compute_kld(...)` | Objective for DK tuning method 2. |
| `getDKparam` / `getKernelParam` | The three DK tuning methods (see above). |
| `getGeneralKernelParam(...)` | Dispatcher: per-block DK parameters for `cw_DK`, shrunk weighted sample covariance otherwise. |
| `getKernelSpec(idx)` | Kernel index → code name, tuning method, paper name, result variable and file name. |
| `getPopulationMetrics(pop)` | AAR, ESS, `time/ESS` and the shrinkage summary for one run. |
| `logsumexp(x)` | Numerically stable `log(Σ exp(x))`, used for importance-weight normalisation. |
| `mappingFxn_beta(x)` | Forward Betancourt stick-breaking map (simplex → `D-1` Betas). |
| `whistogram(...)` | Weighted histogram (DK_tuning only; the other folders define a local `wHist`). |
| `figStyle()` | Shared look of the assembled figures: palette, Arial point sizes, kernel → colour map. Sampled from the final `Fig3.png` so that `Fig4.png`–`Fig6.png` print at the same type size and in the same colours. |

Bayes-net-specific helpers:

| Function | Purpose |
|---|---|
| `initializePriors_final_sim(xdata)` | Dirichlet priors for the 15 network nodes, in the flat layout consumed by the sampler. |
| `getParamIndex(priors)` | Start/end column index of each Dirichlet block inside the flat parameter vector. |
| `getTransformedIndex(indexes)` | Same, for the `D-1`-dimensional transformed space. |
| `getCovidOffsets(priors)` | Node offsets that let `simulateCovidSummary` and `predictCovidProb` read `theta` without building CPD tables. |
| `reshapeProbs(...)` | Converts between the flat "sampling" layout and the nested "simulation" (cell array of CPDs) layout handed to BNT. |
| `predictCovidProb(...)` | Closed-form `P(COV = + | symptoms, odours)`, marginalising the hidden olfactory node. |
| `youdenThreshold(p, y)` | Threshold maximising Youden's J; shared by `classifier_performance` and `make_figures_covid` so Table 1 and the ROC panel mark the same operating point. |

BNT is used **only** by `build_net.m` and by the verification step inside `classifier_performance.m`. `simulateCovidSummary.m` simulates and summarises a cohort in one pass, and `predictCovidProb.m` does inference in closed form.

`build_net.m` is a one-off setup script: it writes `Bayes_net.mat`, which `classifier_performance.m` loads for the cross-check. Nothing *calls* it, but deleting it would make that check unreproducible.

---

## Output: `population` structure

| Field | Type | Contents |
|---|---|---|
| `SMCparams` | `N × D` | final particle population **on the simplex** |
| `HMMparams` | `N × 5` | pair-HMM only: the transition parameters derived from `SMCparams` |
| `weights` | cell array | importance weights, one vector per ASMC iteration; `weights{end}` is the posterior |
| `rho` | cell array | discrepancy of each particle, per iteration |
| `tolSMC` | vector | tolerance schedule `{ε_t}` |
| `acceptRate` | vector | acceptance rate per iteration (the paper's *AAR* is the mean over the SMC stage) |
| `ESS` | vector | effective sample size per iteration, `1/Σ w²` |
| `Kparams` | cell | kernel hyper-parameter per iteration (`κ_t` for DK, covariance for MGK) |
| `shrinkDelta` | vector | block-diagonal shrinkage intensity applied per iteration (0 if none) |
| `rejectBudget` | `n × 4` | per iteration: `[proposals, outside simplex, above tolerance, accepted]` |
| `converged` | logical | `true` only if `tolFinal` was reached |
| `finalTol` | scalar | tolerance actually attained |
| `nStats` | scalar | size of the summary vector — stamps the experimental design |
| `time` | scalar | wall-clock sampling time, seconds |
| `kernel`, `method` | string / int | strategy used |

The main performance metric is **`time / ESS(end)`**.

`rejectBudget` separates the two ways a proposal can fail — leaving the simplex versus merely landing too far — which distinguishes an *inadmissible* kernel from an *imprecise* one.

**All posterior summaries are weighted.** The particles carry importance weights and the ESS can be a small fraction of `N`, so quantiles, medians, means and kernel densities are computed considering the weights.

---

## Results files and naming convention

The variable *inside* each file is named after its kernel (`population_cwDK`, `population_MGK`, `population_ilrtMGK`, `population_beta_logitMGK`, `population_gammaMGK`); `getKernelSpec` returns that name in its `varName` field, so scripts never hard-code it.

Alongside them, each `make_figures_*` writes `metrics_<model>.{mat,csv}` and the assembled figure as `.fig` and `.png`.

`Fig3.png`–`Fig6.png` are additionally written as **`.svg`** for typesetting — `fig3A_performance`, `fig3B_concentration` and `fig3C_marginals` in `DK_tuning/results/`, `Fig4` in `multinomial_posterior/results/`, `Fig5` in `hmm_global_alignment/results/` and `Fig6` in `COVID-19 model/results/`. The SVGs are true vector output (`-dsvg -painters`): labels set with the TeX interpreter are converted to glyph outlines, so they render identically without depending on installed fonts. Supplementary Fig. S1 is written as `.png` only (`DK_tuning/results/FigS1_exact_posterior.png`, renamed `FigS1_exact.png` in `figures_final/`).

The assembled figures are copied by hand into `submission/figures_final/`, which is what the LaTeX sources include; nothing in `scripts/` writes there.

For `DK_tuning/results/tunning_results_{1,2,3,4}.mat` (3-, 4-, 5- and 6-component models), variables are `population<method><dataset>` — e.g. `population13` is tuning method 1 on dataset 3.

---

## Requirements

- **MATLAB R2021a** or later.
- **Statistics and Machine Learning Toolbox** — `randg`, `mvnrnd`, `cholcov`, `betapdf`, `betacdf`, `betainv`, `tiedrank`, `prctile`.
- **Optimization Toolbox** — `fmincon`, for DK tuning method 1.
- **Parallel Computing Toolbox** — `parfor` and `parpool`, for the HMM and COVID samplers. Optional: with no pool open, `parfor` runs serially and only the wall time changes.
- **Bioinformatics Toolbox** — `hmmdecode`, for `test_hmm_predictions.m` only. The pair-HMM *sampler* has its own ancestral sampler in `simulateHmmAlignments.m`.
- **[Bayes Net Toolbox (BNT)](https://github.com/bayesnet/bnt)** — for `build_net.m` and for the verification step in `classifier_performance.m`. A copy is vendored at `COVID-19 model/bnt-master/` and added to the path by those scripts.

Reported benchmarks were run on an Intel Core Ultra 7 258V (8 physical cores, 32 GB RAM) with a pool of 8 workers.

---

## Citation

If you use this code, please cite the paper. Data for the COVID-19 Bayes net is available separately (refer to the *Data and code availability* section of the manuscript).