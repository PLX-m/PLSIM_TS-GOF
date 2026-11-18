# Goodness-of-Fit Testing for High-Dimensional Partial Linear Single-Index Models under Weak Dependence

This repository provides the R/Python code and numerical results accompanying the manuscript **“Goodness-of-Fit Testing for High-Dimensional Partial Linear Single-Index Models under Weak Dependence”** (working title).

The study develops goodness-of-fit tests for **partial linear single-index models (PLSIMs)** under weak dependence, including high-dimensional settings. The repository contains simulation scripts, saved Monte Carlo p-values, plotting code, and a real-data example. The implementation combines R-based statistical procedures with Python/PyTorch models for conditional-expectation estimation.

## Repository structure

- `code/`: single-index estimation, orthogonalization, goodness-of-fit testing, and neural-network utilities.
- `Example123.R`: low-dimensional experiments (Examples 1–3).
- `Example_HD123.R`: high-dimensional experiments (Models HD1–HD3).
- `Simulation A.R`: temporal and cross-covariate dependence experiments.
- `Simulation B.R`: experiments with serially correlated errors.
- `real data.R`: analysis of the Hong Kong hospital admissions dataset from `VGAMextra`.
- `results/`: 82 saved `.RData` files containing simulation-level p-values (30 low-dimensional, 12 high-dimensional, 32 Simulation A, and 8 Simulation B).
- `plot/`: plotting scripts and saved PDF figures for the low-dimensional experiments.

## Dependencies and installation

R and Python are both required. The R scripts call `code/python_ts.py` using the [`reticulate`](https://rstudio.github.io/reticulate/) package.

**Python:** Install the directly used packages into the environment you intend to use:

```bash
python -m pip install -r requirements.txt
```

`requirements.txt` lists `numpy`, `torch`, and `scikit-learn`. Select a suitable CPU or CUDA-enabled PyTorch build for your machine.

**R:** From the repository root, run:

```r
source("install_R_packages.R")
```

The R installer checks for the required packages and installs those that are missing. The principal R packages are `splines2`, `MASS`, `openxlsx`, `caret`, `glmnet`, `randomForest`, `reticulate`, `doParallel`, `foreach`, `doRNG`, `ggplot2`, `patchwork`, and `VGAMextra`.

**Before running an experiment:**

1. Open the selected R script and change any machine-specific `setwd()`, `use_python()`, and `py_bin` paths. Run with the project root as the working directory so that `code/...` and `results/...` paths resolve.
2. Make sure that `reticulate` points to the **same Python environment** in which the Python packages were installed. Configure it before Python is initialized in the R session.
3. Select the intended model, dimensions, hypothesis, and simulation parameters. The code contains illustrative run settings; they do not automatically execute every configuration reported in the paper.
4. Use a new output filename/directory when rerunning experiments: `save()` can overwrite an existing `.RData` file. The original results should be treated as read-only.

For example, to verify the Python configuration in R:

```r
reticulate::use_python("path/to/your/python", required = TRUE)
reticulate::py_config()
```

## Reproducing the simulation studies

Each reported configuration uses **1,000 Monte Carlo replications**. Empirical rejection rates are computed at the nominal significance level **0.05**. The files in `results/` allow these rates to be checked without rerunning the simulations.

### Examples 1–3 and Models HD1–HD3

The low- and high-dimensional experiments are implemented in `Example123.R` and `Example_HD123.R`, with saved results in `results/different_c/` and `results/high/`, respectively.

| Example / Model | Baseline single-index function `g(t)` | Null scenario | Alternative scenario |
| --- | --- | --- | --- |
| Example 1 / HD1 | `t^2` | `H0` | `case1` |
| Example 2 / HD2 | `0.5*exp(-t^2)` | `H0` | `case2` |
| Example 3 / HD3 | `cos(2*t) + 0.5*t^2` | `H0` | `case3` |

Low-dimensional settings use `N = 800`, `p = q = s ∈ {10, 20}`, and deviation magnitudes `c ∈ {0, 0.25, 0.5, 0.75, 1}`. High-dimensional settings use `N = 1000`, `p = q ∈ {1000, 2000}`, and `s = 20`, with size under `H0` and power under the matching alternative at `c = 1`.

**Important:** For each example/model, manually select the matching `Y_base` expression in `generate_PLSIM_ts()` by commenting or uncommenting the relevant line. Changing `scenario` alone does not switch the baseline function. The low-dimensional plotting scripts are in `plot/`.

### Simulation A: dependence in covariates

`Simulation A.R` (`results/Simulation_A/`) examines temporal dependence `rho ∈ {0.5, 0.8}` and cross-covariate dependence `kappa ∈ {0, 0.3, 0.6}` across `(N, p, q, s) = (800, 10, 10, 10)`, `(1000, 1000, 1000, 20)`, and `(2000, 2000, 2000, 20)`. It uses the Example 3 / HD3 baseline, with `scenario = "H0"` for size and `"case3"` for power; the alternative deviation is selected according to dimension.

In filenames in this folder, `sigma03` and `sigma06` indicate `kappa = 0.3` and `0.6`, **not** the error standard deviation.

### Simulation B: serial correlation in errors

`Simulation B.R` (`results/Simulation_B/`) uses `rho = 0.5`, `kappa = 0.3`, and AR(1) error correlation `phi ∈ {0, 0.5, 0.8}`. It considers the same low-dimensional setting and the `N = p = q = 1000` high-dimensional setting from Simulation A. Use `scenario = "H0"` for size and `"low"` / `"high"` for power. The `phi = 0` reference results are shared with Simulation A; this folder contains the nonzero-`phi` results.

## Reproducing reported rejection rates from saved results

The `.RData` files store an object named `res_size` with `settings` and `raw_results`. The three columns of `raw_results` are `p_12`, `p_21`, and `p_CCT`, corresponding to the two component tests and their Cauchy combination.

To compute the empirical size or power from an existing result file **without overwriting objects in the global R environment**:

```r
result_env <- new.env()
load("results/different_c/t_val^2_n800_d10_s10_nsims1000_size.RData",
     envir = result_env)

rr <- result_env$res_size$raw_results
colMeans(rr < 0.05, na.rm = TRUE)
```

For an alternative-model result, use the corresponding power `.RData` filename. Saved values are simulation-level p-values, so the three rejection rates can be recalculated for other significance thresholds as well.


## Computing environments

Simulations were run on two computers. The following environments were **inspected on 18 November 2026**; they are not claimed to reproduce the exact historical configuration of every simulation run.

| Item | PC1 | PC2 |
| --- | --- | --- |
| OS | Windows 11 x64 | Windows 10 x64 |
| CPU cores (physical / logical) | 6 / 12 | 28 / 56 |
| R | 4.4.1 | 4.5.2 |
| Python | 3.13.5 | 3.10.19 |
| NumPy | 2.3.1 | 2.2.6 |
| PyTorch | 2.8.0+cpu | 2.6.0+cu124 |
| scikit-learn | 1.7.1 | 1.7.2 |
| CUDA available at inspection | No | Yes (12.4; RTX 4090 detected) |

Recorded versions of directly used R packages:

| R package | PC1 | PC2 |
| --- | --- | --- |
| `splines2` | 0.5.4 | 0.5.4 |
| `MASS` | 7.3-60.2 | 7.3-65 |
| `openxlsx` | 4.2.6.1 | 4.2.8.1 |
| `caret` | 7.0-1 | 7.0-1 |
| `glmnet` | 4.1-9 | 4.1-10 |
| `randomForest` | 4.7-1.2 | 4.7-1.2 |
| `reticulate` | 1.43.0 | 1.44.1 |
| `doParallel` | 1.0.17 | 1.0.17 |
| `foreach` | 1.5.2 | 1.5.2 |
| `doRNG` | 1.8.6.2 | 1.8.6.2 |
| `ggplot2` | 4.0.1 | 4.0.1 |
| `patchwork` | 1.3.2 | Not installed |
| `VGAMextra` | 0.0-9 | Not installed |

The R parallel workers run CPU-based statistical procedures. The PyTorch network-training functions select `cuda` when CUDA is available and otherwise use `cpu`. Therefore, the supplied code *would* choose CUDA for the neural-network portion on PC2 when the GPU is available; **actual GPU use in past simulation runs has not been independently verified**. No result file is assigned to a specific computer based solely on these current environment snapshots.

## Reproducibility notes

- The supplied simulation results can be used to verify the reported rejection rates without rerunning expensive Monte Carlo experiments.
- Saved result files record p-values and some run settings; they are **not** full records of the original generated datasets or exact historical computational environments. Some simulation parameters are encoded primarily in filenames.
- The baseline link function for Examples 1–3 and HD1–HD3 must be changed manually as explained above.
- Differences in dependency versions, Python device selection, numerical libraries, and parallel random-number generation can affect rerun results. Matching the saved rejection rates is not the same as guaranteeing bit-for-bit identical reruns.

