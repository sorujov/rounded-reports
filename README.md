# Identification and inference for the shape of a distribution from rounded reports

Samir Orujov, School of Business, ADA University, Baku

Survey incomes are reported in round numbers, and tests of unimodality applied to such reports reject for every unimodal income distribution simulated in the paper. The paper studies what the reports identify. Conditional on the rounding level, the greatest common divisor of a report with the least common multiple of the candidate levels has an explicit distribution involving Euler's totient, so the rounding probabilities are identified and the rounding mechanism can be tested without nuisance parameters. At widths that are common multiples of the active levels the cell probabilities of latent income are identified up to a second-order boundary term, and unimodality of the cell vector is tested with a moment-selection bootstrap whose asymptotic level is at most nominal at every point of the null. Finer resolution requires inverting the reporting map. The methods are applied to American Community Survey and Current Population Survey income data.

Keywords: heaping; rounding; modality; excess mass; shape-constrained inference; moment selection. JEL: C12, C14, C81, D31.

Working paper: [SSRN 6454722](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6454722).

## Contents

| Folder | What it holds |
|---|---|
| `manuscript/` | LaTeX source of the paper (`main.tex`), the online supplement (`supplement.tex`) and the title page |
| `code/` | The procedures of the paper and the simulation study |
| `empirics/` | The applications of Section 7, with their outputs as CSV files |
| `data/` | Where the public-use Census files go (not redistributed; see below) |

### code/

- `amlib.R`, `emcpp.cpp`: fit of the reporting mechanism (Proposition 1), specification test (Theorem 1), cells with boundary splitting (Lemma 1), exact discrete excess mass (C++ via Rcpp), moment-selection bootstrap (Theorem 2), and `modality_test()`, which runs the whole procedure.
- `rounding.R`, `heaping.R`, `benchmarks.R`, `deconv.R`, `localprof.R`: benchmark procedures, the reporting map and its inversion, and the band estimator.
- `dgps_heap.R`, `flatdgp.R`: latent designs and reporting mechanisms of Section 6.
- `sim_main.R` (Table 1 and Tables S3 to S5), `sim_rer.R` (re-rounding column), `sim_aux.R` (Tables S1 and S2), `sim_local2.R` and `sim_local3.R` (Table S6). Each script runs the block of Monte Carlo units given by the environment variables `UNIT_START`, `UNIT_END`, `WORKERS` and `RESULTS_DIR`. The `.sbatch` files are the Slurm array jobs used to run the full study.

### empirics/

- `emp_earn.R` (Table 2), `emp_housing.R` (Table 3), `emp_cps1.R` (Table 4), `cps_robust.R` (coarser reporting widths, Section 7.3), `kappa_acs.R` (condition numbers, Section 5) and `mk_tables.R` (LaTeX rows of the tables).
- `emp_earn.csv`, `emp_housing.csv`, `emp_cps1.csv`: the results reported in the paper.

Run the scripts from inside `empirics/`; they read the library from `../code/` and the data from `../data/`.

## Data

All data are public-use files of the U.S. Census Bureau, available from <https://www2.census.gov>:

- American Community Survey 2022 one-year PUMS person files for California, Mississippi and Puerto Rico, and the California housing file, saved as compressed CSV in `data/emp/` (`acs22_person_ca.csv.gz`, `acs22_housing_ca.csv.gz`, and so on);
- Current Population Survey ASEC 2023 person file, and the ASEC 2017 to 2021 public-use files, from which the one-person household extracts `data/cps/onep_<year>.csv` are built.

## Requirements

R 4.2 or later with the packages Rcpp, data.table, multimode and diptest. The manuscript compiles with pdfLaTeX and BibTeX:

```
cd manuscript
pdflatex supplement && pdflatex main && bibtex main && bibtex supplement
pdflatex supplement && pdflatex main && pdflatex supplement && pdflatex main
```

The main text and the supplement refer to each other through `xr-hyper`, hence the repeated runs.

## Citation

Orujov, S. (2026). Identification and inference for the shape of a distribution from rounded reports. Working paper, ADA University. SSRN 6454722.

## Licence

Code: MIT licence (see `LICENSE`). The text of the paper remains the author's; please cite it rather than reuse it.
