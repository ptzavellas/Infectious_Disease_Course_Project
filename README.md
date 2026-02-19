# Infectious Disease Modeling & Bayesian Change-Point Analysis

A comprehensive epidemiological modeling project analyzing COVID-19 transmission dynamics in the Czech Republic using deterministic and stochastic models, including SEIR modeling and Bayesian chain binomial inference.

---

## Objective

To investigate the transmission dynamics of COVID-19 through:

- Epidemic curve analysis
- Effective reproduction number (Rt) estimation
- Deterministic SEIR modeling
- Age-structured contact modeling
- Bayesian chain binomial modeling
- Change-point detection in transmission rates

---

## Dataset

- COVID-19 daily case data (Czech Republic)
- Population size: 10,693,939
- Age-stratified contact matrix
- Serial interval assumptions:
  - Mean: 4 days
  - SD: 2.8 days

---

## Part 1 – Epidemic Curve & Rt Estimation

- Constructed epidemic curve for first 15 days
- Estimated time-varying reproduction number (Rt)
- Observed Rt between approximately 2.3–3.2 (indicating rapid epidemic growth)

---

## Part 2 – Deterministic SEIR Modeling

Implemented SEIR model with parameters:

- R0 = 2.6
- β = 0.371
- γ = 0.14
- σ = 0.25

Simulated:

- Number of exposed individuals over time
- Cumulative infections
- Final epidemic size
- Scenario comparison (baseline vs 40% contact reduction)

Key result:
Without intervention, ~45% of the population infected within 4 months.

---

## Part 3 – Age-Structured Transmission

Used contact matrix to model transmission among:

- Children (20%)
- Adults (50%)
- Seniors

Observed:
- Faster early spread among children
- Adults driving overall epidemic magnitude
- Slower but prolonged spread among seniors

---

## Part 4 – Bayesian Chain Binomial Model (nimble)

Implemented stochastic transmission model with:

- Mean infectious period: 7 days
- β ≈ 0.155
- MCMC estimation using nimble

Produced:
- Posterior density plots
- Trace plots
- Autocorrelation diagnostics

---

## Part 5 – Change-Point Modeling

### Fixed Change-Points
Analyzed multiple predefined dates:
- 10-10-2020
- 05-12-2020
- 08-01-2021
- 21-02-2021
- 15-03-2021

Estimated separate transmission rates (β₁, β₂, ...)

### Estimated Change-Point (Bayesian)
- Discrete uniform prior over all dates
- Posterior identified change-point: 04-09-2020
- β reduced from ~0.756 to ~0.155

---

## Methods Demonstrated

- Epidemic curve construction
- Rt estimation
- SEIR differential equation modeling
- Age-structured contact modeling
- Bayesian inference with nimble
- MCMC diagnostics (trace plots, ACF, posterior densities)
- Change-point detection
- Scenario analysis

---

## Tools Used

- R
- EpiEstim (Rt estimation)
- deSolve (SEIR differential equations)
- nimble & rstan (Bayesian inference)
- coda & bayesplot (MCMC diagnostics)
- ggplot2 & viridis (visualization)
- dplyr & reshape2 (data manipulation)
- socialmixr (contact matrices)
- R0 (basic reproduction number estimation)
- rootSolve
- shinySIR

---

## Reproducibility

The workflow includes:

1. Data loading and preprocessing
2. Rt estimation
3. SEIR modeling
4. Age-structured modeling
5. Bayesian chain binomial inference
6. Change-point estimation
