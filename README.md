# 📈 Log-Bilal ARMA

This repository presents the **Log-Bilal-ARMA model**, an innovative approach for modeling **time series data bounded in the unit interval (0, 1)** using the **Log-Bilal distribution** as its probabilistic foundation.

## 💡 About the Model

The Log-Bilal-ARMA model combines:

- The flexibility of the **Log-Bilal distribution**, obtained through the transformation \( Y = \exp(-X) \), where \( X \sim \text{Bilal} \);
- An **autoregressive moving average (ARMA)** structure to capture temporal dependence.

This framework enables parsimonious modeling of unit-valued dynamic data, such as rates, proportions, and probabilities observed over time.

## 📂 Repository Contents

- `LogB-functions.r` – Auxiliary functions for the model.
- `LogBarma.fit` – Main function for model fitting.
- Simulation scripts, model fitting routines, and Monte Carlo experiments.
- Numerical results and illustrative plots.

## 📊 Applications

The model is suitable for scenarios involving continuous data in the (0, 1) interval, such as:

- **Rates**;
- **Proportions**;

## 🛠️ Work in Progress

This project is under active development. New features and testing procedures are being added regularly.

