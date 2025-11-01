# Kappa-Weibull Parameter Estimation in MATLAB

This repository contains MATLAB code for estimating the parameters of the **kappa-Weibull distribution**, as well as the standard **Weibull distribution**, using **maximum likelihood estimation (MLE)**. The code was developed by **Dionissios Hristopulos** and is archived on Zenodo with a DOI for citation.

🔗 **Zenodo DOI**: [10.5281/zenodo.6942812](https://doi.org/10.5281/zenodo.6942812)

## 📦 Contents

- `Main_kappa_Weibull.m`: Main function to fit data to Weibull and kappa-Weibull distributions.
- `expk.m`, `lnk.m`, `lnL_wbk.m`, `wblkcdf.m`: Supporting functions for distribution fitting and likelihood calculation.
- Example datasets: `Cairo_wind_speed`, `Carbon_fiber`
- Output: CDF plots, Weibull plots, and quantile-quantile plots.

## ⚙️ Requirements

- MATLAB R2015b or later
- Optimization Toolbox (for `fmincon`)

## 🚀 Usage

```matlab
load Cairo_wind_speed
[beta_w, beta_kw, NLL_W, NLL_KW] = Main_kappa_Weibull(wind, 'Cairo wind');

load Carbon_fiber
[beta_w, beta_kw, NLL_W, NLL_KW] = Main_kappa_Weibull(data_carbon, 'Tensile strength');
