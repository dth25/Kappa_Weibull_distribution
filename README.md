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

The function MAIN_KAPPA_WEIBULL calculates the fit of the data set DATA to the WEIBULL and the KAPPA Weibull distributions.  The empirical CDF of the data is first constructed.

The kappa-Weibull CDF fit is performed using the MLE method by minimization of the Negative Log-Likelihood (NLL) of the KAPPA Weibull distribution.

The optimization uses a two-step approach: First, it employs the function FMINCON using the ACTIVE-SET method and is then followed by FMINSEARCH (Simplex search method), which uses the output of FMINCON as initial condition.  Gradient information is supplied to FMINCON via the objective function lnL_wbk

Plots of the CDF for the empirical distribution as well as the best-fit Weibull and kappa-Weibull distributions are generated. The  functions  Phi(tau) (Weibull plots) for the three CDFs are also plotted. Finally, the quantile-quantile plot (empirical to kappa-Weibull) is generated.

The code has been tested with Matlab R2015b

-----------------------------------------------------------------------------------------------------------------------------------

 INPUT VARIABLES

-----------------------------------------------------------------------------------------------------------------------------------

 DATA:           Array of data values

 DATANAM:        String array containing name of data (used for labeling)

-----------------------------------------------------------------------------------------------------------------------------------

 OUTPUT VARIABLES

-----------------------------------------------------------------------------------------------------------------------------------

 BETA_W:         Parameters of Weibull

 BETA_KW:        Parameters of Kappa-Weibull

 NLL_KW:         Negative Log-Likelihood of Fitted Kappa-Weibull

 NLL_W:          Negative Log-Likelihood of Fitted Weibull

======================================================================

EXAMPLE:

load Cairo_wind_speed;

[beta_w, beta_kw, NLL_W, NLL_KW] = Main_kappa_Weibull(wind, 'Cairo wind');

load Carbon_fiber

[beta_w, beta_kw, NLL_W, NLL_KW] = Main_kappa_Weibull(data_carbon, 'Tensile strength')

-----------------------------------------------------------------------------------------------------------------------------------

EXTERNAL FUNCTIONS USED: expk, lnk, lnL_wbk, wblkcdf

              Email: dchristopoulos@tuc.gr (Dionisis Hristopulos)

     Last Modified: July 30, 2022

-----------------------------------------------------------------------------------------------------------------------------------

 REFERENCES (If you use this code please cite the following)

-----------------------------------------------------------------------------------------------------------------------------------

[1] D.T. Hristopulos and A. Baxevani, “Kaniadakis functions beyond statistical mechanics: weakest-link scaling, power-law tails, and modified lognormal distribution,” Entropy, 2022.

[2] D. T. Hristopulos, M. P. Petrakis, and G. Kaniadakis, “Finite-size effects on return interval distributions for weakest-link-scaling Systems,” Physical Review E 89, 052142, 28 May 2014.

https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.052142

[3] D. T. Hristopulos, M. P. Petrakis, and G. Kaniadakis, “Weakest-link scaling and extreme events in finite-sized systems,” Entropy, 17(3):1103-1122, 2015. https://doi.org/10.3390/e17031103

For more information on the datasets used to test this code see references in [1].

