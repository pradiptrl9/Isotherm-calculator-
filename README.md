# Adsorption Isotherm & Chemical Kinetics Analyzer

A scientific computing toolkit for analyzing **adsorption equilibrium isotherms** (Freundlich & Langmuir) and **first-order chemical decomposition kinetics**, available in both **R** and **Python**.

---

## 🔬 Mathematical Formulations

### 1. Adsorption Isotherms

#### A. Freundlich Isotherm (Multilayer / Heterogeneous Adsorption)
- **Non-linear equation**:
  $$q_e = K_f \cdot C_e^{1/n}$$
- **Linearized log-log form**:
  $$\log_{10}(q_e) = \frac{1}{n}\log_{10}(C_e) + \log_{10}(K_f)$$
- **Interpretation**:
  - $n > 1$: Favorable adsorption process
  - $n = 1$: Linear partition adsorption
  - $n < 1$: Unfavorable adsorption process

#### B. Langmuir Isotherm (Monolayer / Homogeneous Adsorption)
- **Non-linear equation**:
  $$q_e = \frac{q_{\max} K_L C_e}{1 + K_L C_e}$$
- **Linear form**:
  $$\frac{C_e}{q_e} = \frac{1}{q_{\max} K_L} + \frac{1}{q_{\max}} C_e$$
- **Dimensionless Separation Factor ($R_L$)**:
  $$R_L = \frac{1}{1 + K_L C_0}$$
  - $0 < R_L < 1$: Favorable adsorption
  - $R_L > 1$: Unfavorable adsorption
  - $R_L = 1$: Linear adsorption
  - $R_L = 0$: Irreversible adsorption

---

### 2. First-Order Reaction Kinetics

- **Differential Rate Law**: $\frac{d[A]}{dt} = -k [A]$
- **Integrated Rate Equation**: $[A]_t = [A]_0 \cdot e^{-k t}$
- **Linearized Form**: $\ln([A]_t) = -k t + \ln([A]_0)$
- **Half-Life**: $t_{1/2} = \frac{\ln(2)}{k} \approx \frac{0.69315}{k}$
- **Shelf-Life**: $t_{90} = \frac{\ln(1/0.90)}{k} \approx \frac{0.10536}{k}$

---

## 📁 Repository Structure

| File | Description |
|------|-------------|
| `isotherm_calculator.R` | Complete Adsorption Isotherm calculator written in native R (base R graphics, zero external dependencies). |
| `isotherm_calculator.py` | Complete Adsorption Isotherm calculator written in Python with NumPy & Matplotlib. |
| `first_order_kinetics.py` | First-order chemical kinetics rate fitting, half-life calculator, and statistical OLS regression. |
| `kinetics_calculator.html` | Interactive web-based browser UI for chemical kinetics calculations and visualizations. |

---

## 🚀 Usage

### Running in R
Run directly in R or RStudio:
```r
source("isotherm_calculator.R")
```
Or via terminal:
```bash
Rscript isotherm_calculator.R
```

### Running in Python
```bash
python isotherm_calculator.py
python first_order_kinetics.py
```

---

## 📊 Features & Outputs

- **Statistical Rigor**: Computes both linear $R^2$ (in transformed space) and real-space non-linear $R^2$ plus RMSE.
- **Model Comparison**: Automatically compares Freundlich vs. Langmuir to determine the best-fit model for experimental data.
- **Predictions**: Calculates expected adsorption capacity $q_e$ for new arbitrary equilibrium concentrations $C_e$.
- **High-Resolution Plots**: Outputs dual-panel comparison graphs showing experimental points, non-linear fitted curves, and linearized regression plots.
