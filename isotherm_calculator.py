"""
===============================================================================
ADSORPTION ISOTHERM CALCULATOR & DATA ANALYZER (Freundlich & Langmuir)
===============================================================================
Author: Adsorption Chemistry & Separation Modeling
Description:
  Calculates and fits experimental adsorption equilibrium data to standard
  isotherm models:
  1. Freundlich Isotherm:  q_e = K_f * C_e^(1/n)
     Linear form:          log10(q_e) = (1/n) * log10(C_e) + log10(K_f)
  2. Langmuir Isotherm:    q_e = (q_max * K_L * C_e) / (1 + K_L * C_e)
     Linear form:          C_e / q_e = (1 / (q_max * K_L)) + (1 / q_max) * C_e
     Separation factor:    R_L = 1 / (1 + K_L * C_0)
===============================================================================
"""

from typing import Tuple, Dict, Any, Optional
import numpy as np

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def freundlich_fit(C_e: np.ndarray, q_e: np.ndarray) -> Dict[str, Any]:
    """
    Fit Freundlich isotherm parameters using linear regression on log10-log10 data.
    
    Parameters
    ----------
    C_e : array-like
        Equilibrium liquid-phase concentrations (e.g., mg/L or mmol/L)
    q_e : array-like
        Equilibrium solid-phase adsorption capacities (e.g., mg/g or mmol/g)
        
    Returns
    -------
    dict
        Fitted parameters: K_f, n, one_over_n, r2_linear, r2_nonlinear, rmse
    """
    C_e = np.asarray(C_e, dtype=float)
    q_e = np.asarray(q_e, dtype=float)
    
    # Filter out non-positive values (log undefined)
    mask = (C_e > 0) & (q_e > 0)
    C_e_valid = C_e[mask]
    q_e_valid = q_e[mask]
    
    if len(C_e_valid) < 2:
        raise ValueError("Need at least 2 valid data points with C_e > 0 and q_e > 0.")
    
    log_C = np.log10(C_e_valid)
    log_q = np.log10(q_e_valid)
    
    # Linear fit: log_q = (1/n)*log_C + log10(K_f)
    slope, intercept = np.polyfit(log_C, log_q, 1)
    
    if slope == 0:
        raise ValueError("Zero slope in log-log data. Cannot calculate Freundlich exponent n.")
    if slope < 0:
        print("[Warning] Slope is negative (1/n < 0). Capacity decreases with concentration, which is abnormal.")

    one_over_n = slope
    n = 1.0 / slope
    K_f = 10.0 ** intercept
    
    # Linear R² (in log10 space)
    log_q_pred = slope * log_C + intercept
    ss_res_log = np.sum((log_q - log_q_pred) ** 2)
    ss_tot_log = np.sum((log_q - np.mean(log_q)) ** 2)
    r2_linear = 1.0 - (ss_res_log / ss_tot_log) if ss_tot_log != 0 else 0.0

    # Non-linear R² and RMSE (in original q_e space)
    q_pred = predict_freundlich(C_e_valid, K_f, n)
    ss_res = np.sum((q_e_valid - q_pred) ** 2)
    ss_tot = np.sum((q_e_valid - np.mean(q_e_valid)) ** 2)
    r2_nonlinear = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
    rmse = np.sqrt(np.mean((q_e_valid - q_pred) ** 2))
    
    return {
        "model": "Freundlich",
        "K_f": K_f,
        "n": n,
        "1/n": one_over_n,
        "slope": slope,
        "intercept": intercept,
        "r2_linear": r2_linear,
        "r2_nonlinear": r2_nonlinear,
        "rmse": rmse,
        "C_e": C_e_valid,
        "q_e": q_e_valid,
        "q_pred": q_pred
    }


def langmuir_fit(C_e: np.ndarray, q_e: np.ndarray) -> Dict[str, Any]:
    """
    Fit Langmuir isotherm parameters using linear regression:
    C_e / q_e = (1 / (q_max * K_L)) + (1 / q_max) * C_e
    
    Returns
    -------
    dict
        Fitted parameters: q_max, K_L, r2_linear, r2_nonlinear, rmse
    """
    C_e = np.asarray(C_e, dtype=float)
    q_e = np.asarray(q_e, dtype=float)
    
    mask = (C_e > 0) & (q_e > 0)
    C_e_valid = C_e[mask]
    q_e_valid = q_e[mask]
    
    if len(C_e_valid) < 2:
        raise ValueError("Need at least 2 valid data points with C_e > 0 and q_e > 0.")
    
    y = C_e_valid / q_e_valid
    slope, intercept = np.polyfit(C_e_valid, y, 1)
    
    if slope <= 0:
        raise ValueError("Non-positive slope in Langmuir linearization; monolayer saturation capacity cannot be calculated.")
    
    q_max = 1.0 / slope
    K_L = slope / intercept if intercept != 0 else float("inf")
    
    # Linear R²
    y_pred = slope * C_e_valid + intercept
    ss_res_lin = np.sum((y - y_pred) ** 2)
    ss_tot_lin = np.sum((y - np.mean(y)) ** 2)
    r2_linear = 1.0 - (ss_res_lin / ss_tot_lin) if ss_tot_lin != 0 else 0.0
    
    # Non-linear R² & RMSE
    q_pred = predict_langmuir(C_e_valid, q_max, K_L)
    ss_res = np.sum((q_e_valid - q_pred) ** 2)
    ss_tot = np.sum((q_e_valid - np.mean(q_e_valid)) ** 2)
    r2_nonlinear = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
    rmse = np.sqrt(np.mean((q_e_valid - q_pred) ** 2))
    
    return {
        "model": "Langmuir",
        "q_max": q_max,
        "K_L": K_L,
        "slope": slope,
        "intercept": intercept,
        "r2_linear": r2_linear,
        "r2_nonlinear": r2_nonlinear,
        "rmse": rmse,
        "C_e": C_e_valid,
        "q_e": q_e_valid,
        "q_pred": q_pred
    }


def predict_freundlich(C_e: np.ndarray, K_f: float, n: float) -> np.ndarray:
    """Predict q_e using Freundlich equation: q_e = K_f * C_e^(1/n)"""
    C_e = np.asarray(C_e, dtype=float)
    return K_f * np.power(np.maximum(C_e, 0), 1.0 / n)


def predict_langmuir(C_e: np.ndarray, q_max: float, K_L: float) -> np.ndarray:
    """Predict q_e using Langmuir equation: q_e = (q_max * K_L * C_e) / (1 + K_L * C_e)"""
    C_e = np.asarray(C_e, dtype=float)
    return (q_max * K_L * C_e) / (1.0 + K_L * C_e)


def plot_isotherms(
    C_e: np.ndarray, 
    q_e: np.ndarray, 
    freundlich_res: Dict[str, Any], 
    langmuir_res: Optional[Dict[str, Any]] = None,
    save_filename: Optional[str] = None
) -> None:
    """
    Plot experimental adsorption data against fitted Freundlich and Langmuir isotherms.
    """
    if not HAS_MATPLOTLIB:
        print("[Info] Matplotlib is not installed. To view plots, run: pip install matplotlib")
        return

    C_e = np.asarray(C_e, dtype=float)
    q_e = np.asarray(q_e, dtype=float)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    C_smooth = np.linspace(max(C_e.min() * 0.5, 0.01), C_e.max() * 1.15, 250)
    
    # 1. Non-linear Isotherm Curve
    ax1.scatter(C_e, q_e, color="crimson", s=65, edgecolor="black", label="Experimental Data", zorder=5)
    
    q_fr_smooth = predict_freundlich(C_smooth, freundlich_res["K_f"], freundlich_res["n"])
    ax1.plot(
        C_smooth, q_fr_smooth, "b-", linewidth=2.2, 
        label=f"Freundlich ($K_f$={freundlich_res['K_f']:.3f}, $n$={freundlich_res['n']:.3f}, $R^2$={freundlich_res['r2_nonlinear']:.4f})"
    )
    
    if langmuir_res:
        q_lg_smooth = predict_langmuir(C_smooth, langmuir_res["q_max"], langmuir_res["K_L"])
        ax1.plot(
            C_smooth, q_lg_smooth, "g--", linewidth=2.2,
            label=f"Langmuir ($q_{{max}}$={langmuir_res['q_max']:.3f}, $K_L$={langmuir_res['K_L']:.3f}, $R^2$={langmuir_res['r2_nonlinear']:.4f})"
        )

    ax1.set_xlabel(r"Equilibrium Concentration $C_e$ (mg/L)", fontsize=11, fontweight="bold")
    ax1.set_ylabel(r"Adsorption Capacity $q_e$ (mg/g)", fontsize=11, fontweight="bold")
    ax1.set_title("Adsorption Isotherm Fitting ($q_e$ vs $C_e$)", fontsize=13, fontweight="bold")
    ax1.legend(loc="lower right", fontsize=9.5)
    ax1.grid(True, linestyle=":", alpha=0.6)

    # 2. Linearized Freundlich Plot: log10(q_e) vs log10(C_e)
    log_C = np.log10(freundlich_res["C_e"])
    log_q = np.log10(freundlich_res["q_e"])
    log_C_smooth = np.linspace(log_C.min() - 0.1, log_C.max() + 0.1, 100)
    log_q_smooth = freundlich_res["slope"] * log_C_smooth + freundlich_res["intercept"]

    ax2.scatter(log_C, log_q, color="darkmagenta", s=65, edgecolor="black", label=r"Data: $\log_{10}(q_e)$ vs $\log_{10}(C_e)$", zorder=5)
    ax2.plot(
        log_C_smooth, log_q_smooth, "m-", linewidth=2,
        label=rf"Linear Fit: Slope={freundlich_res['slope']:.3f}, $R^2$={freundlich_res['r2_linear']:.4f}"
    )
    ax2.set_xlabel(r"$\log_{10}(C_e)$", fontsize=11, fontweight="bold")
    ax2.set_ylabel(r"$\log_{10}(q_e)$", fontsize=11, fontweight="bold")
    ax2.set_title(r"Linearized Freundlich Plot ($\log_{10}(q_e) = \frac{1}{n}\log_{10}(C_e) + \log_{10}(K_f)$)", fontsize=12, fontweight="bold")
    ax2.legend(loc="upper left", fontsize=9.5)
    ax2.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    if save_filename:
        plt.savefig(save_filename, dpi=300)
        print(f"[Success] Plot saved to '{save_filename}'.")
    else:
        plt.show()


# ===============================================================================
# DEMONSTRATION & VERIFICATION RUN
# ===============================================================================
if __name__ == "__main__":
    # Experimental adsorption data
    # C_e : equilibrium concentration (mg/L)
    # q_e : amount adsorbed (mg/g)
    C_e_data = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0])
    q_e_data = np.array([2.50, 3.80, 4.90, 5.70, 6.40, 7.00, 7.50, 8.00])
    
    print("=" * 65)
    print("      ADSORPTION ISOTHERM ANALYSIS: FREUNDLICH & LANGMUIR")
    print("=" * 65)
    
    # 1. Fit Freundlich
    fr_res = freundlich_fit(C_e_data, q_e_data)
    print("\n[Freundlich Model]")
    print(f"  K_f (adsorption capacity) : {fr_res['K_f']:.6f} (mg/g)(L/mg)^(1/n)")
    print(f"  n (adsorption intensity)  : {fr_res['n']:.6f}")
    print(f"  1/n (heterogeneity factor): {fr_res['1/n']:.6f}")
    print(f"  Linear R² (log space)     : {fr_res['r2_linear']:.6f}")
    print(f"  Non-linear R² (real space): {fr_res['r2_nonlinear']:.6f}")
    print(f"  RMSE                      : {fr_res['rmse']:.6f} mg/g")
    
    # Adsorption feasibility interpretation
    print("\n  Adsorption Evaluation (Freundlich):")
    if fr_res['n'] > 1.0:
        print(f"  -> n = {fr_res['n']:.2f} > 1 : Favorable adsorption process.")
    elif fr_res['n'] == 1.0:
        print("  -> n = 1.0 : Linear partition adsorption.")
    else:
        print(f"  -> n = {fr_res['n']:.2f} < 1 : Unfavorable adsorption process.")

    # 2. Fit Langmuir
    lg_res = langmuir_fit(C_e_data, q_e_data)
    print("\n[Langmuir Model]")
    print(f"  q_max (monolayer capacity): {lg_res['q_max']:.6f} mg/g")
    print(f"  K_L (Langmuir constant)   : {lg_res['K_L']:.6f} L/mg")
    print(f"  Linear R²                 : {lg_res['r2_linear']:.6f}")
    print(f"  Non-linear R² (real space): {lg_res['r2_nonlinear']:.6f}")
    print(f"  RMSE                      : {lg_res['rmse']:.6f} mg/g")

    # Separation factor R_L at initial test concentration C0 = 80 mg/L
    C0_test = 80.0
    R_L = 1.0 / (1.0 + lg_res['K_L'] * C0_test)
    print(f"  Separation Factor (R_L @ C0={C0_test} mg/L): {R_L:.4f} (0 < R_L < 1 indicates favorable)")

    # Model comparison
    print("\n[Model Comparison]")
    better_model = "Freundlich" if fr_res['r2_nonlinear'] > lg_res['r2_nonlinear'] else "Langmuir"
    print(f"  Best fit based on Non-linear R²: {better_model} (Freundlich: {fr_res['r2_nonlinear']:.4f} vs Langmuir: {lg_res['r2_nonlinear']:.4f})")

    # Predictions
    print("\n[Predictions for New C_e Concentrations]")
    new_Ce = np.array([15.0, 25.0, 45.0, 90.0])
    pred_fr = predict_freundlich(new_Ce, fr_res['K_f'], fr_res['n'])
    pred_lg = predict_langmuir(new_Ce, lg_res['q_max'], lg_res['K_L'])
    print(f"  {'C_e (mg/L)':<12} {'Freundlich q_e (mg/g)':<24} {'Langmuir q_e (mg/g)':<24}")
    print("  " + "-" * 56)
    for c, q_f, q_l in zip(new_Ce, pred_fr, pred_lg):
        print(f"  {c:<12.1f} {q_f:<24.4f} {q_l:<24.4f}")

    print("=" * 65)

    # Plot (if matplotlib is available, will save or display)
    plot_isotherms(C_e_data, q_e_data, fr_res, lg_res, save_filename="isotherm_fit_plot.png")
