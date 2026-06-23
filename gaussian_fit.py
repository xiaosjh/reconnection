"""Small Gaussian fitting helper inspired by IDL MPFITPEAK.

The main entry point is fit_gaussian(). It fits a single Gaussian peak with
either a constant or linear background and returns a dictionary for convenient
use in notebooks.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.optimize import curve_fit


def _gaussian_constant_background(
    x: np.ndarray,
    amplitude: float,
    center: float,
    sigma: float,
    background: float,
) -> np.ndarray:
    width = max(abs(float(sigma)), np.finfo(float).tiny)
    return background + amplitude * np.exp(-0.5 * ((x - center) / width) ** 2)


def _gaussian_linear_background(
    x: np.ndarray,
    amplitude: float,
    center: float,
    sigma: float,
    background: float,
    slope: float,
    *,
    x_ref: float,
) -> np.ndarray:
    width = max(abs(float(sigma)), np.finfo(float).tiny)
    return background + slope * (x - x_ref) + amplitude * np.exp(
        -0.5 * ((x - center) / width) ** 2
    )


def _double_gaussian_constant_background(
    x: np.ndarray,
    amplitude1: float,
    center1: float,
    sigma1: float,
    amplitude2: float,
    center2: float,
    sigma2: float,
    background: float,
) -> np.ndarray:
    width1 = max(abs(float(sigma1)), np.finfo(float).tiny)
    width2 = max(abs(float(sigma2)), np.finfo(float).tiny)
    return (
        background
        + amplitude1 * np.exp(-0.5 * ((x - center1) / width1) ** 2)
        + amplitude2 * np.exp(-0.5 * ((x - center2) / width2) ** 2)
    )


def _double_gaussian_linear_background(
    x: np.ndarray,
    amplitude1: float,
    center1: float,
    sigma1: float,
    amplitude2: float,
    center2: float,
    sigma2: float,
    background: float,
    slope: float,
    *,
    x_ref: float,
) -> np.ndarray:
    width1 = max(abs(float(sigma1)), np.finfo(float).tiny)
    width2 = max(abs(float(sigma2)), np.finfo(float).tiny)
    return (
        background
        + slope * (x - x_ref)
        + amplitude1 * np.exp(-0.5 * ((x - center1) / width1) ** 2)
        + amplitude2 * np.exp(-0.5 * ((x - center2) / width2) ** 2)
    )


def _estimate_initial_parameters(
    x: np.ndarray,
    y: np.ndarray,
    *,
    positive: bool | None,
    linear_background: bool,
    x_ref: float,
) -> np.ndarray:
    order = np.argsort(x)
    xs = x[order]
    ys = y[order]

    n_edge = max(1, min(xs.size // 5, 20))
    edge_y = np.concatenate([ys[:n_edge], ys[-n_edge:]])
    background = float(np.nanmedian(edge_y))

    if positive is True:
        peak_index = int(np.nanargmax(y - background))
    elif positive is False:
        peak_index = int(np.nanargmin(y - background))
    else:
        high_index = int(np.nanargmax(y))
        low_index = int(np.nanargmin(y))
        high_distance = abs(float(y[high_index] - background))
        low_distance = abs(float(y[low_index] - background))
        peak_index = high_index if high_distance >= low_distance else low_index

    amplitude = float(y[peak_index] - background)
    center = float(x[peak_index])

    half_level = background + 0.5 * amplitude
    if amplitude >= 0:
        near_peak = y >= half_level
    else:
        near_peak = y <= half_level

    if np.count_nonzero(near_peak) >= 2:
        width_guess = 0.5 * (np.nanmax(x[near_peak]) - np.nanmin(x[near_peak]))
        sigma = float(width_guess / np.sqrt(2.0 * np.log(2.0)))
    else:
        sigma = float((np.nanmax(x) - np.nanmin(x)) / 6.0)

    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.nanmedian(np.abs(np.diff(xs))))
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = 1.0

    if linear_background:
        left_x = float(np.nanmean(xs[:n_edge]))
        right_x = float(np.nanmean(xs[-n_edge:]))
        left_y = float(np.nanmedian(ys[:n_edge]))
        right_y = float(np.nanmedian(ys[-n_edge:]))
        dx = right_x - left_x
        slope = (right_y - left_y) / dx if dx != 0 else 0.0
        background = float(background - slope * (float(np.nanmedian(xs)) - x_ref))
        return np.array([amplitude, center, sigma, background, slope], dtype=float)

    return np.array([amplitude, center, sigma, background], dtype=float)


def _estimate_double_initial_parameters(
    x: np.ndarray,
    y: np.ndarray,
    *,
    positive: bool | None,
    linear_background: bool,
    x_ref: float,
) -> np.ndarray:
    first = _estimate_initial_parameters(
        x,
        y,
        positive=positive,
        linear_background=False,
        x_ref=x_ref,
    )
    first_model = _gaussian_constant_background(x, *first)
    residual = y - first_model + first[3]
    second = _estimate_initial_parameters(
        x,
        residual,
        positive=positive,
        linear_background=False,
        x_ref=x_ref,
    )

    background = first[3]
    params = np.array(
        [
            first[0],
            first[1],
            first[2],
            second[0],
            second[1],
            second[2],
            background,
        ],
        dtype=float,
    )

    order = np.argsort(x)
    xs = x[order]
    ys = y[order]
    n_edge = max(1, min(xs.size // 5, 20))

    if linear_background:
        left_x = float(np.nanmean(xs[:n_edge]))
        right_x = float(np.nanmean(xs[-n_edge:]))
        left_y = float(np.nanmedian(ys[:n_edge]))
        right_y = float(np.nanmedian(ys[-n_edge:]))
        dx = right_x - left_x
        slope = (right_y - left_y) / dx if dx != 0 else 0.0
        params = np.append(params, slope)

    return params


def _clean_inputs(
    x: np.ndarray,
    y: np.ndarray,
    yerr: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()

    if x.size != y.size:
        raise ValueError("x and y must have the same number of elements.")
    if x.size < 4:
        raise ValueError("At least four data points are required.")

    mask = np.isfinite(x) & np.isfinite(y)

    if yerr is not None:
        yerr = np.asarray(yerr, dtype=float).ravel()
        if yerr.size != y.size:
            raise ValueError("yerr must have the same number of elements as y.")
        mask &= np.isfinite(yerr) & (yerr > 0)
        yerr = yerr[mask]

    x = x[mask]
    y = y[mask]

    if x.size < 4:
        raise ValueError("At least four finite data points are required.")

    return x, y, yerr


def fit_gaussian(
    x: np.ndarray,
    y: np.ndarray,
    yerr: np.ndarray | None = None,
    p0: list[float] | tuple[float, ...] | np.ndarray | None = None,
    *,
    amplitude_guess: float | None = None,
    center_guess: float | None = None,
    sigma_guess: float | None = None,
    background_guess: float | None = None,
    slope_guess: float | None = None,
    positive: bool | None = None,
    linear_background: bool = False,
    maxfev: int = 10000,
) -> dict[str, Any]:
    """Fit a single Gaussian peak and return a result dictionary.

    Parameters
    ----------
    x, y:
        One-dimensional data arrays.
    yerr:
        Optional 1-sigma uncertainty for y. If supplied, points with non-finite
        or non-positive uncertainty are ignored.
    p0:
        Optional initial guess. Use [amplitude, center, sigma, background] for
        a constant background, or add slope as the fifth value for a linear
        background.
    amplitude_guess, center_guess, sigma_guess, background_guess, slope_guess:
        Optional per-parameter initial guesses. Any value left as None is
        estimated from the data. Do not combine these with p0.
    positive:
        True for an emission peak, False for an absorption peak, None to infer.
    linear_background:
        If True, fit background + slope * (x - x_ref) instead of a constant
        background. The returned slope is in y-units per x-unit.
    maxfev:
        Maximum number of function evaluations for scipy.optimize.curve_fit.
    """
    x_clean, y_clean, yerr_clean = _clean_inputs(x, y, yerr)
    x_ref = float(np.nanmean(x_clean))

    n_params = 5 if linear_background else 4
    individual_guesses = [
        amplitude_guess,
        center_guess,
        sigma_guess,
        background_guess,
        slope_guess,
    ]
    if p0 is None:
        p0_array = _estimate_initial_parameters(
            x_clean,
            y_clean,
            positive=positive,
            linear_background=linear_background,
            x_ref=x_ref,
        )
        guess_overrides = individual_guesses[:n_params]
        for index, value in enumerate(guess_overrides):
            if value is not None:
                p0_array[index] = float(value)
    else:
        if any(value is not None for value in individual_guesses):
            raise ValueError("Use either p0 or individual *_guess values, not both.")
        p0_array = np.asarray(p0, dtype=float).ravel()
        if p0_array.size != n_params:
            raise ValueError(f"p0 must contain {n_params} values.")

    lower = np.full(n_params, -np.inf)
    upper = np.full(n_params, np.inf)
    lower[2] = np.finfo(float).tiny
    if positive is True:
        lower[0] = 0.0
    elif positive is False:
        upper[0] = 0.0

    if p0_array[2] <= 0:
        p0_array[2] = abs(p0_array[2]) if p0_array[2] != 0 else 1.0

    if linear_background:
        model = lambda x_model, amp, cen, sig, bg, slp: _gaussian_linear_background(
            x_model, amp, cen, sig, bg, slp, x_ref=x_ref
        )
    else:
        model = _gaussian_constant_background

    popt, pcov = curve_fit(
        model,
        x_clean,
        y_clean,
        p0=p0_array,
        sigma=yerr_clean,
        absolute_sigma=yerr_clean is not None,
        bounds=(lower, upper),
        maxfev=maxfev,
    )

    popt = np.asarray(popt, dtype=float)
    popt[2] = abs(popt[2])
    yfit = model(x_clean, *popt)
    residuals = y_clean - yfit
    dof = max(0, x_clean.size - popt.size)

    if yerr_clean is not None:
        chi2 = float(np.sum((residuals / yerr_clean) ** 2))
    else:
        chi2 = float(np.sum(residuals**2))

    if pcov is None or not np.all(np.isfinite(pcov)):
        param_errors = np.full_like(popt, np.nan, dtype=float)
    else:
        param_errors = np.sqrt(np.diag(pcov))

    result: dict[str, Any] = {
        "success": True,
        "message": "Fit converged.",
        "amplitude": float(popt[0]),
        "center": float(popt[1]),
        "sigma": float(popt[2]),
        "background": float(popt[3]),
        "slope": float(popt[4]) if linear_background else 0.0,
        "x_ref": x_ref,
        "area": float(popt[0] * popt[2] * np.sqrt(2.0 * np.pi)),
        "params": popt,
        "param_errors": param_errors,
        "covariance": pcov,
        "x": x_clean,
        "y": y_clean,
        "yfit": yfit,
        "residuals": residuals,
        "chi2": chi2,
        "dof": dof,
        "linear_background": linear_background,
    }
    return result


def fit_double_gaussian(
    x: np.ndarray,
    y: np.ndarray,
    yerr: np.ndarray | None = None,
    p0: list[float] | tuple[float, ...] | np.ndarray | None = None,
    *,
    amplitude1_guess: float | None = None,
    center1_guess: float | None = None,
    sigma1_guess: float | None = None,
    amplitude2_guess: float | None = None,
    center2_guess: float | None = None,
    sigma2_guess: float | None = None,
    background_guess: float | None = None,
    slope_guess: float | None = None,
    positive: bool | None = True,
    linear_background: bool = False,
    maxfev: int = 20000,
) -> dict[str, Any]:
    """Fit two Gaussian components with a shared background.

    Parameters
    ----------
    x, y:
        One-dimensional data arrays.
    yerr:
        Optional 1-sigma uncertainty for y.
    p0:
        Optional initial guess. Use
        [amp1, center1, sigma1, amp2, center2, sigma2, background] for a
        constant background, or add slope as the eighth value for a linear
        background.
    amplitude1_guess, center1_guess, sigma1_guess:
        Optional initial guesses for the first Gaussian component.
    amplitude2_guess, center2_guess, sigma2_guess:
        Optional initial guesses for the second Gaussian component.
    background_guess, slope_guess:
        Optional initial guesses for the shared background. slope_guess is
        used only when linear_background=True.
    positive:
        True constrains both components to emission peaks, False constrains
        both to absorption peaks, and None leaves amplitudes unconstrained.
    linear_background:
        If True, fit background + slope * (x - x_ref) plus two Gaussians.
    maxfev:
        Maximum number of function evaluations for scipy.optimize.curve_fit.
    """
    x_clean, y_clean, yerr_clean = _clean_inputs(x, y, yerr)
    x_ref = float(np.nanmean(x_clean))

    n_params = 8 if linear_background else 7
    individual_guesses = [
        amplitude1_guess,
        center1_guess,
        sigma1_guess,
        amplitude2_guess,
        center2_guess,
        sigma2_guess,
        background_guess,
        slope_guess,
    ]
    if p0 is None:
        p0_array = _estimate_double_initial_parameters(
            x_clean,
            y_clean,
            positive=positive,
            linear_background=linear_background,
            x_ref=x_ref,
        )
        guess_overrides = individual_guesses[:n_params]
        for index, value in enumerate(guess_overrides):
            if value is not None:
                p0_array[index] = float(value)
    else:
        if any(value is not None for value in individual_guesses):
            raise ValueError("Use either p0 or individual *_guess values, not both.")
        p0_array = np.asarray(p0, dtype=float).ravel()
        if p0_array.size != n_params:
            raise ValueError(f"p0 must contain {n_params} values.")

    lower = np.full(n_params, -np.inf)
    upper = np.full(n_params, np.inf)
    lower[2] = np.finfo(float).tiny
    lower[5] = np.finfo(float).tiny

    if positive is True:
        lower[0] = 0.0
        lower[3] = 0.0
    elif positive is False:
        upper[0] = 0.0
        upper[3] = 0.0

    for sigma_index in (2, 5):
        if p0_array[sigma_index] <= 0:
            p0_array[sigma_index] = (
                abs(p0_array[sigma_index]) if p0_array[sigma_index] != 0 else 1.0
            )

    if linear_background:
        model = (
            lambda x_model, amp1, cen1, sig1, amp2, cen2, sig2, bg, slp: (
                _double_gaussian_linear_background(
                    x_model,
                    amp1,
                    cen1,
                    sig1,
                    amp2,
                    cen2,
                    sig2,
                    bg,
                    slp,
                    x_ref=x_ref,
                )
            )
        )
    else:
        model = _double_gaussian_constant_background

    popt, pcov = curve_fit(
        model,
        x_clean,
        y_clean,
        p0=p0_array,
        sigma=yerr_clean,
        absolute_sigma=yerr_clean is not None,
        bounds=(lower, upper),
        maxfev=maxfev,
    )

    popt = np.asarray(popt, dtype=float)
    popt[2] = abs(popt[2])
    popt[5] = abs(popt[5])

    baseline = popt[6]
    slope = float(popt[7]) if linear_background else 0.0
    background_model = baseline + slope * (x_clean - x_ref)
    component1 = popt[0] * np.exp(-0.5 * ((x_clean - popt[1]) / popt[2]) ** 2)
    component2 = popt[3] * np.exp(-0.5 * ((x_clean - popt[4]) / popt[5]) ** 2)
    yfit = background_model + component1 + component2
    residuals = y_clean - yfit
    dof = max(0, x_clean.size - popt.size)

    if yerr_clean is not None:
        chi2 = float(np.sum((residuals / yerr_clean) ** 2))
    else:
        chi2 = float(np.sum(residuals**2))

    if pcov is None or not np.all(np.isfinite(pcov)):
        param_errors = np.full_like(popt, np.nan, dtype=float)
    else:
        param_errors = np.sqrt(np.diag(pcov))

    area1 = float(popt[0] * popt[2] * np.sqrt(2.0 * np.pi))
    area2 = float(popt[3] * popt[5] * np.sqrt(2.0 * np.pi))

    result: dict[str, Any] = {
        "success": True,
        "message": "Fit converged.",
        "amplitude1": float(popt[0]),
        "center1": float(popt[1]),
        "sigma1": float(popt[2]),
        "area1": area1,
        "amplitude2": float(popt[3]),
        "center2": float(popt[4]),
        "sigma2": float(popt[5]),
        "area2": area2,
        "background": float(popt[6]),
        "slope": slope,
        "x_ref": x_ref,
        "total_area": float(area1 + area2),
        "component1": component1,
        "component2": component2,
        "background_model": background_model,
        "params": popt,
        "param_errors": param_errors,
        "covariance": pcov,
        "x": x_clean,
        "y": y_clean,
        "yfit": yfit,
        "residuals": residuals,
        "chi2": chi2,
        "dof": dof,
        "linear_background": linear_background,
    }
    return result
