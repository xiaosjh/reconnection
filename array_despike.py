"""Python translation of the SolarSoft/IDL ``array_despike.pro`` routine.

The core algorithm iteratively replaces isolated high-count pixels with a
local weighted mean.  It can optionally repeat the same procedure for local
dark dips (the IDL ``NOLOW`` keyword disables that second stage).

Notes
-----
* The defaults follow the executable IDL source, not its older header comment:
  ``sigma=8``, ``threshold=6``, and ``itmax=20``.
* IDL ``CONVOL(..., /NORMALIZE, /EDGE_TRUNCATE)`` is reproduced with
  ``scipy.ndimage.convolve(..., mode='nearest') / kernel.sum()`` when every
  pixel is valid.  If ``valid_mask`` is supplied, invalid neighbours are
  excluded by a normalized convolution.
* The original IDL loop can execute ``itmax + 1`` passes because it stops when
  ``iter > itmax``.  This implementation preserves that behaviour.

使用示范
from array_despike import array_despike
#IRIS会将一些NaN设置为-200
valid_mask = (
    np.isfinite(spectrum)
    & (spectrum != -200)
)

spectrum_clean, info = array_despike(
    spectrum,
    sigma=8.0,
    threshold=6.0,
    itmax=20,
    no_low=True,#True是表示只去除尖峰，不处理异常暗点，GPT说对IRIS更安全，但我改成False也没啥区别
    valid_mask=valid_mask,
    return_info=True,
    verbose=True,
)
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import convolve


FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


def _default_kernel(ndim: int) -> FloatArray:
    """Return the default IDL kernel for a 1-D or 2-D input."""
    if ndim == 1:
        return np.array([1.0, 0.0, 1.0], dtype=float)
    if ndim == 2:
        return np.array(
            [[1.0, 1.0, 1.0],
             [1.0, 0.0, 1.0],
             [1.0, 1.0, 1.0]],
            dtype=float,
        )
    raise ValueError("array_despike supports only 1-D or 2-D arrays.")


def _validate_kernel(kernel: ArrayLike | None, ndim: int) -> FloatArray:
    """Validate a local-mean kernel."""
    if kernel is None:
        return _default_kernel(ndim)

    out = np.asarray(kernel, dtype=float)
    if out.ndim != ndim:
        raise ValueError(
            f"kernel must have {ndim} dimensions for an input with {ndim} dimensions."
        )
    if any(size % 2 == 0 for size in out.shape):
        raise ValueError("Every kernel dimension must have odd length.")
    if not np.all(np.isfinite(out)):
        raise ValueError("kernel must contain only finite values.")
    if np.any(out < 0):
        raise ValueError("kernel must be non-negative because it defines a local mean.")
    if out.sum() <= 0:
        raise ValueError("kernel must have a positive sum.")
    return out


def _local_mean(
    array: FloatArray,
    kernel: FloatArray,
    valid_mask: BoolArray,
    scale: FloatArray | None,
) -> FloatArray:
    """Compute the IDL-style local mean, optionally excluding invalid pixels."""
    if scale is None:
        normalized = array
    else:
        normalized = array / scale

    # With an all-True mask, the denominator is exactly kernel.sum(), including
    # at the boundary because mode='nearest' reproduces IDL /EDGE_TRUNCATE.
    values = np.where(valid_mask, normalized, 0.0)
    weights = valid_mask.astype(float)

    numerator = convolve(values, kernel, mode="nearest")
    denominator = convolve(weights, kernel, mode="nearest")

    mean = np.full(array.shape, np.nan, dtype=float)
    np.divide(numerator, denominator, out=mean, where=denominator > 0)

    if scale is not None:
        mean *= scale
    return mean


def array_despike(
    array: ArrayLike,
    *,
    kernel: ArrayLike | None = None,
    sigma: float = 8.0,
    threshold: float = 6.0,
    itmax: int = 20,
    scale: ArrayLike | None = None,
    no_low: bool = False,
    valid_mask: ArrayLike | None = None,
    verbose: bool = False,
    return_info: bool = False,
) -> FloatArray | tuple[FloatArray, dict[str, Any]]:
    """Remove isolated spikes from a 1-D or 2-D count array.

    This is a Python translation of ``array_despike.pro``.

    Parameters
    ----------
    array
        Input count array.  For IRIS spectra this should normally be the
        Level-2 data in DN, before division by exposure time.
    kernel
        Local-mean kernel.  The default is ``[1, 0, 1]`` for 1-D input or the
        eight-neighbour 3x3 kernel for 2-D input.
    sigma
        Required excess above the local Poisson scale.  The executable IDL
        source uses 8.0 by default.
    threshold
        Minimum count value for a bright point to be replaced.  In the dark
        stage, the local mean must exceed this value.
    itmax
        IDL iteration limit.  To preserve the original ``iter > itmax`` test,
        at most ``itmax + 1`` passes are made in each stage.
    scale
        Optional array with the same shape as ``array``.  The local mean is
        then computed from ``array / scale`` and multiplied by the centre
        pixel's scale, matching the IDL routine.  Use only positive finite
        values.
    no_low
        If True, remove only bright spikes.  This corresponds to the IDL
        ``/NOLOW`` keyword and is generally the safer choice for IRIS spectra.
    valid_mask
        Boolean mask with True for valid pixels.  NaN and Inf are always
        treated as invalid.  For IRIS, use this to exclude explicit fill
        values such as -200 without clipping ordinary negative noise.
    verbose
        Print the number of replacements in each pass.
    return_info
        If True, also return masks and per-pass replacement counts.

    Returns
    -------
    cleaned
        Floating-point array with detected spikes replaced by the local mean.
    info
        Returned only when ``return_info=True``.  It contains the union mask,
        separate bright/dark masks, and replacement counts.

    Notes
    -----
    The bright-spike criterion is

    ``(array - local_mean) / max(sqrt(local_mean), 1) > sigma``

    together with ``array >= threshold``.  The optional dark-dip stage uses

    ``(local_mean - array) / max(sqrt(array), 1) > sigma``

    together with ``local_mean >= threshold``.
    """
    data = np.asarray(array, dtype=float)
    if data.ndim not in (1, 2):
        raise ValueError(
            f"array must be 1-D or 2-D; received shape {data.shape}. "
            "Apply the function frame by frame to a data cube."
        )
    if sigma <= 0 or not np.isfinite(sigma):
        raise ValueError("sigma must be a positive finite number.")
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite.")
    if not isinstance(itmax, (int, np.integer)) or itmax < 0:
        raise ValueError("itmax must be a non-negative integer.")

    local_kernel = _validate_kernel(kernel, data.ndim)
    cleaned = data.copy()

    finite = np.isfinite(cleaned)
    if valid_mask is None:
        valid = finite
    else:
        supplied_mask = np.asarray(valid_mask, dtype=bool)
        if supplied_mask.shape != cleaned.shape:
            raise ValueError("valid_mask must have the same shape as array.")
        valid = supplied_mask & finite

    scale_array: FloatArray | None
    if scale is None:
        scale_array = None
    else:
        scale_array = np.asarray(scale, dtype=float)
        if scale_array.ndim == 0:
            # A scalar scale has no effect in the original IDL implementation.
            scale_array = None
        else:
            if scale_array.shape != cleaned.shape:
                raise ValueError("scale must be scalar or have the same shape as array.")
            if np.any(~np.isfinite(scale_array[valid])) or np.any(scale_array[valid] <= 0):
                raise ValueError("scale must be positive and finite at every valid pixel.")
            valid &= np.isfinite(scale_array) & (scale_array > 0)

    bright_mask = np.zeros(cleaned.shape, dtype=bool)
    dark_mask = np.zeros(cleaned.shape, dtype=bool)
    bright_counts: list[int] = []
    dark_counts: list[int] = []

    # The original IDL code stops after `iter > itmax`, hence itmax + 1 passes.
    for iteration in range(itmax + 1):
        local = _local_mean(cleaned, local_kernel, valid, scale_array)

        # denominator = np.full(cleaned.shape, np.nan, dtype=float)
        # nonnegative_local = np.isfinite(local) & (local >= 0)
        # denominator[nonnegative_local] = np.maximum(
        #     np.sqrt(local[nonnegative_local]), 1.0
        # )

        # significance = np.full(cleaned.shape, np.nan, dtype=float)
        # np.divide(
        #     cleaned - local,
        #     denominator,
        #     out=significance,
        #     where=np.isfinite(denominator),
        # )
        denominator = np.maximum(
            np.sqrt(np.maximum(local, 0.0)),
            1.0
        )

        significance = (cleaned - local) / denominator

        replace = (
            valid
            & np.isfinite(local)
            & (significance > sigma)
            & (cleaned >= threshold)
        )
        count = int(np.count_nonzero(replace))
        bright_counts.append(count)

        if count:
            cleaned[replace] = local[replace]
            bright_mask |= replace

        if verbose:
            print(f"bright pass {iteration}: replaced {count} pixels")
        if count == 0:
            break

    if not no_low:
        for iteration in range(itmax + 1):
            local = _local_mean(cleaned, local_kernel, valid, scale_array)

            denominator = np.full(cleaned.shape, np.nan, dtype=float)
            nonnegative_data = valid & (cleaned >= 0)
            denominator[nonnegative_data] = np.maximum(
                np.sqrt(cleaned[nonnegative_data]), 1.0
            )

            significance = np.full(cleaned.shape, np.nan, dtype=float)
            np.divide(
                local - cleaned,
                denominator,
                out=significance,
                where=np.isfinite(denominator),
            )

            replace = (
                valid
                & np.isfinite(local)
                & (significance > sigma)
                & (local >= threshold)
            )
            count = int(np.count_nonzero(replace))
            dark_counts.append(count)

            if count:
                cleaned[replace] = local[replace]
                dark_mask |= replace

            if verbose:
                print(f"dark pass {iteration}: replaced {count} pixels")
            if count == 0:
                break

    # Preserve invalid input samples exactly rather than replacing them.
    cleaned[~valid] = data[~valid]

    if not return_info:
        return cleaned

    changed_mask = bright_mask | dark_mask
    info: dict[str, Any] = {
        "mask": changed_mask,
        "bright_mask": bright_mask,
        "dark_mask": dark_mask,
        "bright_counts": tuple(bright_counts),
        "dark_counts": tuple(dark_counts),
        "n_changed_unique": int(np.count_nonzero(changed_mask)),
        "n_bright_unique": int(np.count_nonzero(bright_mask)),
        "n_dark_unique": int(np.count_nonzero(dark_mask)),
        "n_replacements_total": int(sum(bright_counts) + sum(dark_counts)),
    }
    return cleaned, info


def despike_stack(
    stack: ArrayLike,
    *,
    stack_axis: int = 0,
    return_info: bool = False,
    **kwargs: Any,
) -> FloatArray | tuple[FloatArray, list[dict[str, Any]]]:
    """Apply :func:`array_despike` independently to every 2-D frame in a 3-D cube.

    Parameters
    ----------
    stack
        Three-dimensional array.
    stack_axis
        Axis indexing independent frames, commonly the time/raster axis.
    return_info
        Return one information dictionary per frame.
    **kwargs
        Forwarded to :func:`array_despike`.
    """
    cube = np.asarray(stack, dtype=float)
    if cube.ndim != 3:
        raise ValueError("stack must be a 3-D array.")

    moved = np.moveaxis(cube, stack_axis, 0)
    cleaned = np.empty_like(moved, dtype=float)
    all_info: list[dict[str, Any]] = []

    # A cube-wide valid_mask may be supplied and is sliced consistently.
    cube_mask = kwargs.pop("valid_mask", None)
    if cube_mask is not None:
        cube_mask = np.asarray(cube_mask, dtype=bool)
        if cube_mask.shape != cube.shape:
            raise ValueError("A cube-wide valid_mask must match stack.shape.")
        moved_mask = np.moveaxis(cube_mask, stack_axis, 0)
    else:
        moved_mask = None

    cube_scale = kwargs.pop("scale", None)
    if cube_scale is not None and np.asarray(cube_scale).ndim > 0:
        cube_scale = np.asarray(cube_scale, dtype=float)
        if cube_scale.shape != cube.shape:
            raise ValueError("A cube-wide scale array must match stack.shape.")
        moved_scale = np.moveaxis(cube_scale, stack_axis, 0)
    else:
        moved_scale = cube_scale

    for index, frame in enumerate(moved):
        frame_kwargs = dict(kwargs)
        if moved_mask is not None:
            frame_kwargs["valid_mask"] = moved_mask[index]
        if isinstance(moved_scale, np.ndarray):
            frame_kwargs["scale"] = moved_scale[index]
        elif moved_scale is not None:
            frame_kwargs["scale"] = moved_scale

        result = array_despike(frame, return_info=return_info, **frame_kwargs)
        if return_info:
            cleaned[index], frame_info = result
            all_info.append(frame_info)
        else:
            cleaned[index] = result

    cleaned = np.moveaxis(cleaned, 0, stack_axis)
    if return_info:
        return cleaned, all_info
    return cleaned
