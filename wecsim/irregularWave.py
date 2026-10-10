"""Directional irregular-wave excitation from WEC-Sim hydrodynamic data.

The formulas follow the current MATLAB WEC-Sim ``irregWaveSpectrum``,
``waveElevIrreg``, and ``irregExcF`` path for a PM spectrum. Callers may supply
the phase matrix to replay a MATLAB realization or generate a reproducible
Python realization with an integer seed. ``phase_generator="matlab"`` selects
the pinned WEC-Sim Threefry substream for independent seeded replay.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from randomgen import ThreeFry
from scipy.interpolate import CubicSpline, RegularGridInterpolator
from scipy.io import loadmat

from .bodyClass import BodyClass


@dataclass(frozen=True)
class IrregularComponents:
    omega: np.ndarray
    spectral_amplitude: np.ndarray
    d_omega: np.ndarray
    directions: np.ndarray
    spreading: np.ndarray
    phase: np.ndarray


@dataclass(frozen=True)
class FullDirectionalComponents:
    """Frequency-by-heading imported spectrum with angular bin widths."""

    omega: np.ndarray
    spectral_amplitude: np.ndarray  # (frequency, heading), m² s/rad²
    d_omega: np.ndarray
    directions: np.ndarray  # degrees
    d_theta: np.ndarray  # radians
    phase: np.ndarray  # (frequency, heading)


@dataclass(frozen=True)
class IrregularResponse:
    time: np.ndarray
    elevation: np.ndarray
    excitation_force: np.ndarray


def _random_phases(shape: tuple[int, int], seed: int | None,
                   phase_generator: str) -> np.ndarray:
    if phase_generator == "numpy":
        return 2 * np.pi * np.random.default_rng(seed).random(shape)
    if phase_generator != "matlab":
        raise ValueError("phase_generator must be 'numpy' or 'matlab'")
    if (not isinstance(seed, int) or isinstance(seed, bool)
            or not 1 <= seed < 2**32 - 4):
        raise ValueError("MATLAB phase generation needs a positive integer substream seed")
    # Pinned waveClass: RandStream('Threefry','Seed',1), Substream=phaseSeed.
    # MATLAB's 17-word State was paired with R2025b for substreams 1-3 and
    # the published 500-by-3 Morison sea (substream 5). randomgen starts at
    # the matching block when these four uint64 counter words are supplied.
    counter = np.array([2 << 32, (4 << 32) + 3,
                        (6 << 32) + seed + 4, (8 << 32) + 7], dtype=np.uint64)
    return 2 * np.pi * np.random.Generator(
        ThreeFry(key=0, counter=counter)).random(shape)


def imported_full_directional_components(
    h5_file: str | Path,
    spectrum_file: str | Path,
    *,
    phase: np.ndarray | None = None,
    seed: int | None = None,
    phase_generator: str = "numpy",
) -> FullDirectionalComponents:
    """Load WEC-Sim's frequency-resolved directional MAT spectrum.

    Imported ``spread`` is a density per heading radian, so its row integral
    uses ``d_theta``. A supplied phase matrix replays a MATLAB realization;
    ``seed`` makes an independent reproducible Python realization.
    """
    source = loadmat(spectrum_file)
    required = ("frequencies", "spectrum", "spread", "directions")
    if any(name not in source for name in required):
        raise ValueError("full-directional MAT file needs frequencies, spectrum, spread, and directions")
    frequency = np.asarray(source["frequencies"], dtype=float).ravel()
    density = np.asarray(source["spectrum"], dtype=float).ravel()
    heading = np.asarray(source["directions"], dtype=float).ravel()
    spread = np.asarray(source["spread"], dtype=float)
    if (len(frequency) < 2 or len(heading) < 2
            or spread.shape != (len(frequency), len(heading))
            or not all(np.isfinite(values).all()
                       for values in (frequency, density, heading, spread))
            or density.shape != frequency.shape
            or np.any(density < 0) or np.any(spread < 0)
            or np.any(np.diff(frequency) <= 0)
            or np.any(np.diff(heading) <= 0)):
        raise ValueError("full-directional spectrum has invalid frequency, direction, or density")
    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.readH5file()
    bem_omega = np.asarray(body.hydroData["simulation_parameters"]["w"]).ravel()
    selected = ((frequency * 2 * np.pi >= bem_omega.min())
                & (frequency * 2 * np.pi <= bem_omega.max()))
    omega = frequency[selected] * (2 * np.pi)
    if len(omega) < 2:
        raise ValueError("full-directional spectrum needs two BEM-range frequencies")
    width = np.empty(len(omega))
    width[0] = omega[1] - omega[0]
    width[-1] = omega[-1] - omega[-2]
    width[1:-1] = (omega[2:] - omega[:-2]) / 2
    direction_radians = np.deg2rad(heading)
    angle_width = np.empty(len(heading))
    angle_width[0] = direction_radians[1] - direction_radians[0]
    angle_width[-1] = direction_radians[-1] - direction_radians[-2]
    angle_width[1:-1] = (direction_radians[2:] - direction_radians[:-2]) / 2
    if (np.any(width <= 0) or np.any(angle_width <= 0)
            or not np.allclose(spread[selected] @ angle_width, 1,
                               rtol=0, atol=1e-6)):
        raise ValueError("directional spread must integrate to one at each frequency")
    if phase is not None and seed is not None:
        raise ValueError("supply either phase or seed")
    shape = (len(omega), len(heading))
    if phase is None:
        phases = _random_phases(shape, seed, phase_generator)
    else:
        phases = np.asarray(phase, dtype=float)
        if phases.shape != shape or not np.isfinite(phases).all():
            raise ValueError("phase must match the selected frequency-by-heading grid")
    return FullDirectionalComponents(
        omega=omega,
        spectral_amplitude=density[selected, None] * spread[selected] / np.pi,
        d_omega=width,
        directions=heading,
        d_theta=angle_width,
        phase=phases,
    )


def imported_spectrum_components(
    h5_file: str | Path,
    spectrum_file: str | Path,
) -> IrregularComponents:
    """Read a three-column WEC-Sim spectrumImport MAT file with saved phases.

    The columns are frequency in Hz, spectral density in m²/Hz, and phase in
    radians. The current MATLAB wave class keeps only BEM-range frequencies,
    computes midpoint bin widths in rad/s, and converts density to m²/(rad/s).
    Its phase column makes the realized sea state deterministic.
    """
    source = loadmat(spectrum_file)
    if "spectrumData" not in source:
        raise ValueError("imported spectrum MAT file needs spectrumData")
    values = np.asarray(source["spectrumData"], dtype=float)
    if (values.ndim != 2 or values.shape[1] != 3
            or not np.isfinite(values).all()):
        raise ValueError("spectrumData must have three finite columns")

    body = BodyClass(str(h5_file))
    body.bodyNumber = 1
    body.readH5file()
    bem_omega = np.asarray(body.hydroData["simulation_parameters"]["w"]).ravel()
    if len(bem_omega) < 2 or not np.isfinite(bem_omega).all():
        raise ValueError("hydrodynamic frequency range is invalid")
    frequency = values[:, 0]
    keep = ((frequency >= bem_omega.min() / (2 * np.pi))
            & (frequency <= bem_omega.max() / (2 * np.pi)))
    selected = values[keep]
    if (len(selected) < 2 or np.any(np.diff(selected[:, 0]) <= 0)
            or np.any(selected[:, 1] < 0)):
        raise ValueError("imported BEM-range spectrum needs increasing frequencies and nonnegative density")
    omega = selected[:, 0] * (2 * np.pi)
    width = np.empty(len(omega))
    width[0] = omega[1] - omega[0]
    width[-1] = omega[-1] - omega[-2]
    width[1:-1] = (omega[2:] - omega[:-2]) / 2
    return IrregularComponents(
        omega=omega,
        spectral_amplitude=selected[:, 1] / np.pi,
        d_omega=width,
        directions=np.array([0.0]),
        spreading=np.array([1.0]),
        phase=selected[:, 2, None],
    )


def pm_equal_energy_components(
    h5_file: str | Path | None,
    *,
    significant_height: float,
    peak_period: float,
    directions: np.ndarray,
    spreading: np.ndarray,
    count: int | None = None,
    seed: int | None = None,
    phase_generator: str = "numpy",
    phase: np.ndarray | None = None,
    frequency_range: Sequence[float] | None = None,
    discretization: str = "equal_energy",
) -> IrregularComponents:
    """Build current WEC-Sim PM bins from a frequency range.

    ``phase`` overrides random generation for replay. The default seed uses
    NumPy; ``phase_generator="matlab"`` uses the pinned WEC-Sim Threefry
    substream. ``frequency_range`` supplies the WEC-Sim ``bem.range``
    used when a Morison-only body has no HDF5 file.
    """
    return _equal_energy_components(
        h5_file, significant_height=significant_height,
        peak_period=peak_period, directions=directions, spreading=spreading,
        count=count, seed=seed, phase=phase, gamma=None,
        phase_generator=phase_generator,
        frequency_range=frequency_range,
        discretization=discretization,
    )


def jonswap_equal_energy_components(
    h5_file: str | Path,
    *,
    significant_height: float,
    peak_period: float,
    directions: np.ndarray,
    spreading: np.ndarray,
    count: int | None = None,
    seed: int | None = None,
    phase_generator: str = "numpy",
    phase: np.ndarray | None = None,
    gamma: float | None = None,
    frequency_range: Sequence[float] | None = None,
    discretization: str = "equal_energy",
) -> IrregularComponents:
    """Build WEC-Sim's IEC JONSWAP bins.

    ``gamma=None`` uses WEC-Sim's height/period-dependent peak factor.
    Supplied phases replay a MATLAB sea realization.
    """
    if (not np.isfinite([significant_height, peak_period]).all()
            or significant_height <= 0 or peak_period <= 0):
        raise ValueError("significant_height and peak_period must be positive and finite")
    if gamma is None:
        ratio = peak_period / np.sqrt(significant_height)
        gamma = (5.0 if ratio <= 3.6 else 1.0 if ratio > 5.0
                 else np.exp(5.75 - 1.15 * ratio))
    if not np.isfinite(gamma) or gamma <= 0:
        raise ValueError("JONSWAP gamma must be positive and finite")
    return _equal_energy_components(
        h5_file, significant_height=significant_height,
        peak_period=peak_period, directions=directions, spreading=spreading,
        count=count, seed=seed, phase=phase, gamma=gamma,
        phase_generator=phase_generator,
        frequency_range=frequency_range,
        discretization=discretization,
    )


def _equal_energy_components(
    h5_file: str | Path | None, *, significant_height: float, peak_period: float,
    directions: np.ndarray, spreading: np.ndarray, count: int | None,
    seed: int | None, phase: np.ndarray | None, gamma: float | None,
    phase_generator: str = "numpy",
    frequency_range: Sequence[float] | None = None,
    discretization: str = "equal_energy",
) -> IrregularComponents:
    if not isinstance(discretization, str):
        raise ValueError("discretization must be 'equal_energy' or 'traditional'")
    scheme = discretization.lower().replace("_", "")
    if scheme not in ("equalenergy", "traditional"):
        raise ValueError("discretization must be 'equal_energy' or 'traditional'")
    if count is None:
        count = 1000 if scheme == "traditional" else 500
    if (not np.isfinite([significant_height, peak_period]).all()
            or significant_height <= 0 or peak_period <= 0):
        raise ValueError("significant_height and peak_period must be positive and finite")
    if not isinstance(count, int) or count < 2:
        raise ValueError("count must be an integer of at least two")
    direction, spread = _directions_and_spread(directions, spreading)
    if phase is not None and seed is not None:
        raise ValueError("supply either phase or seed")

    if h5_file is not None:
        body = BodyClass(str(h5_file))
        body.bodyNumber = 1
        body.readH5file()
        bem_omega = np.asarray(body.hydroData["simulation_parameters"]["w"]).ravel()
        if len(bem_omega) < 2 or not np.isfinite(bem_omega).all():
            raise ValueError("hydrodynamic frequency range is invalid")
        bem_min, bem_max = float(bem_omega.min()), float(bem_omega.max())
    elif frequency_range is None:
        raise ValueError("a BEM file or explicit frequency_range is required")
    if frequency_range is None:
        omega_min, omega_max = bem_min, bem_max
    else:
        limits = np.asarray(frequency_range, dtype=float)
        if (limits.shape != (2,) or not np.isfinite(limits).all()
                or limits[0] <= 0 or limits[1] <= limits[0]):
            raise ValueError("frequency_range needs two increasing positive rad/s values")
        omega_min, omega_max = map(float, limits)
        if h5_file is not None:
            # Current MATLAB waveClass replaces out-of-BEM endpoints with
            # the corresponding BEM limit before building EqualEnergy bins.
            if omega_min < bem_min or omega_min > bem_max:
                omega_min = bem_min
            if omega_max < bem_min or omega_max > bem_max:
                omega_max = bem_max
            if omega_max <= omega_min:
                raise ValueError("frequency_range has no interval inside the BEM range")
    # Traditional uses a uniform grid; EqualEnergy integrates on 500,000
    # intervals before locating each cumulative-energy boundary.
    dense_omega = np.linspace(
        omega_min, omega_max, count if scheme == "traditional" else 500_001,
    )
    frequency = dense_omega / (2 * np.pi)
    b_pm = 1.25 * (1 / peak_period)**4
    a_pm = b_pm * (significant_height / 2)**2
    spectrum_hz = a_pm * frequency**-5 * np.exp(-b_pm * frequency**-4)
    if gamma is not None:
        peak_frequency = 1 / peak_period
        sigma = np.where(frequency <= peak_frequency, .07, .09)
        peak = gamma**np.exp(
            -(frequency - peak_frequency)**2
            / (2 * sigma**2 * peak_frequency**2)
        )
        spectrum_hz *= (1 - .287 * np.log(gamma)) * peak
    if scheme == "traditional":
        omega = dense_omega
        d_omega = np.full(count, (omega_max - omega_min) / (count - 1))
        spectral_amplitude = spectrum_hz / np.pi
    else:
        df = frequency[1] - frequency[0]
        integrated = np.empty(len(frequency))
        integrated[0] = 0.0
        integrated[1:] = np.cumsum((spectrum_hz[:-1] + spectrum_hz[1:]) * df / 2)
        energy_per_bin = integrated[-1] / (count + 1)
        boundaries = np.zeros(count + 2, dtype=int)
        for k in range(1, count + 2):
            candidate = int(np.searchsorted(integrated, k * energy_per_bin))
            candidate = min(max(candidate, boundaries[k - 1] + 1), len(integrated) - 1)
            earlier = candidate - 1
            if earlier > boundaries[k - 1] and (
                abs(integrated[earlier] - k * energy_per_bin)
                <= abs(integrated[candidate] - k * energy_per_bin)
            ):
                candidate = earlier
            # MATLAB's ``wn(k+1) = wn(k) + wna(k)`` moves one grid point beyond
            # the nearest cumulative-energy sample.
            boundaries[k] = min(candidate + 1, len(integrated) - 1)
        indices = boundaries[1:-1]
        omega = dense_omega[indices]
        d_omega = np.diff(np.r_[dense_omega[0], omega])
        spectral_amplitude = spectrum_hz[indices] / np.pi
    if phase is None:
        phase_array = _random_phases((count, len(direction)), seed, phase_generator)
    else:
        phase_array = np.asarray(phase, dtype=float)
        if phase_array.shape != (count, len(direction)) or not np.isfinite(phase_array).all():
            raise ValueError("phase must have shape (count, number of directions)")
    return IrregularComponents(
        omega=omega, spectral_amplitude=spectral_amplitude,
        d_omega=d_omega, directions=direction, spreading=spread,
        phase=phase_array,
    )


def synthesize_irregular_response(
    h5_file: str | Path,
    components: IrregularComponents,
    *,
    dt: float,
    end_time: float,
    ramp_time: float,
    body_number: int = 1,
    rho: float = 1000.0,
    g: float = 9.81,
    excitation_interpolation: str = "linear",
) -> IrregularResponse:
    """Return wave elevation and all body-DOF excitation at uniform times."""
    if not np.isfinite([dt, end_time, ramp_time, rho, g]).all():
        raise ValueError("time and fluid parameters must be finite")
    if dt <= 0 or end_time < 0 or ramp_time < 0 or rho <= 0 or g <= 0:
        raise ValueError("dt, rho, and g must be positive; times must be nonnegative")
    if excitation_interpolation not in ("linear", "spline"):
        raise ValueError("excitation_interpolation must be linear or spline")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    omega = np.asarray(components.omega, dtype=float).ravel()
    amplitude = np.asarray(components.spectral_amplitude, dtype=float).ravel()
    d_omega = np.asarray(components.d_omega, dtype=float).ravel()
    direction, spread = _directions_and_spread(
        components.directions, components.spreading,
    )
    phase = np.asarray(components.phase, dtype=float)
    if (len(omega) < 2 or amplitude.shape != omega.shape
            or d_omega.shape != omega.shape
            or phase.shape != (len(omega), len(direction))
            or not np.isfinite(omega).all() or not np.isfinite(amplitude).all()
            or not np.isfinite(d_omega).all() or not np.isfinite(phase).all()
            or np.any(np.diff(omega) <= 0) or np.any(amplitude < 0)
            or np.any(d_omega <= 0)):
        raise ValueError("wave components have inconsistent or invalid values")

    body = BodyClass(str(h5_file))
    body.bodyNumber = body_number
    body.readH5file()
    body_dof = int(np.asarray(body.dof).item())
    if excitation_interpolation == "linear":
        body.irrExcitation(omega, len(omega), direction, rho, g)
        real = np.transpose(body.hydroForce["fExt"]["re"], (1, 0, 2))
        imaginary = np.transpose(body.hydroForce["fExt"]["im"], (1, 0, 2))
        mean_drift = np.transpose(body.hydroForce["fExt"]["md"], (1, 0, 2))
    else:
        # MATLAB's passive-yaw preprocessing uses spline interpolation in
        # frequency. At zero yaw, a tabulated BEM direction is selected
        # exactly; varying yaw needs a separate direction-interpolation path.
        bem = body.hydroData["simulation_parameters"]
        bem_frequency = np.asarray(bem["w"], dtype=float).ravel()
        bem_direction = np.asarray(bem["wave_dir"], dtype=float).ravel()
        hydro = body.hydroData["hydro_coeffs"]
        interpolated = []
        for field in (hydro["excitation"]["re"],
                      hydro["excitation"]["im"],
                      hydro["mean_drift"]):
            values = np.asarray(field, dtype=float)
            samples = np.empty((len(omega), len(direction), body_dof))
            for index, heading in enumerate(direction):
                matches = np.flatnonzero(np.isclose(
                    bem_direction, heading, rtol=0, atol=1e-10,
                ))
                if len(matches) != 1:
                    raise ValueError(
                        "spline excitation currently needs tabulated BEM directions"
                    )
                samples[:, index] = CubicSpline(
                    bem_frequency, values[:, matches[0], :], axis=1,
                )(omega).T * rho * g
            interpolated.append(samples)
        real, imaginary, mean_drift = interpolated
    wave_energy = amplitude[:, None] * d_omega[:, None] * spread[None, :]
    height = np.sqrt(wave_energy)
    drift_force = np.einsum("fd,fdc->c", wave_energy, mean_drift)

    time = np.arange(steps + 1) * dt
    elevation = np.zeros(len(time))
    excitation = np.zeros((len(time), body_dof))
    for start in range(0, len(time), 128):
        stop = min(start + 128, len(time))
        t = time[start:stop]
        phase_angle = t[:, None, None] * omega[None, :, None] + phase[None, :, :]
        cosine = np.cos(phase_angle) * height[None, :, :]
        sine = np.sin(phase_angle) * height[None, :, :]
        elevation[start:stop] = cosine.sum(axis=(1, 2))
        excitation[start:stop] = (
            drift_force + np.einsum("tfd,fdc->tc", cosine, real)
            - np.einsum("tfd,fdc->tc", sine, imaginary)
        )
    ramp = np.ones(len(time))
    if ramp_time > 0:
        early = time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    return IrregularResponse(
        time=time, elevation=ramp * elevation,
        excitation_force=ramp[:, None] * excitation,
    )


def synthesize_full_directional_response(
    h5_file: str | Path,
    components: FullDirectionalComponents,
    *,
    dt: float,
    end_time: float,
    ramp_time: float,
    body_number: int = 1,
    rho: float = 1000.0,
    g: float = 9.81,
    excitation_interpolation: str = "linear",
    force_quadrature: str = "integrated",
) -> IrregularResponse:
    """Synthesize elevation and excitation from a full directional spectrum.

    Both wave and force quadrature use frequency and heading bin widths by
    default. ``force_quadrature="matlab_omitted"`` reproduces the pinned
    MATLAB full-directional force block's omission of heading width; it is a
    source diagnostic, not the physical default. The optional
    ``spline_frequency`` interpolation also reproduces its BEM preprocessing.
    """
    if not isinstance(components, FullDirectionalComponents):
        raise TypeError("components must be FullDirectionalComponents")
    if excitation_interpolation not in ("linear", "spline_frequency"):
        raise ValueError("excitation_interpolation must be linear or spline_frequency")
    if force_quadrature not in ("integrated", "matlab_omitted"):
        raise ValueError("force_quadrature must be integrated or matlab_omitted")
    if (not np.isfinite([dt, end_time, ramp_time, rho, g]).all()
            or dt <= 0 or end_time < 0 or ramp_time < 0 or rho <= 0 or g <= 0):
        raise ValueError("time and fluid parameters must be finite and valid")
    steps = round(end_time / dt)
    if not np.isclose(steps * dt, end_time, rtol=0, atol=1e-10):
        raise ValueError("end_time must be an integer multiple of dt")
    omega = np.asarray(components.omega, dtype=float).ravel()
    heading = np.asarray(components.directions, dtype=float).ravel()
    spectrum = np.asarray(components.spectral_amplitude, dtype=float)
    d_omega = np.asarray(components.d_omega, dtype=float).ravel()
    d_theta = np.asarray(components.d_theta, dtype=float).ravel()
    phase = np.asarray(components.phase, dtype=float)
    shape = (len(omega), len(heading))
    if (len(omega) < 2 or len(heading) < 2
            or spectrum.shape != shape or phase.shape != shape
            or d_omega.shape != (len(omega),)
            or d_theta.shape != (len(heading),)
            or not all(np.isfinite(values).all()
                       for values in (omega, heading, spectrum, d_omega,
                                      d_theta, phase))
            or np.any(spectrum < 0) or np.any(d_omega <= 0)
            or np.any(d_theta <= 0) or np.any(np.diff(omega) <= 0)
            or np.any(np.diff(heading) <= 0)):
        raise ValueError("full-directional components have invalid dimensions or values")

    body = BodyClass(str(h5_file))
    body.bodyNumber = body_number
    body.readH5file()
    hydro = body.hydroData["hydro_coeffs"]
    bem = body.hydroData["simulation_parameters"]
    bem_omega = np.asarray(bem["w"], dtype=float).ravel()
    bem_heading = np.mod(np.asarray(bem["wave_dir"], dtype=float).ravel(), 360)
    order = np.argsort(bem_heading)
    bem_heading = bem_heading[order]
    if (len(bem_heading) < 2 or np.any(np.diff(bem_heading) <= 0)
            or np.any(np.diff(bem_omega) <= 0)):
        raise ValueError("BEM frequency and heading grids must increase")
    query_heading = ((heading - bem_heading[0]) % 360) + bem_heading[0]
    query_frequency, query_direction = np.meshgrid(omega, query_heading,
                                                   indexing="ij")
    query = np.column_stack((query_direction.ravel(),
                             query_frequency.ravel()))

    def interpolate(field):
        values = np.asarray(field, dtype=float)
        if excitation_interpolation == "spline_frequency":
            values = CubicSpline(bem_omega, values, axis=2)(omega)
            frequency_grid = omega
        else:
            frequency_grid = bem_omega
        values = np.transpose(values, (1, 2, 0))[order]
        extended = np.concatenate((values, values[:1]), axis=0)
        grid = np.r_[bem_heading, bem_heading[0] + 360]
        return RegularGridInterpolator(
            (grid, frequency_grid), extended, bounds_error=True,
        )(query).reshape(*shape, 6) * rho * g

    real = interpolate(hydro["excitation"]["re"])
    imaginary = interpolate(hydro["excitation"]["im"])
    mean_drift = interpolate(hydro["mean_drift"])
    energy = spectrum * d_omega[:, None] * d_theta[None, :]
    height = np.sqrt(energy).ravel()
    force_energy = (energy if force_quadrature == "integrated"
                    else spectrum * d_omega[:, None])
    force_height = np.sqrt(force_energy).ravel()
    drift_force = np.einsum("fd,fdc->c", force_energy, mean_drift)
    omega_flat = np.repeat(omega, len(heading))
    phase_flat = phase.ravel()
    real_flat = real.reshape(-1, 6)
    imaginary_flat = imaginary.reshape(-1, 6)

    time = np.arange(steps + 1) * dt
    elevation = np.zeros(len(time))
    excitation = np.zeros((len(time), 6))
    for start in range(0, len(time), 128):
        stop = min(start + 128, len(time))
        angle = time[start:stop, None] * omega_flat[None] + phase_flat[None]
        cosine = np.cos(angle)
        sine = np.sin(angle)
        elevation[start:stop] = (cosine * height[None]).sum(axis=1)
        excitation[start:stop] = (
            drift_force + (cosine * force_height[None]) @ real_flat
            - (sine * force_height[None]) @ imaginary_flat
        )
    ramp = np.ones(len(time))
    if ramp_time > 0:
        early = time < ramp_time
        ramp[early] = (1 - np.cos(np.pi * time[early] / ramp_time)) / 2
    return IrregularResponse(time, ramp * elevation,
                             ramp[:, None] * excitation)


def synthesize_multiple_irregular_response(
    h5_file: str | Path,
    seas: Sequence[IrregularComponents],
    *,
    dt: float,
    end_time: float,
    ramp_time: float,
    body_number: int = 1,
    rho: float = 1000.0,
    g: float = 9.81,
    excitation_interpolation: str = "linear",
) -> IrregularResponse:
    """Sum independent realized sea elevations and linear excitation forces.

    Each sea keeps its own spectrum, incident directions, and phase matrix.
    WEC-Sim's multiple-wave input combines their forcing before the body
    dynamics, so radiation and inertia are applied once to the total motion.
    """
    seas = tuple(seas)
    if not seas or any(not isinstance(sea, IrregularComponents) for sea in seas):
        raise ValueError("seas must contain at least one irregular-wave realization")
    time = elevation = excitation_force = None
    for sea in seas:
        response = synthesize_irregular_response(
            h5_file, sea, dt=dt, end_time=end_time, ramp_time=ramp_time,
            body_number=body_number, rho=rho, g=g,
            excitation_interpolation=excitation_interpolation,
        )
        if time is None:
            time = response.time
            elevation = response.elevation.copy()
            excitation_force = response.excitation_force.copy()
        else:
            elevation += response.elevation
            excitation_force += response.excitation_force
    return IrregularResponse(time, elevation, excitation_force)


def _directions_and_spread(directions, spreading):
    direction = np.asarray(directions, dtype=float).ravel()
    spread = np.asarray(spreading, dtype=float).ravel()
    if (len(direction) == 0 or direction.shape != spread.shape
            or not np.isfinite(direction).all() or not np.isfinite(spread).all()
            or np.any(spread < 0) or not np.isclose(spread.sum(), 1.0, atol=1e-10)):
        raise ValueError("directions and spreading must be finite, same-sized, and sum to one")
    return direction, spread
