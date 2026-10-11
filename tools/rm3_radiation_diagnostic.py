"""Compare source BEM and fitted radiation damping for RM3 common surge.

The fitted response deliberately omits HDF5 feedthrough D because the pinned
MATLAB WEC-Sim body class uses zero D for the state-space radiation block.
This diagnoses the model used by Cases 5–6; it does not certify passivity.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np


def common_surge_curves(hydro_file):
    """Return frequencies and effective common-surge damping in N s/m.

    Row zero omits cross-body radiation, as in Case 5. Row one includes it,
    as in Case 6. The value is the total surge force divided by a shared
    surge velocity, so both body outputs contribute.
    """
    with h5py.File(hydro_file) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        if not np.isfinite(rho) or rho <= 0:
            raise ValueError("hydrodynamic density must be positive")
        source_frequency = np.asarray(h5["simulation_parameters/w"]).ravel()
        if (len(source_frequency) < 2
                or not np.all(np.diff(source_frequency) > 0)):
            raise ValueError("hydrodynamic frequencies must increase")
        fitted_frequency = np.r_[0.0, source_frequency]
        source = np.zeros((2, len(source_frequency)))
        fitted = np.zeros((2, len(fitted_frequency)))

        for output_body in range(2):
            prefix = f"body{output_body + 1}/hydro_coeffs/radiation_damping"
            source_blocks = np.asarray(h5[f"{prefix}/all"])
            fit = h5[f"{prefix}/state_space"]
            orders = np.asarray(fit["it"], dtype=int)
            A = np.asarray(fit["A/all"])
            B = np.asarray(fit["B/all"])
            C = np.asarray(fit["C/all"])
            if source_blocks.shape != (6, 12, len(source_frequency)):
                raise ValueError("expected two six-DOF RM3 hydrodynamic bodies")

            for input_body in range(2):
                column = 6 * input_body
                variants = ((0, 1) if input_body == output_body else (1,))
                for coupled in variants:
                    source[coupled] += rho * source_blocks[0, column]
                    order = orders[0, column]
                    if order < 0 or order > A.shape[-1]:
                        raise ValueError("invalid fitted radiation order")
                    if order == 0:
                        continue
                    state_A = A[0, column, :order, :order]
                    state_B = B[0, column, :order, 0]
                    state_C = rho * C[0, column, 0, :order]
                    identity = np.eye(order)
                    for index, frequency in enumerate(fitted_frequency):
                        fitted[coupled, index] += np.real(
                            state_C @ np.linalg.solve(
                                1j * frequency * identity - state_A,
                                state_B,
                            )
                        )

    if not (np.isfinite(source).all() and np.isfinite(fitted).all()):
        raise ValueError("radiation damping contains nonfinite values")
    return source_frequency, source, fitted_frequency, fitted


def common_surge_feedthrough(hydro_file):
    """Sum the saved HDF5 direct terms for the two common-surge variants.

    The pinned MATLAB body class instead assigns a zero direct term when it
    builds the state-space radiation block. These sums are diagnostic only.
    """
    with h5py.File(hydro_file) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        feedthrough = np.zeros(2)
        for output_body in range(2):
            direct = np.asarray(h5[
                f"body{output_body + 1}/hydro_coeffs/radiation_damping/"
                "state_space/D/all"
            ])
            if direct.shape != (6, 12):
                raise ValueError("expected two six-DOF RM3 direct-term blocks")
            for input_body in range(2):
                term = rho * direct[0, 6 * input_body]
                if input_body == output_body:
                    feedthrough[0] += term
                feedthrough[1] += term
    return feedthrough


def active_mode_minimum_damping(hydro_file, *, return_modes=False):
    """Minimum radiation damping on RM3's four moving joint coordinates.

    The coordinates are common surge, float heave, spar heave, and shared
    pitch. Pitch is expressed as equivalent travel at a 20 m lever, so each
    projected coordinate has units of metres. A negative eigenvalue of the
    Hermitian transfer matrix permits a harmonic motion that receives energy
    from radiation. The published source fit omits the saved HDF5 D terms.
    With ``return_modes=True``, append the fitted least-damped eigenvector
    in the same four energy-conjugate coordinates at each frequency.
    """
    with h5py.File(hydro_file) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        frequencies = np.asarray(h5["simulation_parameters/w"]).ravel()
        if rho <= 0 or not np.isfinite(rho) or not np.all(np.diff(frequencies) > 0):
            raise ValueError("invalid hydrodynamic density or frequencies")
        fitted_frequencies = np.r_[0.0, frequencies]
        source = np.zeros((len(frequencies), 12, 12))
        fitted = np.zeros((len(fitted_frequencies), 12, 12), dtype=complex)
        motion = np.zeros((12, 4))
        for body in range(2):
            center_z = float(np.asarray(h5[f"body{body + 1}/properties/cg"]).ravel()[2])
            row = 6 * body
            motion[row, 0] = 1
            motion[row, 3] = center_z / 20
            motion[row + 2, body + 1] = 1
            motion[row + 4, 3] = 1 / 20
            radiation = h5[f"body{body + 1}/hydro_coeffs/radiation_damping"]
            source[:, row:row + 6] = rho * np.moveaxis(
                np.asarray(radiation["all"]), -1, 0,
            )
            fit = radiation["state_space"]
            orders = np.asarray(fit["it"], dtype=int)
            A = np.asarray(fit["A/all"])
            B = np.asarray(fit["B/all"])
            C = np.asarray(fit["C/all"])
            for output in range(6):
                for input_dof in range(12):
                    order = orders[output, input_dof]
                    if order < 0 or order > A.shape[-1]:
                        raise ValueError("invalid fitted radiation order")
                    if order == 0:
                        continue
                    state_A = A[output, input_dof, :order, :order]
                    state_B = B[output, input_dof, :order, 0]
                    state_C = rho * C[output, input_dof, 0, :order]
                    identity = np.eye(order)
                    fitted[:, row + output, input_dof] = [
                        state_C @ np.linalg.solve(
                            1j * frequency * identity - state_A, state_B,
                        )
                        for frequency in fitted_frequencies
                    ]

    def minimum_eigenvalue(matrix):
        projected = motion.T @ matrix @ motion
        hermitian = (projected + projected.conj().swapaxes(-1, -2)) / 2
        return np.linalg.eigvalsh(hermitian)[:, 0]

    def minimum_eigenpair(matrix):
        projected = motion.T @ matrix @ motion
        hermitian = (projected + projected.conj().swapaxes(-1, -2)) / 2
        values, vectors = np.linalg.eigh(hermitian)
        return values[:, 0], vectors[:, :, 0]

    source_minimum = np.empty((2, len(frequencies)))
    fitted_minimum = np.empty((2, len(fitted_frequencies)))
    if return_modes:
        fitted_modes = np.empty((2, len(fitted_frequencies), 4), dtype=complex)
    for coupled in (False, True):
        source_variant = source.copy()
        fitted_variant = fitted.copy()
        if not coupled:
            source_variant[:, :6, 6:] = 0
            source_variant[:, 6:, :6] = 0
            fitted_variant[:, :6, 6:] = 0
            fitted_variant[:, 6:, :6] = 0
        source_minimum[int(coupled)] = minimum_eigenvalue(source_variant)
        if return_modes:
            fitted_minimum[int(coupled)], fitted_modes[int(coupled)] = (
                minimum_eigenpair(fitted_variant))
        else:
            fitted_minimum[int(coupled)] = minimum_eigenvalue(fitted_variant)
    result = (frequencies, source_minimum, fitted_frequencies, fitted_minimum)
    return (*result, fitted_modes) if return_modes else result


def plot_curves(source_frequency, source, fitted_frequency, fitted, output,
                feedthrough=None):
    """Save a two-panel diagnostic plot for the published RM3 Cases 5–6."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        1, 2, figsize=(11, 4.5), sharey=True, layout="constrained",
    )
    for coupled, axis in enumerate(axes):
        axis.axhline(0, color="#303030", linewidth=0.8)
        axis.axvline(2 * np.pi / 8, color="#777777", linestyle=":",
                     linewidth=1.5)
        axis.plot(source_frequency, source[coupled] / 1000,
                  color="#236b82", linewidth=2, label="Source BEM damping")
        axis.plot(fitted_frequency, fitted[coupled] / 1000,
                  color="#c74558", linewidth=2, label="MATLAB state-space fit")
        if feedthrough is not None:
            axis.plot(fitted_frequency,
                      (fitted[coupled] + feedthrough[coupled]) / 1000,
                      color="#c74558", linestyle="--", linewidth=1.4,
                      label="Fit with saved HDF5 direct term")
        axis.fill_between(fitted_frequency, fitted[coupled] / 1000, 0,
                          where=fitted[coupled] < 0,
                          color="#c74558", alpha=0.12)
        axis.set_xlim(0, 1.5)
        axis.set_ylim(-30, 220)
        axis.set_xlabel("Angular frequency (rad/s)")
        axis.set_title("Cross-body radiation " + ("on" if coupled else "off"))
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Effective common-surge damping (kN s/m)")
    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle("RM3 Cases 5–6: fitted damping changes sign near zero frequency")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_active_modes(source_frequency, source, fitted_frequency, fitted, output):
    """Plot the least-damped physical joint motion at each frequency."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True,
                             layout="constrained")
    for coupled, axis in enumerate(axes):
        axis.axhline(0, color="#303030", linewidth=0.8)
        axis.axvline(2 * np.pi / 8, color="#777777", linestyle=":",
                     linewidth=1.5)
        axis.plot(source_frequency, source[coupled] / 1000,
                  color="#236b82", linewidth=2, label="Source BEM")
        axis.plot(fitted_frequency, fitted[coupled] / 1000,
                  color="#c74558", linewidth=2, label="MATLAB fit")
        axis.fill_between(fitted_frequency, fitted[coupled] / 1000, 0,
                          where=fitted[coupled] < 0,
                          color="#c74558", alpha=0.12)
        axis.set_xlim(0, 1.5)
        axis.set_ylim(-25, 80)
        axis.set_xlabel("Angular frequency (rad/s)")
        axis.set_title("Cross-body radiation " + ("on" if coupled else "off"))
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Least active-mode damping (kN s/m)")
    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle("RM3: radiation damping in the four moving joint coordinates")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hydro_file", type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--active-plot", type=Path)
    args = parser.parse_args()
    source_frequency, source, fitted_frequency, fitted = common_surge_curves(
        args.hydro_file,
    )
    feedthrough = common_surge_feedthrough(args.hydro_file)
    for coupled in (0, 1):
        print(f"cross-body {'on' if coupled else 'off'}: "
              f"fit DC {fitted[coupled, 0] / 1000:.4f} kN s/m; "
              f"fit with saved D "
              f"{(fitted[coupled, 0] + feedthrough[coupled]) / 1000:.4f} "
              "kN s/m; "
              f"source minimum {source[coupled].min() / 1000:.4g} kN s/m")
    if args.plot:
        plot_curves(source_frequency, source, fitted_frequency, fitted,
                    args.plot, feedthrough)
        print(args.plot)
    if args.active_plot:
        active = active_mode_minimum_damping(args.hydro_file)
        for coupled in (0, 1):
            print(f"active modes, cross-body {'on' if coupled else 'off'}: "
                  f"source minimum {active[1][coupled].min() / 1000:.4f} "
                  "kN s/m; "
                  f"fit DC minimum {active[3][coupled, 0] / 1000:.4f} kN s/m")
        plot_active_modes(*active, args.active_plot)
        print(args.active_plot)


if __name__ == "__main__":
    main()
