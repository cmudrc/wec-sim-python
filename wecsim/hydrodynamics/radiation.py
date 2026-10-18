"""Source-convention radiation diagnostics for pinned WEC-Sim HDF5 data."""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import lsim


@dataclass(frozen=True)
class RadiationStateSpace:
    """Continuous fitted radiation: x' = A x + B v, F = C x.

    The fitted force is a resisting force; the dynamics subtracts it. The
    pinned MATLAB body block discards the HDF5 direct term D. This model is
    intended for explicit source comparisons, not as a passive replacement
    for convolution radiation.
    """

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray

    def __post_init__(self):
        for name in ("A", "B", "C"):
            object.__setattr__(self, name,
                               np.asarray(getattr(self, name), dtype=float))
        if self.A.ndim != 2:
            raise ValueError("radiation state-space matrices need finite two-body shapes")
        count = self.A.shape[0]
        if (count < 1 or self.A.shape != (count, count)
                or np.shape(self.B) != (count, 12)
                or np.shape(self.C) != (12, count)
                or not all(np.isfinite(matrix).all()
                           for matrix in (self.A, self.B, self.C))):
            raise ValueError("radiation state-space matrices need finite two-body shapes")

    @classmethod
    def from_rm3_hdf5(cls, hydro_file: str | Path, *,
                      body_to_body: bool = False):
        """Load only channels active in RM3's surge/heave/pitch joint.

        Inactive body velocities and forces are identically zero under this
        joint. Removing their independent fitted states does not change the
        four-coordinate motion. The source sets D=0 even when HDF5 stores it.
        """
        if not isinstance(body_to_body, bool):
            raise TypeError("body_to_body must be a boolean")
        channels = []
        with h5py.File(hydro_file) as h5:
            rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
            if not np.isfinite(rho) or rho <= 0:
                raise ValueError("hydrodynamic density must be positive")
            for output_body in range(2):
                fit = h5[
                    f"body{output_body + 1}/hydro_coeffs/"
                    "radiation_damping/state_space"
                ]
                raw_orders = np.asarray(fit["it"])
                all_a = np.asarray(fit["A/all"])
                all_b = np.asarray(fit["B/all"])
                all_c = np.asarray(fit["C/all"])
                if (raw_orders.shape != (6, 12)
                        or not np.isfinite(raw_orders).all()
                        or not np.all(raw_orders == np.floor(raw_orders))
                        or all_a.ndim != 4 or all_a.shape[:2] != (6, 12)
                        or all_a.shape[2] != all_a.shape[3]
                        or all_b.shape != (6, 12, all_a.shape[2], 1)
                        or all_c.shape != (6, 12, 1, all_a.shape[2])
                        or not all(np.isfinite(matrix).all()
                                   for matrix in (all_a, all_b, all_c))):
                    raise ValueError("expected a two-body six-DOF RM3 radiation fit")
                orders = raw_orders.astype(int)
                for output_dof in (0, 2, 4):
                    for input_body in range(2):
                        if output_body != input_body and not body_to_body:
                            continue
                        for input_dof in (0, 2, 4):
                            column = 6 * input_body + input_dof
                            order = int(orders[output_dof, column])
                            if order < 0 or order > all_a.shape[-1]:
                                raise ValueError("invalid fitted radiation order")
                            if order:
                                channels.append((
                                    6 * output_body + output_dof, column,
                                    all_a[output_dof, column, :order, :order],
                                    all_b[output_dof, column, :order, 0],
                                    rho * all_c[output_dof, column, 0, :order],
                                ))
        count = sum(len(channel[3]) for channel in channels)
        a = np.zeros((count, count))
        b = np.zeros((count, 12))
        c = np.zeros((12, count))
        start = 0
        for output, input_dof, block_a, block_b, block_c in channels:
            stop = start + len(block_b)
            a[start:stop, start:stop] = block_a
            b[start:stop, input_dof] = block_b
            c[output, start:stop] = block_c
            start = stop
        return cls(a, b, c)


def replay_rm3_fitted_radiation(
    hydro_file: str | Path, time: np.ndarray, body_velocity: np.ndarray, *,
    body_to_body: bool = False,
) -> np.ndarray:
    """Replay RM3's MATLAB state-space radiation on prescribed body motion.

    Return force with shape ``(time, body, six DOFs)``. The pinned MATLAB
    block multiplies HDF5 C by water density and sets D to zero, even though
    HDF5 contains a direct term. Each fitted channel starts with zero state.
    Linear interpolation between saved velocities approximates the source's
    internal solver input; this function does not advance a WEC trajectory.
    The published fit has negative damping in active joint coordinates, so
    this replay is separate from the default convolution dynamics.
    """
    t = np.asarray(time, dtype=float)
    velocity = np.asarray(body_velocity, dtype=float)
    if (t.ndim != 1 or len(t) < 2 or t[0] != 0
            or not np.isfinite(t).all() or not np.all(np.diff(t) > 0)
            or not np.allclose(np.diff(t), t[1] - t[0], rtol=1e-8, atol=1e-12)
            or velocity.shape != (len(t), 2, 6)
            or not np.isfinite(velocity).all()):
        raise ValueError("radiation replay needs uniform time from zero and two finite six-DOF velocity histories")
    if not isinstance(body_to_body, bool):
        raise TypeError("body_to_body must be a boolean")

    force = np.zeros_like(velocity)
    with h5py.File(hydro_file) as h5:
        rho = float(np.asarray(h5["simulation_parameters/rho"]).item())
        if not np.isfinite(rho) or rho <= 0:
            raise ValueError("hydrodynamic density must be positive")
        for output_body in range(2):
            fit = h5[f"body{output_body + 1}/hydro_coeffs/radiation_damping/state_space"]
            orders = np.asarray(fit["it"])
            A = np.asarray(fit["A/all"])
            B = np.asarray(fit["B/all"])
            C = np.asarray(fit["C/all"])
            if (orders.shape != (6, 12)
                    or not np.isfinite(orders).all()
                    or not np.all(orders == np.floor(orders))
                    or A.ndim != 4 or A.shape[:2] != (6, 12)
                    or A.shape[2] != A.shape[3]
                    or B.shape != (6, 12, A.shape[2], 1)
                    or C.shape != (6, 12, 1, A.shape[2])
                    or not all(np.isfinite(matrix).all()
                               for matrix in (A, B, C))):
                raise ValueError("expected a two-body six-DOF radiation fit")
            for output_dof in range(6):
                for input_body in range(2):
                    if input_body != output_body and not body_to_body:
                        continue
                    for input_dof in range(6):
                        column = 6 * input_body + input_dof
                        order = int(orders[output_dof, column])
                        if order < 0 or order > A.shape[-1]:
                            raise ValueError("invalid fitted radiation order")
                        if order == 0:
                            continue
                        system = (
                            A[output_dof, column, :order, :order],
                            B[output_dof, column, :order, :1],
                            C[output_dof, column, :1, :order],
                            np.zeros((1, 1)),
                        )
                        _, response, _ = lsim(
                            system, U=velocity[:, input_body, input_dof], T=t,
                        )
                        force[:, output_body, output_dof] += rho * response
    return force
