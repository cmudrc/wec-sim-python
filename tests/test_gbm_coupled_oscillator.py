"""Analytic rigid-heave/flexible-mode check without MATLAB or a BEM file."""

import numpy as np

import wecsim.dynamics.gbmFloating as gbm


class _TwoCoordinateHydro:
    """Supply a positive coupled mass and no damping or restoring."""

    def __init__(self, _path):
        self.dof = np.array([7])
        self.dof_gbm = np.array([1])
        self.cg = np.zeros(3)
        self.dispVol = np.array([.1])

    def readH5file(self):
        pass

    def hydroForcePre(self, *_args):
        added = np.zeros((7, 7))
        added[2, 6] = added[6, 2] = 2
        excitation = np.zeros(7)
        excitation[6] = 2
        self.hydroForce = {
            "fAddedMass": added,
            "fDamping": np.zeros((7, 7)),
            "linearHydroRestCoef": np.zeros((7, 7)),
            "fExt": {"re": excitation, "im": np.zeros(7)},
            "gbm": {
                "mass_ff": np.array([[4.]]),
                "damping": np.zeros((1, 1)),
                "stiffness": np.zeros((1, 1)),
            },
        }


def test_regular_gbm_solver_respects_coupled_mass_and_wave_work(monkeypatch):
    # Isolate the integrator from HDF5 preprocessing. With 0.1 m³ displaced
    # volume, equilibrium mass is 100 kg at the solver's 1,000 kg/m³ density.
    monkeypatch.setattr(gbm, "BodyClass", _TwoCoordinateHydro)
    response = gbm.solve_floating_gbm_regular(
        "synthetic.h5", dt=.05, end_time=2, height=2, period=8,
        pitch_inertia=20, mass="equilibrium", ramp_time=0,
    )
    frequency = 2 * np.pi / 8
    mass = np.array([[100., 2.], [2., 4.]])
    force = np.array([0., 2.])
    acceleration_amplitude = np.linalg.solve(mass, force)
    time = response.time
    expected_position = ((1 - np.cos(frequency * time))[:, None]
                         * acceleration_amplitude / frequency**2)
    expected_speed = (np.sin(frequency * time)[:, None]
                      * acceleration_amplitude / frequency)

    position = np.column_stack((response.body_position[:, 2],
                                response.mode_position[:, 0]))
    speed = np.column_stack((response.body_velocity[:, 2],
                             response.mode_velocity[:, 0]))
    np.testing.assert_allclose(position, expected_position, rtol=0, atol=1e-8)
    np.testing.assert_allclose(speed, expected_speed, rtol=0, atol=1e-8)
    np.testing.assert_allclose(
        response.mode_acceleration[:, 0],
        acceleration_amplitude[1] * np.cos(frequency * time),
        rtol=0, atol=1e-8,
    )
    np.testing.assert_allclose(speed @ mass.T,
                               np.sin(frequency * time)[:, None]
                               * force / frequency, rtol=0, atol=1e-7)
    wave_work = np.trapezoid(
        2 * np.cos(frequency * time) * speed[:, 1], time,
    )
    kinetic_energy = speed[-1] @ mass @ speed[-1] / 2
    assert abs(wave_work - kinetic_energy) < 1e-3
    np.testing.assert_allclose(response.wave_elevation,
                               np.cos(frequency * time), rtol=0, atol=1e-12)
    np.testing.assert_array_equal(response.body_position[:, [0, 4]], 0)
