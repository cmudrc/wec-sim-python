"""Portable published-source gate for public RM3 linear-generator motion."""

from pathlib import Path

import numpy as np

from wecsim import DirectLinearGenerator, RegularWave, WEC


FIXTURE = Path(__file__).parent / "fixtures/matlab/rm3_direct_generator"


def _gate(actual, source, limit, label):
    assert actual.shape == source.shape, label
    error = float(np.max(np.abs(actual - source)))
    assert np.isfinite(error) and error < limit, (
        f"{label}: {error:.6g} exceeds {limit}"
    )


def test_public_rm3_generator_matches_published_first_20_seconds():
    with np.load(FIXTURE / "published_20s.npz") as source:
        wec = WEC("RM3 direct linear generator")
        floating = wec.body("float", FIXTURE / "rm3.h5")
        spar = wec.body("spar", FIXTURE / "rm3.h5")
        wec.coordinate("float_heave", floating.move("heave"))
        wec.coordinate("spar_heave", spar.move("heave"))
        wec.pto(
            "PTO1", spar.at(0, 0, 0), floating.at(0, 0, 0),
            axis=(0, 0, 1),
            linear_generator=DirectLinearGenerator(
                stator_resistance=4.58, friction=-100, pole_pitch=.072,
                magnet_flux=8, inductance=.285, load_resistance=-117.6471,
            ),
        )
        # The source saves at 0.01 s, but advances the electrical and body
        # states at 0.0005 s. Keep that integration step in the public path.
        result = wec.run(RegularWave(2.5, 8), dt=.0005,
                         end_time=20, ramp_time=100)
        sample = slice(None, None, 20)
        pto = result.ptos["PTO1"]
        generator = pto.linear_generator
        assert generator is not None
        assert np.ptp(source["pto_stroke"]) > .15
        assert np.ptp(source["force"]) > 100

        _gate(result.time[sample], source["time"], 1e-9, "time")
        for index, name in enumerate(("float", "spar")):
            body = result.bodies[name]
            _gate(body.position[sample, 2], source["body_heave"][:, index],
                  1e-5, f"{name} heave")
            _gate(body.velocity[sample, 2], source["body_speed"][:, index],
                  1e-5, f"{name} speed")
        _gate(pto.stroke[sample], source["pto_stroke"], 1e-5, "PTO stroke")
        _gate(pto.velocity[sample], source["pto_speed"], 1e-5, "PTO speed")
        _gate(pto.force[sample], source["force"], .005, "generator force")
        _gate(pto.absorbed_power[sample], source["absorbed_power"],
              .005, "absorbed power")
        _gate(generator.friction_force[sample], source["friction_force"],
              .001, "friction force")
        _gate(generator.electrical_power[sample], source["electrical_power"],
              .005, "electrical power")
        _gate(generator.phase_current[sample], source["phase_current"],
              .0001, "three-phase current")
        _gate(generator.phase_voltage[sample], source["phase_voltage"],
              .01, "three-phase voltage")
