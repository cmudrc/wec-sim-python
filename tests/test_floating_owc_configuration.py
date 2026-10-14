"""Check the public configuration of the published floating OWC layout."""

import pytest

from wecsim import MoorDyn, PMWave, RegularWave, WEC


def _layout(tmp_path):
    # Case construction stores these paths; it does not open the native model.
    hydro = tmp_path / "floating_owc.h5"
    library = tmp_path / "libmoordyn.so"
    lines = tmp_path / "lines.txt"
    for path in (hydro, library, lines):
        path.touch()
    mooring = MoorDyn(library, lines)
    wec = WEC("Floating OWC")
    floater = wec.body(
        "floater", hydro, hydro_body=1,
        inertia=(1.531e9, 1.531e9, 0.1118e9),
    )
    column = wec.body("column", hydro, hydro_body=2, mass=4_493_450)
    return wec, floater, column, mooring


def test_public_floating_owc_preserves_source_attachment_and_pto(tmp_path):
    wec, floater, column, mooring = _layout(tmp_path)
    wec.floating_owc(
        floater, column, moordyn=mooring,
        moordyn_point=floater.at(0, 0, 31.945),
        column_height=50.69, column_diameter=5.89,
        pto_name="turbine_slider", pto_stiffness=15_000,
        pto_damping=60_000,
    )
    case = wec.to_case(RegularWave(4.5, 11.2), dt=.01, end_time=.1)
    assert case["simulation"]["ramp_time"] == 50
    assert [body["hydro_body"] for body in case["bodies"]] == [1, 2]
    assert [body["name"] for body in case["bodies"]] == ["floater", "column"]
    assert case["constraint"]["kind"] == "floating_owc"
    assert case["constraint"]["column_height"] == 50.69
    assert case["constraint"]["column_diameter"] == 5.89
    assert case["pto"] == {
        "kind": "column_slider", "name": "turbine_slider",
        "stiffness": 15_000, "damping": 60_000,
    }
    assert case["mooring"]["session"] is mooring
    assert case["mooring"]["point"] == [0, 0, 31.945]


def test_public_floating_owc_rejects_column_local_mooring_point(tmp_path):
    wec, floater, column, mooring = _layout(tmp_path)
    with pytest.raises(ValueError, match="floater-local"):
        wec.floating_owc(
            floater, column, moordyn=mooring,
            moordyn_point=column.at(0, 0, 31.945),
            column_height=50.69, column_diameter=5.89,
        )


def test_public_floating_owc_rejects_incompatible_sea_and_coordinates(tmp_path):
    wec, floater, column, mooring = _layout(tmp_path)
    wec.floating_owc(
        floater, column, moordyn=mooring,
        column_height=50.69, column_diameter=5.89,
    )
    with pytest.raises(ValueError, match="zero-heading regular waves"):
        wec.to_case(PMWave(4.5, 11.2, seed=1), dt=.01, end_time=.1)
    with pytest.raises(ValueError, match="no radiation memory"):
        wec.to_case(RegularWave(4.5, 11.2), dt=.01, end_time=.1,
                    radiation_memory=.1)
    wec.coordinate("heave", floater.move("heave"))
    with pytest.raises(ValueError, match="zero-heading regular waves"):
        wec.to_case(RegularWave(4.5, 11.2), dt=.01, end_time=.1)
