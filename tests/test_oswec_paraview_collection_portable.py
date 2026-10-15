"""Pair the published OSWEC ParaView collection in ordinary Python CI."""

from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

import numpy as np

from wecsim.visualization.paraviewClass import ParaviewClass


REFERENCE = Path(__file__).parent / "fixtures/matlab/oswec_nonlinear_viz_collection.zip"


def test_published_oswec_paraview_collection(tmp_path):
    """Check every published time and path, plus the body display manifest."""
    with ZipFile(REFERENCE) as source:
        source_root = ElementTree.fromstring(source.read("OSWEC.pvd"))
        source_bodies = source.read("bodies.txt")

    times = np.arange(1201) * 0.1
    actual = ParaviewClass(None).write_paraview_response(
        times, tmp_path, model_name="OSWEC",
        bodies=[("flap", (1, 1, 0), 1), ("base", (1, 1, 0), 1)],
    )
    assert (tmp_path / "bodies.txt").read_bytes() == source_bodies

    actual_root = ElementTree.parse(actual).getroot()
    assert actual_root.tag == source_root.tag == "VTKFile"
    assert actual_root.attrib == source_root.attrib == {
        "type": "Collection", "version": "0.1",
    }
    source_entries = [item.attrib for item in source_root.findall("./Collection/DataSet")]
    actual_entries = [item.attrib for item in actual_root.findall("./Collection/DataSet")]
    assert len(source_entries) == len(actual_entries) == 3 * len(times) == 3603
    assert actual_entries == source_entries
