import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "ustruct")


@pytest.mark.parametrize("ele", ["P1P1_VMS"])
def test_block_compression_ustruct(ele, n_proc):
    folder = os.path.join(base_folder, "block_compression_ustruct", ele)
    field = ["Displacement", "Pressure", "Stress", "Divergence"]
    run_with_reference(folder, field, n_proc)


def test_tensile_adventitia_HGO(n_proc):
    folder = os.path.join(base_folder, "tensile_adventitia_HGO")
    field = ["Displacement", "Velocity", "Stress", "VonMises_stress"]
    run_with_reference(folder, field, n_proc)


def test_LV_Guccione_active(n_proc):
    folder = os.path.join(base_folder, "LV_Guccione_active")
    field = [
        "Displacement",
        "Velocity",
        "Pressure",
        "VonMises_stress",
        "Cauchy_stress",
        "Strain",
        "Jacobian",
    ]
    run_with_reference(folder, field, n_proc)
