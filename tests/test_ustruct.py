import os
import pytest

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "ustruct"

# Fields to test
fields = [
    "Displacement",
    "Velocity",
    "Pressure",
    "VonMises_stress",
    "Cauchy_stress",
    "Strain",
    "Jacobian",
]


@pytest.mark.parametrize("ele", ["P1P1_VMS"])
def test_block_compression(ele, n_proc):
    test_folder = os.path.join("block_compression", ele)
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_tensile_adventitia_HGO(n_proc):
    test_folder = "tensile_adventitia_HGO"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_LV_Guccione_active(n_proc):
    test_folder = "LV_Guccione_active"
    run_with_reference(base_folder, test_folder, fields, n_proc)
