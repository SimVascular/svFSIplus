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

def test_LV_NeoHookean_passive_genBC(n_proc):
    test_folder = "LV_NeoHookean_passive_genBC"

    # Remove old genBC output
    os.chdir(os.path.join("cases", base_folder, test_folder))
    os.system("rm -f -r AllData InitialData GenBC.int ")

    # Compile genBC
    os.chdir( "genBC_svFSIplus")
    os.system("make clean")
    os.system("make")

    # Change back to original directory
    os.chdir("../../../../")

    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=3)
