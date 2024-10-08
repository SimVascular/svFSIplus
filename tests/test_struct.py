from .conftest import run_with_reference
import os
import subprocess

# Common folder for all tests in this file
base_folder = "struct"

# Fields to test
fields = [
    "Cauchy_stress",
    "Def_grad",
    "Displacement",
    "Jacobian",
    "Jacobian",
    "Stress",
    "Strain",
    "Velocity",
    "VonMises_stress",
]


def test_LV_Guccione_passive(n_proc):
    test_folder = "LV_Guccione_passive"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_LV_HolzapfelOgden_passive(n_proc):
    test_folder = "LV_HolzapfelOgden_passive"
    run_with_reference(base_folder, test_folder, fields, n_proc)

def test_LV_HolzapfelOgdenModifiedAnisotropy_passive(n_proc):
    test_folder = "LV_HolzapfelOgdenModifiedAnisotropy_passive"
    run_with_reference(base_folder, test_folder, fields, n_proc)

def test_block_compression(n_proc):
    test_folder = "block_compression"
    run_with_reference(base_folder, test_folder, fields, n_proc)

def test_robin(n_proc):
    test_folder = "robin"
    run_with_reference(base_folder, test_folder, fields, n_proc)

def test_LV_NeoHookean_passive(n_proc):
    test_folder = "LV_NeoHookean_passive"
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=5)
    
def test_LV_NeoHookean_passive_genBC(n_proc):
    test_folder = "LV_NeoHookean_passive_genBC"

    # Remove old genBC output
    os.chdir(os.path.join("cases", base_folder, test_folder))
    for name in ["AllData", "InitialData", "GenBC.int"]:
        if os.path.isfile(name):
            os.remove(name)

    # Compile genBC
    os.chdir("genBC_svFSIplus")
    subprocess.run(["make", "clean"], check=True)
    subprocess.run(["make"], check=True)

    # Change back to original directory
    os.chdir("../../../../")

    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=3)

def test_LV_NeoHookean_passive_sv0D(n_proc):
    test_folder = "LV_NeoHookean_passive_sv0D"

    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=3)

def test_tensile_adventitia_Guccione_active(n_proc):
    test_folder = "tensile_adventitia_Guccione_active"
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=2)

def test_tensile_adventitia_Newtonian_viscosity(n_proc):
    test_folder = "tensile_adventitia_Newtonian_viscosity"
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=1)

def test_tensile_adventitia_Potential_viscosity(n_proc):
    test_folder = "tensile_adventitia_Potential_viscosity"
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=1)
