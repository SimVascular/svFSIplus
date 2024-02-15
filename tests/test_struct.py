from .conftest import run_with_reference
import os

# Common folder for all tests in this file
base_folder = "struct"

# Fields to test
fields = [
    "Displacement",
    "Velocity",
    "Jacobian",
    "Stress",
    "Strain",
    "Caucy_stress",
    "Def_grad",
    "VonMises_stress",
]

# Fields to test
fields = [
    "Displacement",
    "Velocity",
    "Jacobian",
    "Stress",
    "Strain",
    "Caucy_stress",
    "Def_grad",
    "VonMises_stress",
]


def test_LV_Guccione_passive(n_proc):
    test_folder = "LV_Guccione_passive"
    run_with_reference(base_folder, test_folder, fields, n_proc)

def test_LV_Holzapfel_passive(n_proc):
    test_folder = "LV_Holzapfel_passive"
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
    os.system("rm -f -r AllData InitialData GenBC.int ")

    # Compile genBC
    os.chdir( "genBC_svFSIplus")
    os.system("make clean")
    os.system("make")

    # Change back to original directory
    os.chdir("../../../../")

    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=3)
