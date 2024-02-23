import pytest

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "shell"

# Fields to test
fields = ["Displacement", "Velocity", "Stress", "Strain"]


@pytest.mark.parametrize("n_proc", [1])
def test_plate(n_proc):
    test_folder = "plate"
    t_max = 10
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)


@pytest.mark.parametrize("n_proc", [1])
def test_valve(n_proc):
    test_folder = "valve"
    run_with_reference(base_folder, test_folder, fields, n_proc)
