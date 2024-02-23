import os
import pytest

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "stokes"

# Fields to test
fields = ["Pressure", "Velocity", "Traction", "WSS", "Vorticity", "Divergence"]


@pytest.mark.parametrize("mesh", ["N" + str(2**i).zfill(3) for i in range(2, 3)])
@pytest.mark.parametrize("ele", ["P1P1"])
def test_manufactured_solution(ele, mesh, n_proc):
    test_folder = os.path.join("manufactured_solution", ele, mesh)
    t_max = {"P1P1": 250, "P2P1": 50}
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max[ele])
