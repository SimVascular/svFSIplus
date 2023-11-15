import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "stokes")


@pytest.mark.parametrize("mesh", ["N" + str(2**i).zfill(3) for i in range(2, 3)])
@pytest.mark.parametrize("ele", ["P1P1"])
def test_manufactured_solution(ele, mesh, n_proc):
    folder = os.path.join(base_folder, "manufactured_solution", ele, mesh)
    fields = ["Pressure", "Velocity"]
    t_max = {"P1P1": 250, "P2P1": 50}
    run_with_reference(folder, fields, n_proc, t_max[ele])
