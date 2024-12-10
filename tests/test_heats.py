import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "heats"

# Fields to test
fields = ["Temperature"]


@pytest.mark.parametrize("linear_solver", ["CG", "BICG", "GMRES"])
def test_diffusion_line_source(linear_solver, n_proc):
    test_folder = "diffusion_line_source"
    name_inp = "solver_" + linear_solver + ".xml"
    run_with_reference(base_folder, test_folder, fields, n_proc, 2, name_inp=name_inp)
