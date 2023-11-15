import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "heats")


@pytest.mark.parametrize(
    "name_inp", ["svFSI_CG.xml", "svFSI_BICG.xml", "svFSI_GMRES.xml"]
)
def test_diffusion_line_source(name_inp, n_proc):
    folder = os.path.join(base_folder, "diffusion_line_source")
    field = ["Temperature"]
    run_with_reference(folder, field, n_proc, 2, name_inp=name_inp)
