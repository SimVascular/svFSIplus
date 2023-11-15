import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "fsi")


def test_pipe_3d(n_proc):
    folder = os.path.join(base_folder, "pipe_3d")
    fields = ["Displacement", "Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 5)
