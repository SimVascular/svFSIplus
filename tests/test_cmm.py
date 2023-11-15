import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "cmm")


def test_pipe_3d(n_proc):
    folder = os.path.join(base_folder, "pipe_3d")
    inflate_folder = os.path.join(folder, "2a-inflate")
    fields = ["Displacement"]
    t_max = 3
    run_with_reference(inflate_folder, fields, 1, t_max)

    inflate_cmm_folder = os.path.join(folder, "3a-inflate-cmm")
    fields = ["Displacement", "Pressure", "Velocity"]
    t_max = 5
    run_with_reference(inflate_cmm_folder, fields, n_proc, t_max)

    prestress_folder = os.path.join(folder, "2b-prestress")
    fields = ["Stress"]
    t_max = 3
    run_with_reference(prestress_folder, fields, 1, t_max)

    prestress_cmm_folder = os.path.join(folder, "3b-prestress-cmm")
    fields = ["Displacement", "Stress", "Pressure", "Velocity"]
    t_max = 5
    run_with_reference(prestress_cmm_folder, fields, n_proc, t_max)
