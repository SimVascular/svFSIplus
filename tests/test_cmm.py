import os

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "cmm"


def test_pipe_3d(n_proc):
    folder = "pipe_3d"
    inflate_folder = os.path.join(folder, "2a-inflate")
    fields = ["Displacement"]
    t_max = 3
    run_with_reference(base_folder, inflate_folder, fields, 1, t_max)

    inflate_cmm_folder = os.path.join(folder, "3a-inflate-cmm")
    fields = ["Displacement", "Pressure", "Velocity"]
    t_max = 5
    run_with_reference(base_folder, inflate_cmm_folder, fields, n_proc, t_max)

    prestress_folder = os.path.join(folder, "2b-prestress")
    fields = ["Stress"]
    t_max = 3
    run_with_reference(base_folder, prestress_folder, fields, 1, t_max)

    prestress_cmm_folder = os.path.join(folder, "3b-prestress-cmm")
    t_max = 5
    run_with_reference(base_folder, prestress_cmm_folder, fields, n_proc, t_max)
