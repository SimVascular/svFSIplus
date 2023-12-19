import os

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "cmm"

# Fields to test
fields = ["Stress", "Displacement",  "Pressure", "Velocity", "Traction", "WSS"]


def test_pipe_3d(n_proc):
    folder = "pipe_3d"
    inflate_folder = os.path.join(folder, "2a-inflate")
    t_max = 3
    run_with_reference(base_folder, inflate_folder, fields[1:2], 1, t_max)

    inflate_cmm_folder = os.path.join(folder, "3a-inflate-cmm")
    t_max = 5
    run_with_reference(base_folder, inflate_cmm_folder, fields[1::], n_proc, t_max)

    prestress_folder = os.path.join(folder, "2b-prestress")
    t_max = 3
    run_with_reference(base_folder, prestress_folder, fields[0:2], 1, t_max)

    prestress_cmm_folder = os.path.join(folder, "3b-prestress-cmm")
    t_max = 5
    run_with_reference(base_folder, prestress_cmm_folder, fields[1::], n_proc, t_max)



def test_iliac_artery_variable_wall_props(n_proc):
    folder = "iliac_artery_variable_wall_props"
    inflate_folder = os.path.join(folder, "1-rigid-solution")
    t_max = 3
    run_with_reference(base_folder, inflate_folder, fields[2::], n_proc, t_max)

    inflate_cmm_folder = os.path.join(folder, "2-inflate")
    t_max = 3
    run_with_reference(base_folder, inflate_cmm_folder, fields[1:2], 1, t_max)

    prestress_cmm_folder = os.path.join(folder, "3-inflate-cmm")
    t_max = 3
    run_with_reference(base_folder, prestress_cmm_folder, fields[1::], n_proc, t_max)

