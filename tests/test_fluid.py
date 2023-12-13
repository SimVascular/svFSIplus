from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "fluid"

# Fields to test
fields = ["Velocity", "Pressure", "Traction", "WSS"]


def test_pipe_RCR_3d(n_proc):
    test_folder = "pipe_RCR_3d"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)


def test_cavity_2d(n_proc):
    test_folder = "driven_cavity_2d"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)


def test_dye_AD(n_proc):
    test_folder = "dye_AD"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_newtonian(n_proc):
    test_folder = "newtonian"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_casson(n_proc):
    test_folder = "casson"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_carreau_yasuda(n_proc):
    test_folder = "carreau_yasuda"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_iliac_artery(n_proc):
    test_folder = "iliac_artery"
    run_with_reference(base_folder, test_folder, fields, n_proc)
