import os
import pytest

import pandas as pd

from .conftest import run_with_reference, DEFAULT_TOL

# Common folder for all tests in this file
base_folder = os.path.join("cases", "cep")


def test_cable_TTP_1d(n_proc):
    folder = os.path.join(base_folder, "cable_TTP_1d")
    field = ["Action_potential"]
    run_with_reference(folder, field, n_proc)


def test_spiral_BO_2d(n_proc):
    folder = os.path.join(base_folder, "spiral_BO_2d")
    field = ["Action_potential"]
    run_with_reference(folder, field, n_proc)


def test_square_AP_2d(n_proc):
    folder = os.path.join(base_folder, "square_AP_2d")
    field = ["Action_potential"]
    run_with_reference(folder, field, n_proc)


def test_purkinje(n_proc):
    folder = os.path.join(base_folder, "purkinje")
    field = ["Action_potential"]
    run_with_reference(folder, field, n_proc)


@pytest.mark.parametrize(
    "confs_ecgs",
    [
        ["BICG_CN_epicardium_BO", -0.0786707, 0.0786707, 0.00891599],
        ["CG_RK4_myocardium_BO", -0.0781115, 0.0781115, 0.00885261],
        ["GMRES_FE_epicardium_TTP", -0.0786707, 0.0786707, 0.00891599],
        ["GMRES_FE_pfib_AP", 0.0786707, -0.0786707, -0.00891599],
    ],
)
def test_niederer_benchmark_ECGs_quadrature(confs_ecgs, n_proc):
    folder = os.path.join(base_folder, "niederer_benchmark_ECGs_quadrature")
    field = ["Action_potential"]
    t_max = 1
    name_inp = "svFSI_" + confs_ecgs[0] + ".xml"
    name_ref = "result_" + confs_ecgs[0] + "_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, field, n_proc, t_max, name_ref, name_inp)

    for jj in range(0, 3):
        ecg_trace = pd.read_csv(
            folder + "/" + str(n_proc) + "-procs/ecglead_" + str(jj + 1) + ".txt",
            header=None,
        )
        assert (
            abs((ecg_trace.iloc[-1, 1] - confs_ecgs[jj + 1]) / confs_ecgs[jj + 1])
            < DEFAULT_TOL
        ), (
            "Results in field ecglead_"
            + str(jj + 1)
            + ".txt differ by more than rtol="
            + str(DEFAULT_TOL)
            + " for test case "
            + confs_ecgs[0]
        )
