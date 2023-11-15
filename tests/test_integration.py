import os
import pytest

import pandas as pd

from .conftest import run_with_reference, DEFAULT_TOL


@pytest.mark.parametrize("mesh", ["N" + str(2**i).zfill(3) for i in range(2, 3)])
@pytest.mark.parametrize("ele", ["P1P1"])
def test_stokes_manufactured_solution(ele, mesh, n_proc):
    folder = os.path.join("cases", "stokes_manufactured_solution", ele, mesh)
    fields = ["Pressure", "Velocity"]
    t_max = {"P1P1": 250, "P2P1": 50}
    run_with_reference(folder, fields, n_proc, t_max[ele])


def test_1Dcable_TTP(n_proc):
    folder = os.path.join("cases", "1Dcable_TTP")
    field = ["Action_potential"]
    run_with_reference(folder, field, n_proc)


def test_2Dspiral_BO(n_proc):
    folder = os.path.join("cases", "2Dspiral_BO")
    field = ["Action_potential"]
    run_with_reference(folder, field, n_proc)


def test_2Dsquare_AP(n_proc):
    folder = os.path.join("cases", "2Dsquare_AP")
    field = ["Action_potential"]
    run_with_reference(folder, field, n_proc)


def test_purkinje(n_proc):
    folder = os.path.join("cases", "purkinje")
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
    folder = os.path.join("cases", "niederer_benchmark_ECGs_quadrature")
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


@pytest.mark.parametrize(
    "name_inp", ["svFSI_CG.xml", "svFSI_BICG.xml", "svFSI_GMRES.xml"]
)
def test_diffusion_line_source(name_inp, n_proc):
    folder = os.path.join("cases", "diffusion_line_source")
    field = ["Temperature"]
    run_with_reference(folder, field, n_proc, 2)


def test_pipe3D_RCR(n_proc):
    folder = os.path.join("cases", "pipe3D_RCR")
    fields = ["Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 2)


def test_cavity_2d(n_proc):
    folder = os.path.join("cases", "driven_cavity_2D")
    fields = ["Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 2)


def test_ale_3d_pipe(n_proc):
    folder = os.path.join("cases", "ale_3d_pipe")
    fields = ["Displacement", "Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 5)


@pytest.mark.parametrize("ele", ["P1P1_VMS"])
def test_block_compression_ustruct(ele, n_proc):
    folder = os.path.join("cases", "block_compression_ustruct", ele)
    field = ["Displacement", "Pressure", "Stress", "Divergence"]
    run_with_reference(folder, field, n_proc)


def test_tensile_adventitia_HGO(n_proc):
    folder = os.path.join("cases", "tensile_adventitia_HGO")
    field = ["Displacement", "Velocity", "Stress", "VonMises_stress"]
    run_with_reference(folder, field, n_proc)


def test_LV_Guccione_active(n_proc):
    folder = os.path.join("cases", "LV_Guccione_active")
    field = [
        "Displacement",
        "Velocity",
        "Pressure",
        "VonMises_stress",
        "Cauchy_stress",
        "Strain",
        "Jacobian",
    ]
    run_with_reference(folder, field, n_proc)


def test_LV_Guccione_passive(n_proc):
    folder = os.path.join("cases", "LV_Guccione_passive")
    fields = ["Displacement", "Velocity", "Jacobian"]
    run_with_reference(folder, fields, n_proc)


def test_block_compression_struct(n_proc):
    folder = os.path.join("cases", "block_compression_struct")
    fields = [
        "Displacement",
        "Velocity",
        "Jacobian",
        "Stress",
        "Strain",
        "Caucy_stress",
        "Def_grad",
        "VonMises_stress",
    ]
    run_with_reference(folder, fields, n_proc)


def test_dye_AD(n_proc):
    folder = os.path.join("cases", "dye_AD")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_newtonian_flow(n_proc):
    folder = os.path.join("cases", "newtonian_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_casson_flow(n_proc):
    folder = os.path.join("cases", "casson_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_carreau_yasuda_flow(n_proc):
    folder = os.path.join("cases", "carreau_yasuda_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_cmm_3d_pipe(n_proc):
    folder = os.path.join("cases", "cmm_3d_pipe")
    inflate_folder = os.path.join(folder, "2a-inflate")
    fields = ["Displacement"]
    t_max = 3
    run_with_reference(inflate_folder, fields, t_max, 1)

    inflate_cmm_folder = os.path.join(folder, "3a-inflate-cmm")
    fields = ["Displacement", "Pressure", "Velocity"]
    t_max = 5
    run_with_reference(inflate_cmm_folder, fields, n_proc, t_max)

    prestress_folder = os.path.join(folder, "2b-prestress")
    fields = ["Stress"]
    t_max = 3
    run_with_reference(prestress_folder, fields, t_max, 1)

    prestress_cmm_folder = os.path.join(folder, "3b-prestress-cmm")
    fields = ["Displacement", "Stress", "Pressure", "Velocity"]
    t_max = 5
    run_with_reference(prestress_cmm_folder, fields, n_proc, t_max)
