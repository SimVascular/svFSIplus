import json
import os
import pdb
import subprocess
import meshio
import pdb
import pytest

import numpy as np
import pandas as pd

this_file_dir = os.path.abspath(os.path.dirname(__file__))
cpp_exec = os.path.join(this_file_dir, "..", "build", "svFSI-build", "bin", "svFSI")

# relative tolerances for tested results
RTOL = {'Pressure': 1.0e-12, 'Velocity': 1.0e-12, 'Action_potential': 1.0e-12, 'Temperature': 1.0e-12, 'ECG': 1.0e-12, 'Displacement': 1.0e-12, \
        'Stress': 1.0e-12, 'VonMises_stress': 1.0e-12, 'Cauchy_stress' : 1.0e-12, 'Strain' : 1.0e-12, 'Jacobian' : 1.0e-12, 'Divergence' : 1.0e-12, \
        'Def_grad': 1.0e-12, 'Traction':1.0e-12, 'WSS':1.0e-12}

# number of processors to test
procs = [1, 3, 4]


def run_by_name(folder, name, t_max, n_proc=1):
    """
    Run a test case and return results
    Args:
        folder: location from which test will be executed
        name: name of svFSIplus input file (.xml)
        t_max: time step to compare
        n_proc: number of processors

    Returns:
    Simulation results
    """
    # run simulation
    cmd = " ".join(["mpirun", "--oversubscribe" if n_proc>1 else "", "-np", str(n_proc), cpp_exec, name])
    subprocess.call(cmd, cwd=folder, shell=True)

    # read results
    fname = os.path.join(folder, str(n_proc) + "-procs", "result_" + str(t_max).zfill(3) + "_cpp.vtu")
    assert os.path.exists(fname), "no svFSIplus output: " + fname
    return meshio.read(fname)


def run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc=1):
    """
    Run a test case and compare it to a stored reference solution
    Args:
        folder: location from which test will be executed
        name_inp: name of svFSIplus input file (.xml)
        name_ref: name of refence file (.vtu)
        fields: array fields to compare (e.g. ["Pressure", "Velocity"])
        t_max: time step to compare
        n_proc: number of processors
    """
    # run simulation
    res = run_by_name(folder, name_inp, t_max, n_proc)

    # read reference
    fname = os.path.join(folder, name_ref)
    ref = meshio.read(fname)

    # check results
    for f in fields:
        a = res.point_data[f]
        b = ref.point_data[f]

        # truncate last dimension if solution is 2D but reference is 3D
        if len(a.shape) == 2:
            if a.shape[1] == 2 and b.shape[1] == 3:
                assert not np.any(b[:, 2])
                b = b[:, :2]

        # compare solution to reference
        close = np.isclose(a, b, rtol=RTOL[f])
        if np.all(close):
            return
        else:
            msg = "Test failed!"
            msg += "\nResults in field " + f + " differ by more than rtol=" + str(RTOL[f])
            msg += " in " + str(np.sum(close)) + " out of " + str(close.size) + " results."
            msg += " Max. abs. difference is " + "{:.1e}".format(np.max(np.abs(a-b)))
            raise ValueError(msg)


@pytest.mark.parametrize("mesh", ["N" + str(2**i).zfill(3) for i in range(2, 3)])
@pytest.mark.parametrize("ele", ["P1P1"])
@pytest.mark.parametrize("n_proc", procs)
def test_stokes_manufactured_solution(ele, mesh, n_proc):
    folder = os.path.join("cases", "stokes_manufactured_solution", ele, mesh)
    fields = ["Pressure", "Velocity"]
    t_max = {"P1P1": 250, "P2P1": 50}
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max[ele]).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max[ele], n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_1Dcable_TTP(n_proc):
    folder = os.path.join("cases", "1Dcable_TTP")
    field = ["Action_potential"]
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_2Dspiral_BO(n_proc):
    folder = os.path.join("cases", "2Dspiral_BO")
    field = ["Action_potential"]
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_2Dsquare_AP(n_proc):
    folder = os.path.join("cases", "2Dsquare_AP")
    field = ["Action_potential"]
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_purkinje(n_proc):
    folder = os.path.join("cases", "purkinje")
    field = ["Action_potential"]
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)

@pytest.mark.parametrize("confs_ecgs", [["BICG_CN_myocardium_BO"  , -0.0781125,  0.0781125,  0.00885273],
                                        ["CG_RK4_endocardium_BO"  , -0.0780188,  0.0780188,  0.00884210],
                                        ["GMRES_FE_epicardium_TTP", -0.0786707,  0.0786707,  0.00891599],
                                        ["GMRES_FE_pfib_AP"       ,  0.0786707, -0.0786707, -0.00891599]])
@pytest.mark.parametrize("n_proc", procs)
def test_niederer_benchmark_ECGs_quadrature(confs_ecgs, n_proc):
    folder = os.path.join("cases", "niederer_benchmark_ECGs_quadrature")
    field = ["Action_potential"]
    t_max = 1
    name_inp = "svFSI_" + confs_ecgs[0] + ".xml"
    name_ref = "result_" + confs_ecgs[0] + "_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)

    for jj in range(0, 3):
        ecg_trace = pd.read_csv(folder + "/" + str(n_proc) + "-procs/ecglead_" + str(jj + 1) + ".txt", header = None)
        assert abs((ecg_trace.iloc[-1, 1] - confs_ecgs[jj + 1]) / confs_ecgs[jj + 1]) < RTOL['ECG'], \
               "Results in field ecglead_" + str(jj + 1) + ".txt differ by more than rtol=" + str(RTOL['ECG']) + " for test case " + confs_ecgs[0] 

@pytest.mark.parametrize("name_inp", ["svFSI_CG.xml", "svFSI_BICG.xml", "svFSI_GMRES.xml"])
@pytest.mark.parametrize("n_proc", procs)
def test_diffusion_line_source(name_inp, n_proc):
    folder = os.path.join("cases", "diffusion_line_source")
    field = ["Temperature"]
    t_max = 2
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)


@pytest.mark.parametrize("n_proc", procs)
def test_pipe3D_RCR(n_proc):
    folder = os.path.join("cases", "pipe3D_RCR")
    fields = ["Pressure", "Velocity"]
    t_max = 2
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)

    
@pytest.mark.parametrize("n_proc", procs)
def test_cavity_2d(n_proc):
    folder = os.path.join("cases", "driven_cavity_2D")
    fields = ["Pressure", "Velocity"]
    t_max = 2
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)


@pytest.mark.parametrize("n_proc", procs)
def test_ale_3d_pipe(n_proc):
    folder = os.path.join("cases", "ale_3d_pipe")
    fields = ["Displacement", "Pressure", "Velocity"]
    t_max = 5
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)

@pytest.mark.parametrize("ele", ["P1P1_VMS"])
@pytest.mark.parametrize("n_proc", procs)
def test_block_compression_ustruct(ele, n_proc):
    folder = os.path.join("cases", "block_compression_ustruct", ele)
    field = ["Displacement", "Pressure", "Stress", "Divergence"]
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_tensile_adventitia_HGO(n_proc):
    folder = os.path.join("cases", "tensile_adventitia_HGO")
    field = ["Displacement", "Velocity", "Stress", "VonMises_stress"]
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_LV_Guccione_active(n_proc):
    folder = os.path.join("cases", "LV_Guccione_active")
    field = ["Displacement", "Velocity", "Pressure", "VonMises_stress", "Cauchy_stress", "Strain", "Jacobian"]
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_LV_Guccione_passive(n_proc):
    folder = os.path.join("cases", "LV_Guccione_passive")
    fields = ["Displacement", "Velocity", "Jacobian"] 
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_block_compression_struct(n_proc):
    folder = os.path.join("cases", "block_compression_struct")
    fields = ["Displacement", "Velocity", "Jacobian", "Stress", "Strain", "Caucy_stress", "Def_grad", "VonMises_stress"] 
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_dye_AD(n_proc):
    folder = os.path.join("cases", "dye_AD")
    fields = ["Velocity", "Pressure", "Traction", "WSS"] 
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_newtonian_flow(n_proc):
    folder = os.path.join("cases", "newtonian_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"] 
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_casson_flow(n_proc):
    folder = os.path.join("cases", "casson_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"] 
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_carreau_yasuda_flow(n_proc):
    folder = os.path.join("cases", "carreau_yasuda_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"] 
    t_max = 1
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max, n_proc)

@pytest.mark.parametrize("n_proc", procs)
def test_cmm_3d_pipe(n_proc):
    folder = os.path.join("cases", "cmm_3d_pipe")
    inflate_folder = os.path.join(folder,"2a-inflate")
    fields = ["Displacement"]
    t_max = 3
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(inflate_folder, name_inp, name_ref, fields, t_max, 1)

    inflate_cmm_folder = os.path.join(folder,"3a-inflate-cmm")
    fields = ["Displacement", "Pressure", "Velocity"]
    t_max = 5
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(inflate_cmm_folder, name_inp, name_ref, fields, t_max, n_proc)

    prestress_folder = os.path.join(folder,"2b-prestress")
    fields = ["Stress"]
    t_max = 3
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(prestress_folder, name_inp, name_ref, fields, t_max, 1)

    prestress_cmm_folder = os.path.join(folder,"3b-prestress-cmm")
    fields = ["Displacement", "Stress", "Pressure", "Velocity"]
    t_max = 5
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(prestress_cmm_folder, name_inp, name_ref, fields, t_max, n_proc)   


