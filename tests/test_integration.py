import json
import os
import pdb
import subprocess
import meshio
import pdb
import pytest

import numpy as np

this_file_dir = os.path.abspath(os.path.dirname(__file__))
cpp_exec = os.path.join(this_file_dir, "..", "build", "svFSI-build", "bin", "svFSI")
# todo: add second executable for "classic" svFSI and compare results

RTOL = {'Pressure': 1.0e-12, 'Velocity': 1.0e-12, 'Action_potential': 1.0e-12, 'Temperature': 1.0e-12}


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
    # subprocess.call(["mpirun -np " + str(n_proc), cpp_exec, name], cwd=folder)
    subprocess.call([cpp_exec, name], cwd=folder)

    # read results
    fname = os.path.join(folder, str(n_proc) + "-procs", "result_" + str(t_max).zfill(3) + "_cpp.vtu")
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
    res = run_by_name(folder, name_inp, t_max)

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
        assert np.all(np.isclose(a, b, rtol=RTOL[f]))


@pytest.mark.parametrize("mesh", ["N" + str(2**i).zfill(3) for i in range(2, 3)])
@pytest.mark.parametrize("ele", ["P1P1"])
def test_stokes_manufactured_solution(ele, mesh):
    folder = os.path.join("cases", "stokes_manufactured_solution", ele, mesh)
    fields = ["Pressure", "Velocity"]
    t_max = {"P1P1": 250, "P2P1": 50}
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max[ele]).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max[ele])

def test_niederer_benchmark():
    folder = os.path.join("cases", "niederer_benchmark")
    field = ["Action_potential"]
    t_max = 30
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max)

def test_diffusion_line_source():
    folder = os.path.join("cases", "diffusion_line_source")
    field = ["Temperature"]
    t_max = 20
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, field, t_max)
    
def test_pipe3D_RCR():
    folder = os.path.join("cases", "pipe3D_RCR")
    fields = ["Pressure", "Velocity"]
    t_max = 5
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max)

def test_cavity_2d():
    folder = os.path.join("cases", "driven_cavity_2D")
    fields = ["Pressure", "Velocity"]
    t_max = 10
    name_inp = "svFSI.xml"
    name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, name_inp, name_ref, fields, t_max)
