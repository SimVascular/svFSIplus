import numpy as np

import pytest
import os
import subprocess
import meshio

this_file_dir = os.path.abspath(os.path.dirname(__file__))
cpp_exec = os.path.join(this_file_dir, "..", "build", "svFSI-build", "bin", "svFSI")

# Default relative tolerances for tested results
DEFAULT_TOL = 1.0e-12

# Dictionary with exceptions from DEFAULT_TOL
RTOL = {}

# Number of processors to test
PROCS = [1, 3, 4]


# Fixture to parametrize the number of processors for all tests
@pytest.fixture(params=PROCS)
def n_proc(request):
    return request.param


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
    cmd = " ".join(
        [
            "mpirun",
            "--oversubscribe" if n_proc > 1 else "",
            "-np",
            str(n_proc),
            cpp_exec,
            name,
        ]
    )
    subprocess.call(cmd, cwd=folder, shell=True)

    # read results
    fname = os.path.join(
        folder, str(n_proc) + "-procs", "result_" + str(t_max).zfill(3) + "_cpp.vtu"
    )
    assert os.path.exists(fname), "no svFSIplus output: " + fname
    return meshio.read(fname)


def run_with_reference(
    folder, fields, n_proc=1, t_max=1, name_ref=None, name_inp="svFSI.xml"
):
    """
    Run a test case and compare it to a stored reference solution
    Args:
        folder: location from which test will be executed
        fields: array fields to compare (e.g. ["Pressure", "Velocity"])
        n_proc: number of processors
        t_max: time step to compare
        name_inp: name of svFSIplus input file (.xml)
        name_ref: name of refence file (.vtu)
    """
    # default reference name
    if not name_ref:
        name_ref = "result_" + str(t_max).zfill(3) + ".vtu"

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
        if f in RTOL:
            rtol = RTOL[f]
        else:
            rtol = DEFAULT_TOL
        close = np.isclose(a, b, rtol=rtol)
        if np.all(close):
            return
        else:
            msg = "Test failed!"
            msg += (
                "\nResults in field " + f + " differ by more than rtol=" + str(RTOL[f])
            )
            msg += (
                " in " + str(np.sum(close)) + " out of " + str(close.size) + " results."
            )
            msg += " Max. abs. difference is " + "{:.1e}".format(np.max(np.abs(a - b)))
            raise ValueError(msg)
