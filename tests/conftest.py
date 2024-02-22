import numpy as np

import pytest
import os
import subprocess
import meshio

this_file_dir = os.path.abspath(os.path.dirname(__file__))
cpp_exec = os.path.join(this_file_dir, "..", "build", "svFSI-build", "bin", "svFSI")

# Relative tolerances for each tested field
RTOL = {
    "Action_potential": 1.0e-10,
    "Cauchy_stress": 1.0e-4,
    "Def_grad": 1.0e-10,
    "Divergence": 1.0e-9,
    "Displacement": 1.0e-10,
    "Jacobian": 1.0e-10,
    "Pressure": 1.0e-6,
    "Stress": 1.0e-4,
    "Strain": 1.0e-10,
    "Temperature": 1.0e-10,
    "Traction": 1.0e-6,
    "Velocity": 1.0e-7,
    "VonMises_stress": 1.0e-3,
    "Vorticity": 1.0e-7,
    "WSS": 1.0e-8,
}

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
    if not os.path.exists(fname):
        raise RuntimeError("No svFSIplus output: " + fname)
    return meshio.read(fname)


def run_with_reference(
    base_folder,
    test_folder,
    fields,
    n_proc=1,
    t_max=1,
    name_ref=None,
    name_inp="svFSI.xml",
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
    folder = os.path.join("cases", base_folder, test_folder)
    res = run_by_name(folder, name_inp, t_max, n_proc)

    # read reference
    fname = os.path.join(folder, name_ref)
    ref = meshio.read(fname)

    # check results
    msg = ""
    for f in fields:
        # extract field
        if f not in res.point_data.keys():
            raise ValueError("Field " + f + " not in simulation result")
        a = res.point_data[f]

        if f not in ref.point_data.keys():
            raise ValueError("Field " + f + " not in reference result")
        b = ref.point_data[f]

        # truncate last dimension if solution is 2D but reference is 3D
        if len(a.shape) == 2:
            if a.shape[1] == 2 and b.shape[1] == 3:
                assert not np.any(b[:, 2])
                b = b[:, :2]

        # pick tolerance for current field
        if f not in RTOL:
            raise ValueError("No tolerance defined for field " + f)
        rtol = RTOL[f]

        # relative difference (as computed in np.isclose)
        # note that we consider rtol as absolute zero (and as relative tolerance)
        a_fl = a.flatten()
        b_fl = b.flatten()
        rel_diff = np.abs(a_fl - b_fl) - rtol - rtol * np.abs(b_fl)

        # throw error if not all results are within relative tolerance
        close = rel_diff <= 0.0
        if not np.all(close):
            # portion of individual results that are above the tolerance
            wrong = 1 - np.sum(close) / close.size

            # location of maximum relative difference
            i_max = rel_diff.argmax()

            # maximum relative difference
            max_rel = rel_diff[i_max]

            # maximum absolute difference at same location
            max_abs = np.abs(a_fl[i_max] - b_fl[i_max])

            # throw error message for pytest
            msg += "Test failed in field " + f + "."
            msg += " Results differ by more than rtol=" + str(rtol)
            msg += " in {:.1%}".format(wrong)
            msg += " of results."
            msg += " Max. rel. difference is"
            msg += " {:.1e}".format(max_rel)
            msg += " (abs. {:.1e}".format(max_abs) + ")\n"
    # check all fields first and then throw error if any failed
    if msg:
        raise AssertionError(msg)
