import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "fluid")


def test_pipe3D_RCR(n_proc):
    folder = os.path.join(base_folder, "pipe3D_RCR")
    fields = ["Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 2)


def test_cavity_2d(n_proc):
    folder = os.path.join(base_folder, "driven_cavity_2D")
    fields = ["Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 2)


def test_dye_AD(n_proc):
    folder = os.path.join(base_folder, "dye_AD")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_newtonian_flow(n_proc):
    folder = os.path.join(base_folder, "newtonian_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_casson_flow(n_proc):
    folder = os.path.join(base_folder, "casson_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_carreau_yasuda_flow(n_proc):
    folder = os.path.join(base_folder, "carreau_yasuda_flow")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)
