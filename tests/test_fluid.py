import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "fluid")


def test_pipe_RCR_3d(n_proc):
    folder = os.path.join(base_folder, "pipe_RCR_3d")
    fields = ["Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 2)


def test_cavity_2d(n_proc):
    folder = os.path.join(base_folder, "driven_cavity_2d")
    fields = ["Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 2)


def test_dye_AD(n_proc):
    folder = os.path.join(base_folder, "dye_AD")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_newtonian(n_proc):
    folder = os.path.join(base_folder, "newtonian")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_casson(n_proc):
    folder = os.path.join(base_folder, "casson")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)


def test_carreau_yasuda(n_proc):
    folder = os.path.join(base_folder, "carreau_yasuda")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)

def test_iliac_artery(n_proc):
    folder = os.path.join(base_folder, "iliac_artery")
    fields = ["Velocity", "Pressure", "Traction", "WSS"]
    run_with_reference(folder, fields, n_proc)