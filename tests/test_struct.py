import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "struct")

# Fields to test
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


def test_LV_Guccione_passive(n_proc):
    folder = os.path.join(base_folder, "LV_Guccione_passive")
    fields = ["Displacement", "Velocity", "Jacobian"]
    run_with_reference(folder, fields, n_proc)


def test_block_compression(n_proc):
    folder = os.path.join(base_folder, "block_compression")
    run_with_reference(folder, fields, n_proc)


def test_robin(n_proc):
    folder = os.path.join(base_folder, "robin")
    run_with_reference(folder, fields, n_proc)
