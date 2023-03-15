from tempfile import TemporaryDirectory

import numpy as np
import pytest


@pytest.fixture
def tempdir():
    """Temporary directory for test purposes."""
    with TemporaryDirectory() as tempdir:
        yield tempdir
