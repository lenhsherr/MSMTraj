# Import package, test suite, and other packages as needed
import sys
from pathlib import Path

import pytest

from  FOAM.trajectory import Trajectory

DATA_DIR = 'FOAM/data/'

def test_load_existing_h5():
    config = {'coordinates_path': DATA_DIR + 'alanine-dipeptide-0-250ns-nowater.h5'}
    traj = Trajectory(config)
    assert traj.coordinates_path is not None

def test_read_existing_h5():
    config = {'coordinates_path': DATA_DIR + 'alanine-dipeptide-0-250ns-nowater.h5'}
    traj = Trajectory(config)
    traj.read()
    assert traj.coordinates is not None