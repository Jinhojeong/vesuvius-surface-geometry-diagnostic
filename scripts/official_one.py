"""Isolated official score: argv = gt.npy pred.npy -> prints 4 floats.
Run in a subprocess so a Betti-matching C++ abort on one patch cannot kill the
parent run (see peel191.py). Score order: blend, topo, surface-dice, voi."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from official_score import official

gt = np.load(sys.argv[1]); pred = np.load(sys.argv[2])
s = official(gt.astype(bool), pred.astype(bool))
print(" ".join(f"{x:.6f}" for x in s))
