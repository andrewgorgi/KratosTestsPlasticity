"""Plot annular geometry, background mesh and element centers for Gap-SBM."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from annulus_study import main

if __name__ == '__main__':
    main('gap_sbm', 'plot')
