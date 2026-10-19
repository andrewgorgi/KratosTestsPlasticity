"""Cyclic annular torsion: gap_sbm, run_and_post."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cyclic_study import main

if __name__ == '__main__':
    main('gap_sbm', 'run_and_post')
