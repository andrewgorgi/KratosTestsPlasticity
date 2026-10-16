"""convergence del caso gap_sbm. Eseguibile da qualunque directory."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from study import main

if __name__ == "__main__":
    main("gap_sbm", "convergence")
