"""Run the supplied fluid scripts in separate output directories."""
from pathlib import Path
import subprocess
import sys

if __name__ == '__main__':
    here=Path(__file__).resolve().parent
    for name in ('fluid_break_to_2','fluid_sources','fluid_pitchfork','fluid_correction'):
        folder=here/'rerun'/name
        folder.mkdir(parents=True,exist_ok=True)
        subprocess.run([sys.executable,'-u',str(here/'source'/f'{name}.py')],cwd=folder,check=True)
