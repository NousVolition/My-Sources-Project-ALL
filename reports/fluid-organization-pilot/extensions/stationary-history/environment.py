import sys,platform,os,importlib.metadata
from common import save,now
save('results/environment.json',dict(recorded_utc=now(),python=sys.version,platform=platform.platform(),logical_cpus=os.cpu_count(),versions={s:importlib.metadata.version(s) for s in ['numpy','scipy','numba','llvmlite','matplotlib','pytest','threadpoolctl']},precision='float64 and complex128',fft_workers_per_run=1,compiled_kernels='Numba, fastmath disabled, serial reductions',calibration_workers=8))
