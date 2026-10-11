import platform,json,os
import numpy,scipy,matplotlib,pytest,threadpoolctl
from common import save
workers=min(8,max(1,(os.cpu_count() or 2)//2))
save('results/environment.json',{'python':platform.python_version(),'platform':platform.platform(),'numpy':numpy.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'pytest':pytest.__version__,'threadpoolctl':threadpoolctl.__version__,'fft_workers':1,'simulation_processes':{'calibration':4,'test':workers,'controls':workers,'maximum_concurrent_test_and_control_processes':2*workers},'execution_note':'The already-running calibration used four processes. Later test/control batches used eight on this host, and independent control trajectories ran alongside test trajectories. Solver operations, seeds and model freeze were unchanged.','arrays':'float64/complex128','observation_units':'Declared model length/time units; target integral and RMSE have model-time units'})
