"""Two bounded refinement batches using the byte-identical original model."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,traceback,shutil,ctypes
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import numpy as np
import scipy,matplotlib,platform
from threadpoolctl import threadpool_info
from model import Model,cases,run
R=Path(__file__).resolve().parent
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def serial(v):return v.tolist() if isinstance(v,np.ndarray) else v.item()
def write(name,data):
    path=R/name;temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,indent=2,default=serial),encoding='utf-8');temp.replace(path)
def verify(p):
    for name,digest in read(R/'freeze.json')['files'].items():assert sha(R/name)==digest,name
    b=Path(p['refinement']['parent_directory']);parent=read(b/'protocol.json')
    for name,digest in p['refinement']['parent_file_hashes'].items():assert sha(b/name)==digest,name
    for key in p['refinement']['same_model_keys']:assert p[key]==parent[key],key
    assert sha(R/'model.py')==p['refinement']['parent_model_sha256']
def free_ram():
    class MemoryStatus(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(name,ctypes.c_ulonglong) for name in ('total_physical','available_physical','total_page','available_page','total_virtual','available_virtual','available_extended')]
    m=MemoryStatus();m.length=ctypes.sizeof(m);assert ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m));return m.available_physical/2**30
def main():
    p=read(R/'protocol.json');verify(p);cs=cases(p)
    assert len(cs)==45 and p['batches']==[dict(n=64,dt=.005),dict(n=64,dt=.0025)]
    assert read(R/'preflight.json')['passed'];pools=threadpool_info();assert all(x['num_threads']==1 for x in pools)
    with (R/'controller.lock').open('x') as out:json.dump(dict(pid=os.getpid(),created_utc=datetime.now(timezone.utc).isoformat()),out)
    (R/'results').mkdir(exist_ok=False)
    handoff=dict(status='running',pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),stdout='run.log',stderr='run.err',protocol='protocol.json',freeze='freeze.json',summary='results/summary.json',parent=str(Path(p['refinement']['parent_directory'])),no_automatic_retry=True)
    write('handoff.json',handoff);write('results/cases.json',cs)
    write('environment.json',dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,platform=platform.platform(),threadpools=pools,fft_workers=1))
    start=time.monotonic();completed=[]
    for batch in p['batches']:
        assert free_ram()>=p['resources']['minimum_free_RAM_GiB'],'Free RAM below entry requirement'
        assert shutil.disk_usage(R).free>2*2**30,'Insufficient free disk'
        n,dt=batch['n'],batch['dt'];label=f'n{n}-dt{dt:g}'
        def progress(t,bounds):write('results/progress.json',dict(status='running',pid=os.getpid(),batch=label,t=t,epoch=time.time(),completed=completed,bounds=bounds))
        progress(0,{})
        result,fields,times=run(p,cs,n,dt,progress)
        np.savez_compressed(R/'results'/f'{label}.npz',fields=fields,time=times)
        write(f'results/{label}.json',result);completed.append(label)
        write('results/progress.json',dict(status='batch_completed',pid=os.getpid(),batch=label,t=4,epoch=time.time(),completed=completed))
        print(label,'complete',result['elapsed_seconds'],'seconds',flush=True)
        assert time.monotonic()-start<=p['resources']['max_total_seconds'],'Total time budget exceeded'
        assert sum(f.stat().st_size for f in R.rglob('*') if f.is_file())/2**30<=p['resources']['max_storage_GiB'],'Storage budget exceeded'
        del fields,result
    verify(p)
    write('results/summary.json',dict(status='completed',batches=completed,conditions=45,trajectories=90,elapsed_seconds=time.monotonic()-start,parent_hashes_matched=True,interpretation='Await comparison of old and new fields; retain all failed screens.'))
    write('results/progress.json',dict(status='completed',pid=os.getpid(),epoch=time.time(),completed=completed))
    handoff['status']='integration_completed';handoff['next']='Analyze once after process exit; preserve previous spatial failures.';write('handoff.json',handoff)
    print('COMPLETE',time.monotonic()-start,flush=True)
if __name__=='__main__':
    try:main()
    except Exception:
        failure=dict(error=traceback.format_exc(),pid=os.getpid(),utc=datetime.now(timezone.utc).isoformat())
        if (R/'results').is_dir():write('results/failure.json',failure)
        print(failure['error'],flush=True);raise
