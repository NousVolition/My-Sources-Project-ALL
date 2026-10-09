"""Downsample saved velocity directions for display. No molecular simulation."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
from scipy import ndimage,fft
ROOT=Path(__file__).resolve().parent
coords=np.array(np.meshgrid(np.linspace(-2.4,2.4,5),np.linspace(-2.4,2.4,5),np.linspace(-2.4,2.4,5),indexing='ij')).reshape(3,-1).T
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def make():
    records=[];sources=[]
    measured=json.loads((ROOT/'measurements.json').read_text())
    for kind,title,n in [('original','Original matched-stretch',64),('aligned','Aligned surroundings',160),('compressive','Compressive surroundings',160)]:
        frames=[]
        if kind=='original':
            folder=ROOT.parent/'runs/baseline-n64-base';paths=[folder/f'field-{i:03d}.npy' for i in (0,2,4,6,8,10,20,30,40)]
        else:
            folder=ROOT.parents[1]/'adversarial-vortex-study/runs'/f'{kind}-fourier-n160-base';paths=sorted(folder.glob('field-t*.npz'))
        for p in paths:
            if kind=='original':h=np.load(p,allow_pickle=False);t=int(p.stem.split('-')[1])*.01;markers=None
            else:
                with np.load(p,allow_pickle=False) as z:h=z['h'];t=float(z['t']);markers=z['markers'][:64]
            u=fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1),workers=1)
            sample=np.stack([ndimage.map_coordinates(ui,((coords+3)*n/6).T,order=1,mode='grid-wrap') for ui in u],axis=1)
            frame={'t':round(t,2),'velocity':np.round(sample,5).tolist()}
            if markers is not None:
                terms=min(measured['local_terms'][f'{kind}-fourier-n160-base'],key=lambda r:abs(r['t']-t))
                frame.update(markers=np.round(markers,5).tolist(),spin=np.round(terms['marker_spin'],4).tolist())
            frames.append(frame);sources.append({'case':kind,'file':p.name,'sha256':digest(p)})
        records.append({'id':kind,'title':title,'grid':n,'frames':frames,'arrow_scale':float(.66/max(np.linalg.norm(f['velocity'],axis=1).max() for f in frames)),
                        'final_winner':None if kind=='original' else measured['local_terms'][f'{kind}-fourier-n160-base'][-1]['marker_max_index']})
    data={'points':coords.tolist(),'cases':records,'source_hashes':sources,'external_force':0,'notes':'Velocity arrows at fixed sample sites; their movement is not a particle trajectory. Saved material dots, where available, are fluid labels, not individual molecules.'}
    (ROOT/'motion-data.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(json.dumps({'source_fields':len(sources),'bytes':(ROOT/'motion-data.json').stat().st_size}))
if __name__=='__main__':make()
