"""Display-only meshes and slices from retained velocity fields. No evolution."""
import base64
import hashlib
import io
import json
from pathlib import Path
import numpy as np
from scipy import ndimage
from PIL import Image
from matplotlib import colormaps
from solver import Solver
from protocol import ROOT, CASES, atomic_json

VERSION = 1
DISPLAY_N = 64
CACHE = ROOT / 'visual-cache'


def encoded(a):
    return base64.b64encode(a.tobytes()).decode('ascii')


def surface(field, level):
    """Marching tetrahedra on a periodic display lattice, including box faces."""
    n = field.shape[0]
    f = np.pad(field, ((0,1),)*3, mode='wrap')
    corners = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                        [0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    values = np.stack([f[x:x+n,y:y+n,z:z+n].ravel() for x,y,z in corners], axis=1)
    active = (values.min(axis=1) < level) & (values.max(axis=1) >= level)
    ids = np.flatnonzero(active)
    base = np.array(np.unravel_index(ids, (n,n,n))).T
    vertices = base[:,None,:] + corners[None,:,:]
    values = values[active]
    triangles = []
    for tet in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
        val = values[:,tet]; pos = vertices[:,tet,:]
        masks = np.sum((val >= level) * (1 << np.arange(4)), axis=1)
        for mask in range(1,15):
            selected = masks == mask
            if not selected.any(): continue
            vv=val[selected]; pp=pos[selected]
            ins=[j for j in range(4) if mask & (1<<j)]
            outs=[j for j in range(4) if not mask & (1<<j)]
            def edge(a,b):
                t=(level-vv[:,a])/(vv[:,b]-vv[:,a])
                return pp[:,a]+t[:,None]*(pp[:,b]-pp[:,a])
            if len(ins) == 1:
                triangles.append(np.stack([edge(ins[0],o) for o in outs],axis=1))
            elif len(outs) == 1:
                triangles.append(np.stack([edge(i,outs[0]) for i in ins],axis=1))
            else:
                a,b=ins; c,d=outs
                p,q,r,s=edge(a,c),edge(a,d),edge(b,c),edge(b,d)
                triangles.extend((np.stack((p,q,r),axis=1),np.stack((q,s,r),axis=1)))
    if not triangles:return {'positions':'','normals':'','count':0,'level':level}
    pts=np.concatenate(triangles).reshape(-1,3)
    grad=np.stack([(np.roll(field,-1,j)-np.roll(field,1,j))/2 for j in range(3)])
    normals=-np.stack([ndimage.map_coordinates(g,pts.T,order=1,mode='grid-wrap') for g in grad],axis=1)
    normals/=np.maximum(np.linalg.norm(normals,axis=1,keepdims=True),1e-12)
    coords=pts*6/n-3
    return {'positions':encoded(np.rint(coords/3*32767).astype('<i2')),
            'normals':encoded(np.rint(normals*127).astype('i1')),
            'count':len(pts),'level':level}


def slice_image(mag, plane):
    n=mag.shape[0]
    a=mag[:,:,n//2].T if plane=='xy' else mag[n//2,:,:].T
    # Native-grid sample, fixed physical range across every frame and case.
    rgba=np.uint8(colormaps['magma'](np.clip(a/200,0,1))*255)
    out=io.BytesIO(); Image.fromarray(rgba[::-1]).save(out,format='PNG')
    return 'data:image/png;base64,'+base64.b64encode(out.getvalue()).decode('ascii')


def prepare(case, records, proto):
    # A single resolution per case: never splice different grids into one movie.
    available=[]
    for n in (48,64,80,112,160):
        rid=f'{case}-fourier-n{n}-base'
        if rid in records and records[rid].get('rows') and (ROOT/'runs'/rid/'field-t0.00.npz').exists():
            available.append((records[rid]['status']=='complete',n,rid))
    if available:
        _,n,rid=max(available)
        # result.json is atomically published after its field copy finishes.
        # Ignore a newer copy that the running solver may still be writing.
        saved_t=records[rid]['rows'][-1]['t']
        paths=[p for p in sorted((ROOT/'runs'/rid).glob('field-t*.npz'))
               if float(p.stem.removeprefix('field-t')) <= saved_t+1e-9]
        frames=[]
        for p in paths:
            stamp=f'{VERSION}:{p.stat().st_mtime_ns}:{p.stat().st_size}'
            cp=CACHE/(rid+'-'+p.stem+'.json')
            cached=json.loads(cp.read_text()) if cp.exists() else None
            if not cached or cached.get('stamp')!=stamp:
                with np.load(p,allow_pickle=False) as saved:
                    h=saved['h']; t=float(saved['t'])
                cached=make_frame(h,n,t)
                cached.update(stamp=stamp,source=str(p.relative_to(ROOT)),
                              source_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                atomic_json(cp,cached)
            frame=dict(cached)
            row=next((r for r in records[rid]['rows'] if abs(r['t']-frame['t'])<1e-8),None)
            if row:
                frame.update(Wmax=row['Wmax'],I=row['I'],energy=row['energy'])
            frames.append(frame)
        return {'case':case,'n':n,'run':rid,'initial_only':False,'frames':frames}
    n=80
    cp=CACHE/f'{case}-initial-v{VERSION}.json'
    if cp.exists(): frame=json.loads(cp.read_text())
    else:
        solver=Solver(n,workers=2)
        h,_,_=solver.initial(case,proto['chosen_parameters']['outer_gamma'])
        frame=make_frame(h,n,0)
        frame.update(source='Initial field reconstructed with the unchanged design-2 constructor.',I=0,
                     energy=solver.inner(h,h)/2)
        atomic_json(cp,frame)
    return {'case':case,'n':n,'run':'Initial construction only','initial_only':True,'frames':[frame]}


def make_frame(h,n,t):
    solver=Solver(n,workers=2)
    omega=solver.real(solver.curl(h,physical=True))
    mag=np.sqrt(np.sum(omega*omega,axis=0))
    del omega
    q=np.arange(DISPLAY_N)*n/DISPLAY_N
    coords=np.array(np.meshgrid(q,q,q,indexing='ij'))
    display=ndimage.map_coordinates(mag,coords,order=1,mode='grid-wrap')
    return {'t':t,'Wmax':float(mag.max()),'meshes':[surface(display,v) for v in (8.,40.)],
            'xy':slice_image(mag,'xy'),'yz':slice_image(mag,'yz')}


def viewer_html(records,proto):
    CACHE.mkdir(exist_ok=True)
    data={'display_n':DISPLAY_N,'cases':[prepare(c,records,proto) for c in CASES]}
    template=(ROOT/'vortex-viewer.html').read_text(encoding='utf-8')
    return template.replace('__VORTEX_DATA__',json.dumps(data,separators=(',',':')).replace('</','<\\/'))
