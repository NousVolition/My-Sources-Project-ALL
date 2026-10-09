"""Record the unchanged supplied breathing update at every fluid step.

Default source location is the parent folder when published under timeline/.
No measured quantity feeds back into the evolution.
"""
from pathlib import Path
import argparse, ast, hashlib, json, time
import numpy as np

def save_json(path, data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-dir',type=Path,default=here.parent)
    ap.add_argument('--out',type=Path,default=here/'rerun')
    args=ap.parse_args();source=args.source_dir.resolve();out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    if (out/'results.json').exists():raise RuntimeError('Output already complete; use a new output folder.')
    script=source/'breathing_rerun.py';tree=ast.parse(script.read_text(encoding='utf-8'))
    stop=next(i for i,node in enumerate(tree.body) if isinstance(node,ast.Assign)
              and any(isinstance(t,ast.Name) and t.id=='out' for t in node.targets))
    tree.body=tree.body[:stop]
    ns={'__file__':str(script),'__name__':'unchanged_breathing_initialization'}
    exec(compile(tree,str(script),'exec'),ns)
    m,ops,phi=ns['m'],ns['ops'],ns['phi'];dot=ns['dot'];mirror=ns['mirror_of']
    kx,ky,kz,k2,keep,dx=ops
    hashes={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in ('breathing_rerun.py','clay_hug.py')}
    assert hashes['breathing_rerun.py']=='af87b8d45b7b3760fb1929bd5d36522dabc4df12cc78206c0b7f692b5fba1a56'
    assert hashes['clay_hug.py']=='422b7160f9e96badfc68fb76958c3cb512982905b73b37f8db92e5fd37379c01'
    refpath=source/'reference-results.json'
    if not refpath.exists():refpath=source/'results.json'
    reference=json.loads(refpath.read_text())
    def measure(u):
        assert all(np.isfinite(v).all() for v in u)
        a=[v-w for v,w in zip(u,mirror(u))]
        U=float(np.sqrt(dot(u,u)));A=float(np.sqrt(dot(a,a)));D=dot(a,phi)
        remainder=[v-D*f for v,f in zip(a,phi)]
        B=float(np.sqrt(dot(remainder,remainder)))
        return dict(U=U,A=A,D=D,B=B,E=A/U,alignment=None if A<1e-12 else D/A,
                    decomposition_error=abs(A*A-D*D-B*B),energy=.5*m.L**3*U*U)
    all_results={};state={'status':'running','cases_complete':[],'source_hashes':hashes,'started_epoch':time.time()}
    save_json(out/'state.json',state)
    for name,u0 in ns['cases'].items():
        slug=name.replace(' ','-');folder=out/slug;folder.mkdir(exist_ok=True)
        print('START '+name,flush=True)
        u=m.project([v.copy() for v in u0],kx,ky,kz,k2,keep)
        p0=float(np.sqrt(sum(c*c for c in u)).max())
        dt=.04*dx/max(p0,1e-6);steps=int(np.ceil(.4/dt));dt=.4/steps
        zindex=ns['n']//2;plane_z=float(ns['z'][0,0,zindex]);series=[];views=[];field_hashes=[]
        for step in range(steps+1):
            kick={}
            if step:
                before=measure(u)
                gap=.15*np.sin(2*np.pi*step*dt/.4)
                kicked=[a+gap*b for a,b in zip(u,ns['breath_u'])]
                after=measure(kicked)
                kick=dict(kick_coefficient=float(gap),kick_A_change=after['A']-before['A'],
                          kick_D_change=after['D']-before['D'],kick_energy_change=after['energy']-before['energy'])
                u=m.advance(kicked,dt,ops)
            stats=measure(u);stats.update(t=step*dt,step=step,**kick)
            stats['divergence']=m.div_max(u,kx,ky,kz,keep)
            stats['speed_cfl']=float(np.sqrt(sum(v*v for v in u)).max()*dt/dx)
            series.append(stats)
            a=[v-w for v,w in zip(u,mirror(u))]
            views.append(np.stack([u[0][:,:,zindex],u[1][:,:,zindex],
                np.sqrt(sum(v*v for v in u))[:,:,zindex],np.sqrt(sum(v*v for v in a))[:,:,zindex]]))
            checkpoint=folder/f'field-{step:03d}.npz';np.savez_compressed(checkpoint,u=np.asarray(u))
            field_hashes.append({'file':checkpoint.name,'sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest()})
        endpoint_errors={k:max(abs(series[i][k]-reference[name][j][k]) for i,j in ((0,0),(-1,1))) for k in ('E','D')}
        assert max(endpoint_errors.values())<1e-12,(name,endpoint_errors)
        assert max(x['decomposition_error'] for x in series)<1e-12
        assert max(x['divergence'] for x in series)<1e-10
        result=dict(name=name,n=33,nu=.01,T=.4,dt=dt,steps=steps,plane_z=plane_z,series=series,
                    endpoint_absolute_errors=endpoint_errors,source_hashes=hashes)
        save_json(folder/'result.json',result);save_json(folder/'field-manifest.json',field_hashes)
        np.savez_compressed(folder/'slices.npz',views=np.asarray(views),axis=ns['x'][:,0,0])
        all_results[name]=result;state['cases_complete'].append(name);save_json(out/'state.json',state)
        print('DONE '+name+' '+json.dumps(series[-1]),flush=True)
    save_json(out/'results.json',all_results)
    state['status']='complete';state['elapsed_seconds']=time.time()-state['started_epoch'];save_json(out/'state.json',state)
    print('ALL SIX COMPLETE',flush=True)

if __name__=='__main__':main()
