"""Remeasure saved breathing fields and source starts without evolving them."""
import argparse,ast,hashlib,json
from pathlib import Path
import numpy as np
from diagnostics import Meter,initial,mirror,norm,dot
from study import ROOT,load_solver,save

def circulation(u,m,ops,cx):
    """Exact line integral of the retained Fourier interpolant around a rectangle."""
    n=u.shape[1];h=np.fft.fftn(u,axes=(1,2,3))*ops[4]/n**3
    k=np.fft.fftfreq(n)*n*2*np.pi/m.L
    xa,xb=cx-.3,cx+.3;ya,yb=-.3,.3
    phase=lambda a:np.exp(1j*k*(a+m.L/2))
    def integral(a,b):
        v=np.empty_like(k,dtype=complex);v[k==0]=b-a
        v[k!=0]=(phase(b)[k!=0]-phase(a)[k!=0])/(1j*k[k!=0]);return v
    vx=integral(xa,xb)[:,None,None]*(phase(ya)-phase(yb))[None,:,None]*phase(0)[None,None,:]
    vy=(phase(xb)-phase(xa))[:,None,None]*integral(ya,yb)[None,:,None]*phase(0)[None,None,:]
    return float(np.sum(h[0]*vx+h[1]*vy).real)

def source_audit(m):
    base,phi,ops,c=initial(m,33);base=(base+mirror(base))*.5
    sources={'no_C':base,'mirrored_CD':base+(c+mirror(c))*.5,'C_only':base+c,
             'mirrored_C':base+mirror(c),'uneven_70_30':base+.7*c+.3*mirror(c)}
    meter=Meter(m,ops,phi);out={}
    for name,u in sources.items():
        u=np.asarray(m.project(u,*ops[:-1]));v=meter.full(u)
        v['circulation_rectangles_z0']={str(x):circulation(u,m,ops,x) for x in (-.9,.9)}
        v['energy_ratio_to_no_C']=v['energy']/m.energy(base)
        out[name]=v
    save(ROOT/'source-initial-audit.json',dict(n=33,L=6,nu=.01,formula='Existing fluid_sources.py initialization; no new evolution.',
        circulation='Counterclockwise rectangles of width and height 0.6 centered at (+/-0.9,0,0), exact line integrals of retained Fourier interpolant.',cases=out))
    print('Five source starts audited',flush=True)

def breathing(m,source,trace):
    tree=ast.parse(source.read_text());stop=next(i for i,node in enumerate(tree.body) if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='out' for t in node.targets))
    tree.body=tree.body[:stop];ns={'__file__':str(source),'__name__':'read_initialization'};exec(compile(tree,str(source),'exec'),ns)
    phi=np.asarray(ns['phi']);ops=ns['ops'];breath=np.asarray(ns['breath_u']);meter=Meter(m,ops,phi)
    full=json.loads((trace/'results.json').read_text());out={}
    for name,r in full.items():
        folder=trace/name.replace(' ','-');manifest=json.loads((folder/'field-manifest.json').read_text())
        rows=[];previous=None;lastloc=None;I=0.;zint=0.;inputE=0.;inputZ=0.;zrateint=0.
        for item in manifest:
            path=folder/item['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256']
            u=np.load(path)['u'];step=len(rows);dt=r['dt'];v=meter.full(u,lastloc);lastloc=v['peak_index']
            v.update(step=step,t=step*dt)
            if previous is not None:
                kicked=previous+.15*np.sin(2*np.pi*step*dt/.4)*breath
                kb=meter.basic(kicked);inputE+=kb['energy']-rows[-1]['energy'];inputZ+=kb['enstrophy']-rows[-1]['enstrophy']
                zint+=dt*.5*(kb['enstrophy']+v['enstrophy']);I+=dt*.5*(rows[-1]['W']+v['W'])
                wh=meter.curlh(meter.fft(kicked));w=meter.real(wh)
                rhs=np.asarray(m.rhs(kicked,*ops[:-1]));wrhs=meter.real(meter.curlh(meter.fft(rhs)))
                zrateint+=dt*.5*(m.L**3*dot(w,wrhs)+v['enstrophy_rate'])
            v.update(I=I,cumulative_kick_energy=inputE,cumulative_kick_enstrophy=inputZ,
                     energy_balance_relative=(v['energy']-rows[0]['energy']-inputE+2*m.NU*zint)/rows[0]['energy'] if rows else 0.,
                     enstrophy_balance_relative=(v['enstrophy']-rows[0]['enstrophy']-inputZ-zrateint)/rows[0]['enstrophy'] if rows else 0.)
            rows.append(v);previous=u
        out[name]=dict(n=33,L=6,nu=.01,T=.4,dt=r['dt'],steps=r['steps'],rows=rows)
        save(ROOT/'breathing-diagnostics'/f"{name.replace(' ','-')}.json",out[name])
        print('BREATHING AUDIT',name,len(rows),flush=True)
    save(ROOT/'breathing-audit-summary.json',{name:{'saved_fields':len(r['rows']),'last':{k:r['rows'][-1][k] for k in ('E','D','W','I','energy','enstrophy','energy_balance_relative','enstrophy_balance_relative')},'max_tailZ':max(x['high_band_enstrophy_fraction'] for x in r['rows'])} for name,r in out.items()})
if __name__=='__main__':
    local=ROOT.parent/'breathing-rerun-intake/v2'
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-file',type=Path,default=local/'received/breathing-rerun/breathing_rerun.py' if local.exists() else ROOT.parent/'breathing-rerun/breathing_rerun.py')
    ap.add_argument('--trace-dir',type=Path,default=local/'timeline' if local.exists() else ROOT.parent/'breathing-rerun/timeline/rerun')
    args=ap.parse_args();m=load_solver();source_audit(m);breathing(m,args.source_file,args.trace_dir)
