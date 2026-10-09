"""Completed-only analysis of the driven fluid matrix; running source is read-only."""
import hashlib,json,math,os
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.ndimage import map_coordinates
from run_fluid import ROOT,solver
from source.diagnostics import Meter,initial,mirror,norm,dot
def save(name,v):(ROOT/name).write_bytes((json.dumps(v,indent=2,allow_nan=False)+'\n').encode())
def width(mag,idx,w,dx):
 n=mag.shape[0];d=w/np.linalg.norm(w);e=np.cross(d,np.eye(3)[np.argmin(abs(d))]);e/=np.linalg.norm(e);f=np.cross(d,e)
 radii=np.arange(0,3+dx/8,dx/4);origin=np.array(idx)*dx;values=[]
 for theta in np.linspace(0,np.pi,12,endpoint=False):
  axis=np.cos(theta)*e+np.sin(theta)*f;cross=[]
  for sign in (-1,1):
   line=map_coordinates(mag,((origin+sign*radii[:,None]*axis)%6/dx).T,order=1,mode='grid-wrap');hits=np.flatnonzero(line<mag[idx]/2)
   if len(hits):
    j=hits[0];cross.append(float(radii[j-1]+(mag[idx]/2-line[j-1])*(radii[j]-radii[j-1])/(line[j]-line[j-1])))
  if len(cross)==2:values.append(sum(cross))
 return dict(physical=min(values) if values else None,cells=min(values)/dx if values else None,censored_directions=12-len(values))
def forcing(m,j,b,p,ops):
 x,y,z,*_=m.grids(j['n']);g=np.exp(-8*(x*x+y*y));raw=np.array([g*np.sign(x),np.zeros_like(g),np.zeros_like(g)])
 even=np.array(m.project(raw,*ops[:5]));even=(even+mirror(even))/2;even-=even.mean(axis=(1,2,3),keepdims=True);even/=norm(even)
 F=even if j['parity']=='even-gap' else p*(1 if j['sign']>=0 else -1)
 F-=F.mean(axis=(1,2,3),keepdims=True);F=np.array(m.project(F,*ops[:5]));F/=norm(F)
 return j['A']*norm(b)/.4*F
def compare(a,b):
 t=np.linspace(0,.4,401);out={}
 for k in ('W','I','E','D','energy','enstrophy'):
  av=np.interp(t,[r['t'] for r in a['dense']],[r[k] for r in a['dense']]);bv=np.interp(t,[r['t'] for r in b['dense']],[r[k] for r in b['dense']])
  out[k]=float(np.linalg.norm(av-bv)/max(np.linalg.norm(bv),1e-14))
 return out
def ball(r,A,P,sign,t):
 tb=20*t;period=20*P
 sol=solve_ivp(lambda tt,y:[y[1],r*y[0]-y[0]**3-.4*y[1]+sign*A*np.sin(2*np.pi*tt/period)],(0,float(tb[-1])),[sign*.1,0],t_eval=tb,method='DOP853',rtol=1e-11,atol=1e-13)
 assert sol.success
 return sol.y[0]
def main():
 m=solver();runs={}
 for p in (ROOT/'fluid-runs').glob('*/result.json'):
  try:r=json.loads(p.read_text())
  except json.JSONDecodeError:continue
  if r['status']=='complete':runs[r['job']['id']]=r
 verification={};old=json.loads((ROOT/'fluid-field-verification.json').read_text()) if (ROOT/'fluid-field-verification.json').exists() else {}
 details={};old_details=json.loads((ROOT/'fluid-extended-diagnostics.json').read_text()) if (ROOT/'fluid-extended-diagnostics.json').exists() else {}
 for name,r in runs.items():
  rp=ROOT/'fluid-runs'/name/'result.json';digest=hashlib.sha256(rp.read_bytes()).hexdigest()
  if name in old and old[name]['result_sha256']==digest and name in old_details:
   verification[name]=old[name];details[name]=old_details[name];continue
  b,p,ops,_=initial(m,r['job']['n']);meter=Meter(m,ops,p);F=forcing(m,r['job'],b,p,ops);curlF=meter.real(meter.curlh(meter.fft(F)))
  vals=[];errors=[]
  for i,(field,row) in enumerate(zip(r['fields'],r['rows'])):
   path=ROOT/'fluid-runs'/name/field['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==field['sha256']
   u=np.load(path)['u'];assert np.isfinite(u).all();w=meter.real(meter.curlh(meter.fft(u)));mag=np.sqrt(np.sum(w*w,axis=0));idx=tuple(row['peak_index']);ix=(slice(None),)+idx;direction=w[ix]/mag[idx];sine=np.sin(2*np.pi*row['t']/r['job']['P'])
   force_term=float(np.dot(curlF[ix]*sine,direction))
   vals.append(dict(t=row['t'],minimum_transverse_width=width(mag,idx,w[ix],ops[-1]),
      vorticity_force_along_peak=force_term,vorticity_total_rhs_along_peak=row['vorticity_actual_rhs_along_peak']+force_term,
      velocity_force_vector=(F[ix]*sine).tolist(),instantaneous_force_power=6**3*dot(u,F*sine),
      direct_signed_D_source=dot(F-mirror(F),p)*sine))
   if i in (0,len(r['rows'])-1):
    v=meter.full(u)
    errors += [abs(v[k]-row[k])/max(1,abs(row[k])) for k in ('W','E','D','energy','enstrophy','divergence_max')]
  assert max(errors)<1e-11
  details[name]=vals;verification[name]=dict(result_sha256=digest,fields_checked=len(r['fields']),max_endpoint_scaled_error=max(errors))
  print('VERIFIED',name,flush=True)
 save('fluid-field-verification.json',verification);save('fluid-extended-diagnostics.json',details)
 pairs={}
 for name,r in runs.items():
  if r['job']['sign']!=1:continue
  other=name.replace('-s1-','-s-1-')
  if other not in runs:continue
  b=runs[other];assert r['dt']==b['dt'];errors=[]
  for f,g in zip(r['fields'],b['fields']):
   assert f['step']==g['step']
   u=np.load(ROOT/'fluid-runs'/name/f['file'])['u'];v=np.load(ROOT/'fluid-runs'/other/g['file'])['u'];errors.append(norm(v-mirror(u))/norm(u))
  pairs[name]=dict(max_field_relative_error=max(errors),max_E_difference=max(abs(x['E']-y['E']) for x,y in zip(r['dense'],b['dense'])),max_D_sum=max(abs(x['D']+y['D']) for x,y in zip(r['dense'],b['dense'])))
 comparisons={}
 for parity in ('even-gap','odd-lean'):
  a=f'n49-{parity}-A0.5-P0.2-s1-base'
  for b in (a.replace('n49','n65'),a.replace('-base','-half')):
   if a in runs and b in runs:comparisons[a+' vs '+b]=compare(runs[a],runs[b])
 train='n49-odd-lean-A0.15-P0.2-s1-base';predictions={};fits={}
 if train in runs:
  tr=runs[train];t=np.linspace(0,.4,401);target=np.interp(t,[v['t'] for v in tr['dense']],[v['D']/tr['initial_norm'] for v in tr['dense']]);d0=target[0]
  slope=float(np.dot(t,target-d0)/np.dot(t,t))
  for rpar in (-1,1):
   x=ball(rpar,.15,.2,1,t);basis=x-.1;k=max(0,float(np.dot(basis,target-d0)/np.dot(basis,basis)));fits[str(rpar)]=dict(k=k,linear_baseline_slope=slope,training=train)
   for name,r in runs.items():
    if r['job']['parity']!='odd-lean':continue
    j=r['job'];y=np.interp(t,[v['t'] for v in r['dense']],[v['D']/r['initial_norm'] for v in r['dense']]);yy=y[0]+k*(ball(rpar,j['A'],j['P'],j['sign'],t)-j['sign']*.1)
    constant=np.full_like(y,y[0]);linear=y[0]+j['sign']*slope*t;den=max(np.linalg.norm(y),1e-14)
    predictions[f'r{rpar} '+name]=dict(run=name,r=rpar,training=name==train,mirror_of_training=name==train.replace('-s1-','-s-1-'),relative_curve_error=float(np.linalg.norm(yy-y)/den),constant_error=float(np.linalg.norm(constant-y)/den),linear_error=float(np.linalg.norm(linear-y)/den),passes_10_percent=bool(np.linalg.norm(yy-y)/den<=.1),t=t.tolist(),measured=y.tolist(),predicted=yy.tolist())
 summaries={}
 for name,r in runs.items():
  widths=[v['minimum_transverse_width']['cells'] for v in details[name]]
  summaries[name]=dict(job=r['job'],final={k:r['rows'][-1][k] for k in ('W','I','E','D','energy','enstrophy','gap_center_line_half_peak')},min_width_cells=min(x for x in widths if x is not None),max_tailZ=max(v['high_band_enstrophy_fraction'] for v in r['rows']),max_divergence=max(v['divergence_max'] for v in r['rows']),max_energy_balance=max(abs(v['energy_balance']) for v in r['dense']),max_direct_D_force=max(abs(v['direct_signed_D_source']) for v in details[name]))
 out=dict(completed=len(runs),planned=17,status='complete' if len(runs)==17 else 'partial_completed_only',pairs=pairs,comparisons=comparisons,runs=summaries,fits=fits,predictions=predictions,force_terms='Raw Meter fields are force-free terms. Extended diagnostics add F and curl(F) explicitly.',prediction_protocol_sha256=hashlib.sha256((ROOT/'prediction-protocol.json').read_bytes()).hexdigest())
 save('fluid-analysis.json',out)
 os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(2,3,figsize=(14,7.5),layout='constrained')
 for name,r in runs.items():
  j=r['job']
  if j['n']!=49 or j['sign']!=1 or j['half']:continue
  label=f"{j['parity']}, A={j['A']}, P={j['P']}"
  for ax,key in zip(axes.flat,('D','E','W','I','energy','enstrophy')):
   ax.plot([v['t'] for v in r['dense']],[v[key] for v in r['dense']],label=label);ax.set_title(key);ax.set_xlabel('Fluid time');ax.grid(alpha=.2)
 axes[0,0].legend(fontsize=7);fig.suptitle('Driven fluid comparisons | completed runs only | 49 cubed, viscosity0.01\nD signed; E unsigned; forcing included at both Heun stages')
 fig.savefig(ROOT/'fluid-curves.png',dpi=150);plt.close(fig)
 if predictions:
  fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
  for ax,rpar in zip(axes,(-1,1)):
   for P in (.2,.1):
    key=f'r{rpar} n49-odd-lean-A0.5-P{P}-s1-base'
    if key not in predictions:continue
    v=predictions[key];line=ax.plot(v['t'],v['measured'],label=f'Fluid P={P}')[0];ax.plot(v['t'],v['predicted'],'--',color=line.get_color(),label=f'Ball prediction P={P}')
   ax.set_title(f'Scalar r={rpar}; calibration fixed');ax.set_xlabel('Fluid time');ax.set_ylabel('D / initial velocity norm');ax.grid(alpha=.2);ax.legend(fontsize=8)
  fig.suptitle('Held-out prediction: larger forcing and changed frequency\nOnly low-amplitude P=0.2 run used for calibration')
  fig.savefig(ROOT/'fluid-predictions.png',dpi=150);plt.close(fig)
 print(json.dumps(dict(completed=len(runs),pairs=len(pairs),predictions=len(predictions))),flush=True)
if __name__=='__main__':main()
