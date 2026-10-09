"""Read completed controls and immutable saved fields; never alters solver data."""
import hashlib,json,os,time
from pathlib import Path
import numpy as np
from diagnostics import Meter,initial,norm,mirror,transform
from study import ROOT,load_solver,save
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'plot-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(p):return json.loads(p.read_text())
def curve(r,key,t):
    series=r['dense'] if key in ('W','I','energy','enstrophy') else r['rows']
    return np.interp(t,[x['t'] for x in series],[x[key] for x in series])
def compare(a,b):
    t=np.linspace(0,min(a['rows'][-1]['t'],b['rows'][-1]['t']),201)
    out={'through':float(t[-1])}
    for k in ('W','I','E','D','energy','enstrophy'):
        av,bv=curve(a,k,t),curve(b,k,t)
        denom=np.linalg.norm(bv)
        if denom<1e-10:
            denom=(1. if k=='E' else b['rows'][0]['U'])*np.sqrt(len(t))
        out[k]=float(np.linalg.norm(av-bv)/denom)
    return out

def main():
    m=load_solver();runs={}
    for p in (ROOT/'runs').glob('*/result.json'):
        try:r=read(p)
        except json.JSONDecodeError:continue
        if r['status']=='complete':runs[r['job']['id']]=r
    validation={};cache=ROOT/'field-verification.json'
    old=read(cache) if cache.exists() else {}
    for name,r in runs.items():
        resultsha=hashlib.sha256((ROOT/'runs'/name/'result.json').read_bytes()).hexdigest()
        if name in old and old[name]['result_sha256']==resultsha:
            validation[name]=old[name];continue
        errors=[]
        for f in r['fields']:
            path=ROOT/'runs'/name/f['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==f['sha256'],path
        _,phi,ops,_=initial(m,r['job']['n'],r['job']['shape']);meter=Meter(m,ops,phi,r['job']['transform'])
        for idx in (0,-1):
            u=np.load(ROOT/'runs'/name/r['fields'][idx]['file'])['u'];v=meter.full(u)
            errors += [abs(v[k]-r['rows'][idx][k])/max(1.,abs(v[k])) for k in ('E','D','W','energy','enstrophy','divergence_max','high_band_energy_fraction','high_band_enstrophy_fraction')]
        assert max(errors)<1e-11,(name,max(errors))
        validation[name]=dict(result_sha256=resultsha,field_hashes_checked=len(r['fields']),max_endpoint_scaled_error=max(errors))
    save(cache,validation)
    pairs={}
    for name,r in runs.items():
        if r['job']['sign']!=1:continue
        other=name.replace('-s1-','-s-1-')
        if other not in runs:continue
        b=runs[other];assert r['dt']==b['dt'];err=[]
        for f,g in zip(r['fields'],b['fields']):
            u=np.load(ROOT/'runs'/name/f['file'])['u'];v=np.load(ROOT/'runs'/other/g['file'])['u']
            u=transform(u,r['job']['transform'],True);v=transform(v,r['job']['transform'],True)
            err.append(norm(v-mirror(u))/norm(u))
        pairs[name]=dict(max_field_relative_error=max(err),max_E_absolute_error=max(abs(x['E']-y['E']) for x,y in zip(r['rows'],b['rows'])),
                         max_D_sum_scaled=max(abs(x['D']+y['D']) for x,y in zip(r['rows'],b['rows']))/r['rows'][0]['U'])
    comparisons={}
    for a,b in ((33,49),(49,65)):
        for s in (0,1,-1):
            ka=f'n{a}-s{s}-shape0-none-base';kb=f'n{b}-s{s}-shape0-none-base'
            if ka in runs and kb in runs:comparisons[ka+' vs '+kb]=compare(runs[ka],runs[kb])
    for s in (0,1,-1):
        a=f'n65-s{s}-shape0-none-base';b=a.replace('-base','-half')
        if a in runs and b in runs:comparisons[a+' vs '+b]=compare(runs[a],runs[b])
    for kind in ('shift1','shifthalf','rotate'):
        a='n33-s1-shape0-none-base';b=a.replace('none',kind)
        if a in runs and b in runs:comparisons[a+' vs '+b]=compare(runs[a],runs[b])
    transform_checks={}
    for kind in ('shift1','shifthalf','rotate'):
        a='n33-s1-shape0-none-base';b=a.replace('none',kind)
        if a not in runs or b not in runs:continue
        ra,rb=runs[a],runs[b];errors=[];peakerrors=[]
        _,phi,ops,_=initial(m,33);meter=Meter(m,ops,phi)
        for fa,fb in zip(ra['fields'],rb['fields']):
            assert fa['step']==fb['step']
            u=np.load(ROOT/'runs'/a/fa['file'])['u'];v=np.load(ROOT/'runs'/b/fb['file'])['u'];back=transform(v,kind,True)
            errors.append(norm(back-u)/norm(u))
            peakerrors.append(abs(meter.basic(back)['W']/meter.basic(u)['W']-1))
        transform_checks[kind]=dict(max_pulled_back_field_relative_error=max(errors),max_pulled_back_W_relative_error=max(peakerrors),raw_sampled_W_curve_relative_difference=comparisons[a+' vs '+b]['W'])
    one=read(ROOT/'one-step-controls.json')
    audit=read(ROOT/'source-initial-audit.json');breath=read(ROOT/'breathing-audit-summary.json')
    budgets={}
    for name,r in runs.items():
        integral=0.;budget=[]
        for i,row in enumerate(r['rows']):
            if i:
                before=r['rows'][i-1];integral+=(row['t']-before['t'])*.5*(row['enstrophy_rate']+before['enstrophy_rate'])
            budget.append({'t':row['t'],'enstrophy_balance_relative':(row['enstrophy']-r['rows'][0]['enstrophy']-integral)/r['rows'][0]['enstrophy']})
        budgets[name]={'quadrature':'Trapezoid on full diagnostic outputs, approximately every 0.02; not every solver step.','series':budget}
    save(ROOT/'enstrophy-budgets.json',budgets)
    out=dict(completed=len(runs),planned=22,built_epoch=time.time(),one_step_max_by_transform={k:max(x['relative_step_residual'] for x in one if x['transform']==k) for k in ('mirror','shift1','shifthalf','rotate')},
        mirror_pairs=pairs,curve_comparisons=comparisons,transformed_fields=transform_checks,
        source_energy_changes_percent={k:100*(v['energy_ratio_to_no_C']-1) for k,v in audit['cases'].items()},
        runs={name:dict(n=r['job']['n'],dt=r['dt'],steps=r['steps'],end_time=r['rows'][-1]['t'],
            initial_energy=r['rows'][0]['energy'],final_W=r['rows'][-1]['W'],final_E=r['rows'][-1]['E'],final_D=r['rows'][-1]['D'],
            max_tail_energy=max(x['high_band_energy_fraction'] for x in r['rows']),max_tail_enstrophy=max(x['high_band_enstrophy_fraction'] for x in r['rows']),
            max_energy_balance_relative=max(abs(x['viscous_energy_balance_residual']) for x in r['dense']),
            max_divergence=max(x['divergence_max'] for x in r['rows'])) for name,r in runs.items()})
    save(ROOT/'analysis.json',out)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    fig,axes=plt.subplots(2,3,figsize=(13,7.8),layout='constrained')
    panels=[('W','Peak spin W'),('I','Accumulated peak spin I'),('E','Unsigned mirror difference E'),('D','Signed projection D'),('energy','Kinetic energy'),('enstrophy','Enstrophy')]
    for n in (33,49,65):
        name=f'n{n}-s1-shape0-none-base'
        if name not in runs:continue
        r=runs[name];t=np.array([x['t'] for x in r['rows']])
        for ax,(key,title) in zip(axes.flat,panels):ax.plot(t,curve(r,key,t),label=f'{n} cubed; dt={r["dt"]:.5f}');ax.set_title(title);ax.set_xlabel('Model time');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Unforced mirror controls | L=6, viscosity=0.01, Heun\nOnly completed runs; common available interval 0 to 0.4')
    fig.savefig(ROOT/'curves.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for n in (33,49,65):
        name=f'n{n}-s1-shape0-none-base'
        if name not in runs:continue
        s=runs[name]['rows'][-1]['spectra'];x=np.array(s['shell_index'])
        for ax,key in zip(axes,('energy','enstrophy')):
            y=np.array(s[key]);line=ax.semilogy(x[y>1e-18],y[y>1e-18],label=f'{n} cubed')[0]
            _,_,ops,_=initial(m,n)
            actual_cutoff=float(np.max(np.abs(ops[0][ops[4]]))/(2*np.pi/m.L))
            ax.axvline(actual_cutoff,color=line.get_color(),alpha=.6,linestyle=':')
            ax.set_title(key.capitalize()+' spectrum at t=0.4');ax.set_xlabel('Radial shell index');ax.grid(alpha=.2);ax.legend()
    fig.suptitle('All velocity components | dotted lines: actual retained component limit\nCubical retained support extends to its diagonal corner')
    fig.savefig(ROOT/'spectra.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.8),layout='constrained')
    for name in ('small lean','gap hug'):
        r=read(ROOT/'breathing-diagnostics'/f"{name.replace(' ','-')}.json");rows=r['rows'];t=[x['t'] for x in rows]
        for ax,key in zip(axes,('W','high_band_enstrophy_fraction','energy_balance_relative')):ax.plot(t,[x[key] for x in rows],label=name)
    for ax,title in zip(axes,('Peak spin W','Enstrophy fraction near cutoff','Energy balance residual / initial energy')):ax.set_title(title,fontsize=11);ax.set_xlabel('Model time');ax.grid(alpha=.2);ax.legend()
    axes[1].axhline(.01,color='black',linestyle=':',label='1% screen')
    fig.suptitle('Saved breathing fields: 33 cubed | L=6 | viscosity=0.01\nSmall lean dt=0.00526316; gap hug dt=0.00298507; saved every step\nPrescribed kicks included in energy accounting')
    fig.savefig(ROOT/'breathing-checks.png',dpi=150);plt.close(fig)
    print(json.dumps({'completed':len(runs),'pairs':len(pairs),'comparisons':comparisons,'mirror_step_max':out['one_step_max_by_transform']['mirror']}),flush=True)
if __name__=='__main__':main()
