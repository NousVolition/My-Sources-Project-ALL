"""Beads and molecules are never treated as independent statistical replicas."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import t as student_t

def dump(path,value):
    path.write_text(json.dumps(value,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))

def mic(d,box):return d-np.rint(d/box)*box

def snapshot(x,box,edges,bond_edges):
    n=len(x);b=mic(x[:,1:]-x[:,0,None],box);length=np.linalg.norm(b,axis=-1)
    angle=np.rad2deg(np.arccos(np.clip(np.sum(b[:,0]*b[:,1],axis=-1)/np.prod(length,axis=-1),-1,1)))
    delta=mic(x[None,:,0]-x[:,None,0],box);dist=np.linalg.norm(delta,axis=-1);np.fill_diagonal(dist,np.inf)
    unit=np.divide(delta,dist[:,:,None],out=np.zeros_like(delta),where=np.isfinite(dist[:,:,None]))
    nearest=np.argsort(dist,axis=1)[:,:4];neighbors=unit[np.arange(n)[:,None],nearest]
    tetra=np.ones(n)
    for j in range(4):
        for k in range(j):tetra-=3/8*(np.sum(neighbors[:,j]*neighbors[:,k],axis=1)+1/3)**2
    hb=0
    for j in range(2):
        cosine=np.sum(unit*(b[:,j]/length[:,j,None])[:,None,:],axis=-1)
        hb+=np.sum((dist<.35)&(cosine>np.cos(np.pi/6)))
    rdf=np.histogram(dist.ravel(),edges)[0]
    bondhist=np.histogram(length.ravel()*10,bond_edges)[0]
    values=np.array([length.mean()*10,angle.mean(),2*hb/n,tetra.mean(),np.mean(np.sum(dist<.35,axis=1)),np.mean(np.min(dist,axis=1))*10])
    return values,rdf,bondhist

NAMES=['bond_mean_A','angle_mean_deg','hydrogen_bond_degree','tetrahedral_order','coordination_035nm','nearest_OO_A']

def analyze_run(folder):
    cfg=json.loads((folder/'config.json').read_text());box=cfg['side_nm'];n=cfg['n'];beads=cfg['beads']
    data=np.load(folder/'trajectory.npz');x=data['positions_nm'].astype(float)
    edges=np.linspace(0,box/2,101);bond_edges=np.linspace(.7,1.35,326)
    stats=[];rdf=np.zeros(100);bh=np.zeros(len(bond_edges)-1)
    for frame in x:
        vals=[]
        for bead in frame:
            v,r,b=snapshot(bead,box,edges,bond_edges);vals.append(v);rdf+=r;bh+=b
        stats.append(np.mean(vals,axis=0))
    stats=np.asarray(stats);means=stats.mean(axis=0)
    summary=dict(tag=folder.name,**cfg,**dict(zip(NAMES,means)),potential_mean_kJmol=float(np.mean(data['potential_kJmol_per_molecule'])))
    for j,key in enumerate(NAMES):summary[key+'_second_minus_first_half']=float(np.mean(stats[len(stats)//2:,j])-np.mean(stats[:len(stats)//2,j]))
    # An RDF peak is not the same estimator as a mean bond length. Keep both.
    br=(bond_edges[1:]+bond_edges[:-1])/2
    intrag=gaussian_filter1d(bh,2)/br**2
    histogram_mean=float(np.sum(bh*br)/np.sum(bh))
    summary['bond_distribution_sd_A']=float(np.sqrt(np.sum(bh*(br-histogram_mean)**2)/np.sum(bh)))
    summary['intramolecular_RDF_peak_A']=float(br[np.argmax(intrag)])
    summary['experimental_peak_A']=.990 if cfg['isotope']=='H2O' else .985
    summary['experimental_peak_uncertainty_A']=.005
    summary['peak_error_A']=summary['intramolecular_RDF_peak_A']-summary['experimental_peak_A']
    volumes=4*np.pi/3*(edges[1:]**3-edges[:-1]**3)
    rdf=rdf/(len(x)*beads*n*(n/box**3)*volumes)
    raw=dict(time_ps=data['times_ps'],structural_metrics=stats,metric_names=np.array(NAMES),rdf_r_nm=(edges[1:]+edges[:-1])/2,rdf_gOO=rdf,bond_r_A=br,bond_histogram=bh,bond_radial_profile=intrag)
    if (folder/'mobility.npz').exists():
        mob=np.load(folder/'mobility.npz');pos=mob['centroid_nm'].astype(float)[:,:,0];times=mob['times_ps'];dts=cfg['sample_ps']
        # Remove COM translation and unwrap each molecule through periodic boundaries.
        increments=mic(np.diff(pos,axis=0),box);unwrapped=np.concatenate([pos[:1],pos[:1]+np.cumsum(increments,axis=0)])
        unwrapped-=np.mean(unwrapped,axis=1)[:,None,:]
        maxlag=min(round(4/dts),len(times)-2);lags=np.arange(1,maxlag+1)
        msd=np.array([np.mean(np.sum((unwrapped[k:]-unwrapped[:-k])**2,axis=-1)) for k in lags]);lagtime=lags*dts
        Ds=[]
        for lo,hi in [(1.,2.),(2.,4.)]:
            select=(lagtime>=lo-1e-8)&(lagtime<=hi+1e-8)
            Ds.append(float(np.polyfit(lagtime[select],msd[select],1)[0]/6*1000)) # nm2/ps -> 10^-9 m2/s
        energy=mob['ring_energy_kJmol']/beads
        drift=float((energy[-1]-energy[0])/(n*.008314462618*cfg['T_K']))
        energy_range=float(np.ptp(energy)/(n*.008314462618*cfg['T_K']))
        energy_slope=float(np.polyfit(times,energy-energy[0],1)[0]/(n*.008314462618*cfg['T_K']))
        nearest=[]
        for snap in pos:
            d=np.linalg.norm(mic(snap[:,None]-snap[None,:],box),axis=-1);np.fill_diagonal(d,np.inf);nearest.append(np.argsort(d,axis=1)[:,:4])
        nearest=np.array(nearest);lag=round(1/dts)
        retained=np.mean([np.mean([len(set(a)&set(b))/4 for a,b in zip(nearest[t],nearest[t+lag])]) for t in range(len(nearest)-lag)])
        disagreement=abs(Ds[1]-Ds[0])/max(abs(Ds[1]),1e-12)
        summary.update(short_D_1_2ps_1e9_m2s=Ds[0],short_D_2_4ps_1e9_m2s=Ds[1],diffusion_window_relative_disagreement=disagreement,diffusion_window_pass=disagreement<.2,NVE_energy_drift_NkBT=drift,NVE_energy_range_NkBT=energy_range,NVE_energy_linear_slope_NkBT_per_ps=energy_slope,NVE_energy_drift_pass=abs(drift)<.005,neighbor_replaced_fraction_1ps=1-retained)
        raw.update(msd_lag_ps=lagtime,msd_nm2=msd)
    np.savez_compressed(folder/'analysis.npz',**raw);dump(folder/'summary.json',summary)
    return summary

def interval(values):
    values=np.asarray(values);n=len(values);mean=float(values.mean());se=float(values.std(ddof=1)/np.sqrt(n)) if n>1 else None
    half=float(student_t.ppf(.975,n-1)*se) if n>1 else None
    return dict(n=n,mean=mean,CI95=[mean-half,mean+half] if n>1 else None,independent_run_values=values.tolist())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',required=True,type=Path);a=ap.parse_args()
    results=[]
    for folder in sorted(a.data.glob('n*')):
        if not (folder/'timing.json').exists():continue
        r=analyze_run(folder);results.append(r);print(folder.name,round(r['bond_mean_A'],5),round(r['hydrogen_bond_degree'],4),flush=True)
    main=[r for r in results if r['n']==64 and r['beads'] in [1,32] and r['T_K']==300 and r['dt_ps']==.0005 and r['mobility_NVE_ps']==5.]
    keys=NAMES+['intramolecular_RDF_peak_A','bond_distribution_sd_A','short_D_2_4ps_1e9_m2s','neighbor_replaced_fraction_1ps']
    groups={}
    for isotope in ['H2O','D2O']:
        for beads in [1,32]:
            group=[r for r in main if r['isotope']==isotope and r['beads']==beads]
            if group:groups[f'{isotope}-P{beads}']={key:interval([r[key] for r in group]) for key in keys if all(key in r for r in group)}
    contrasts={}
    lookup={(r['isotope'],r['beads'],r['seed']):r for r in main}
    complete=[s for s in [1,2,3] if all((i,b,s) in lookup for i in ['H2O','D2O'] for b in [1,32])]
    if complete:
        for key in keys:
            contrasts[key]={}
            for name,terms in [('classical_D_minus_H',[('D2O',1,1),('H2O',1,-1)]),('quantum_D_minus_H',[('D2O',32,1),('H2O',32,-1)]),('NQE_H',[('H2O',32,1),('H2O',1,-1)]),('NQE_D',[('D2O',32,1),('D2O',1,-1)]),('isotope_by_quantum_interaction',[('D2O',32,1),('H2O',32,-1),('D2O',1,-1),('H2O',1,1)])]:
                if all(key in lookup[i,b,s] for i,b,_ in terms for s in complete):
                    contrasts[key][name]=interval([sum(sign*lookup[i,b,s][key] for i,b,sign in terms) for s in complete])
    controls=[]
    ref=lookup.get(('H2O',32,1))
    if ref:
        for r in results:
            if r['T_K']==300 and r['isotope']=='H2O' and r['seed']==1 and (r['n']!=64 or r['beads'] not in [1,32] or r['dt_ps']!=.0005):
                dr=r['bond_mean_A']-ref['bond_mean_A'];dh=r['hydrogen_bond_degree']-ref['hydrogen_bond_degree']
                controls.append(dict(tag=r['tag'],bond_difference_A=dr,hbond_difference=dh,structural_screen_pass=abs(dr)<.002 and abs(dh)<.15))
    dump(a.data/'summary.json',dict(runs=results,groups=groups,paired_contrasts=contrasts,controls=controls,complete_main_seeds=complete))

if __name__=='__main__':main()
