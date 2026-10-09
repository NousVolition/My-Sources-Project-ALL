"""Verify saved fluid states and render measurements from the supplied update."""
from pathlib import Path
import argparse,ast,hashlib,json,os
import numpy as np

def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=here/'rerun')
    ap.add_argument('--source-dir',type=Path,default=here.parent);args=ap.parse_args()
    data=args.data.resolve();source=args.source_dir.resolve()
    os.environ.setdefault('MPLCONFIGDIR',str(data/'plot-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from PIL import Image
    results=json.loads((data/'results.json').read_text())
    p=source/'breathing_rerun.py';tree=ast.parse(p.read_text());stop=next(i for i,n in enumerate(tree.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='out' for t in n.targets));tree.body=tree.body[:stop]
    ns={'__file__':str(p),'__name__':'read_templates'};exec(compile(tree,str(p),'exec'),ns)
    phi=np.asarray(ns['phi']);ref=(-np.arange(33))%33
    phi_norm=np.sqrt(np.mean(np.sum(phi*phi,axis=0)));assert abs(phi_norm-1)<1e-13
    norm=lambda u:float(np.sqrt(np.mean(np.sum(u*u,axis=0))))
    def mirror(u):
        a=np.take(u,ref,axis=1);a[0]*=-1;return a
    rows=[];verified=0;worst=0.;mirror_error=0.
    for name,r in results.items():
        folder=data/name.replace(' ','-');manifest=json.loads((folder/'field-manifest.json').read_text())
        views=np.load(folder/'slices.npz')['views'];errors=[]
        for i,record in enumerate(manifest):
            path=folder/record['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
            u=np.load(path)['u'];assert np.isfinite(u).all();a=u-mirror(u);D=float(np.mean(np.sum(a*phi,axis=0)))
            measured={'U':norm(u),'A':norm(a),'D':D,'B':norm(a-D*phi)};measured['E']=measured['A']/measured['U']
            error=max(abs(measured[k]-r['series'][i][k]) for k in measured);assert error<1e-12;errors.append(error)
            expected=np.stack([u[0,:,:,16],u[1,:,:,16],np.sqrt(np.sum(u*u,axis=0))[:,:,16],np.sqrt(np.sum(a*a,axis=0))[:,:,16]])
            assert np.max(abs(views[i]-expected))<1e-12
            if name=='small lean':
                other=np.load(data/'minus-lean'/record['file'])['u'];mirror_error=max(mirror_error,norm(mirror(u)-other)/norm(other))
            verified+=1
        first,last=r['series'][0],r['series'][-1]
        changes=None if first['A']<1e-12 else {k:100*(last[k]/first[k]-1) for k in ('A','U','D','E')}
        log_parts=None if first['A']<1e-12 else dict(actual_difference=float(np.log(last['A']/first['A'])),flow_normalization=float(-np.log(last['U']/first['U'])),relative_difference=float(np.log(last['E']/first['E'])))
        rows.append(dict(case=name,steps=r['steps'],dt=r['dt'],changes_percent=changes,
            start=first,end=last,log_E_accounting=log_parts,
            final_squared_asymmetry_outside_template=None if last['A']<1e-12 else (last['B']/last['A'])**2,
            max_kick_A_change=max(abs(x.get('kick_A_change',0)) for x in r['series']),
            max_kick_D_change=max(abs(x.get('kick_D_change',0)) for x in r['series']),
            max_E=max(x['E'] for x in r['series']),max_divergence=max(x['divergence'] for x in r['series']),
            max_cfl=max(x['speed_cfl'] for x in r['series']),max_saved_field_remeasurement_error=max(errors)))
        worst=max(worst,max(errors))
    analysis=dict(status='all_steps_verified',saved_full_fields=verified,max_remeasurement_absolute_error=worst,
        maximum_full_field_mirror_pair_relative_error=mirror_error,cases=rows,
        equations=['E=A/U','A^2=D^2+B^2','log(E/E0)=log(A/A0)-log(U/U0)'],
        definitions={'U':'sqrt(mean sum u_i^2)','A':'sqrt(mean sum (u_i-Mu_i)^2)','D':'mean sum (u_i-Mu_i)*phi_i','B':'norm(u-Mu-D*phi)','phi':'Fixed unit-norm antisymmetric C template'},
        scope='Same 33-grid driven test and original case-dependent steps. No grid/time refinement or claim of pressure-triggered breathing.',
        animation='Images are saved fluid slices at z=-0.090909 in domain units; arrows show in-plane velocity. Norm diagnostics use the entire 3D field. Frames use actual saved steps, with fixed colour scales within each movie. GIF timing is illustrative playback.')
    (data/'analysis.json').write_text(json.dumps(analysis,indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for ax,(name,r) in zip(axes.flat,results.items()):
        series=r['series'];t=np.array([v['t'] for v in series]);first=series[0]
        if first['A']<1e-12:
            ax.semilogy(t,[v['E'] for v in series],color='#1d4ed8',label='E (roundoff scale)');ax.set_ylabel('E')
        else:
            for key,label,color,style in [('A','Actual difference A','#ca8a04','-'),('U','Flow strength U','#0f766e','-'),('D','Signed projection D','#7c3aed','-'),('E','Relative difference E','#1d4ed8','--')]:
                ax.plot(t,[v[key]/first[key] for v in series],label=label,color=color,ls=style,lw=1.8)
            ax.set_ylabel('Fraction of initial value');ax.axhline(1,color='#999',lw=.6)
        ax.set_title(name);ax.set_xlabel('Model time');ax.grid(alpha=.15);ax.legend(fontsize=7)
    fig.suptitle('Every fluid step measured | breathing drive unchanged',fontsize=14)
    fig.savefig(data/'all-cases.png',dpi=140);plt.close(fig)

    for name in ('small lean','gap hug'):
        r=results[name];folder=data/name.replace(' ','-');saved=np.load(folder/'slices.npz');views=saved['views'];axis=saved['axis']
        chosen=np.flatnonzero(abs(axis)<=1.65);xy=axis[chosen];xx,yy=np.meshgrid(xy,xy,indexing='ij')
        crop=np.take(np.take(views,chosen,axis=2),chosen,axis=3)
        series=r['series'];times=np.array([v['t'] for v in series]);first=series[0]
        fig=plt.figure(figsize=(9,6.4),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1.2,1])
        a1=fig.add_subplot(gs[0,0]);a2=fig.add_subplot(gs[0,1]);a3=fig.add_subplot(gs[1,0]);a4=fig.add_subplot(gs[1,1])
        extent=[xy[0]-.5*ns['dx'],xy[-1]+.5*ns['dx'],xy[0]-.5*ns['dx'],xy[-1]+.5*ns['dx']]
        imgs=[]
        for ax,c,title in ((a1,2,'Flow speed + velocity arrows'),(a2,3,'Actual local mirror difference')):
            im=ax.imshow(crop[0,c].T,origin='lower',extent=extent,vmin=0,vmax=float(crop[:,c].max()),cmap='viridis' if c==2 else 'magma',interpolation='nearest')
            fig.colorbar(im,ax=ax,shrink=.85,label='Velocity units');ax.set(title=title,xlabel='x',ylabel='y');imgs.append(im)
        stride=2;arrow_scale=float(np.sqrt(crop[:,0]**2+crop[:,1]**2).max())*9
        quiver=a1.quiver(xx[::stride,::stride],yy[::stride,::stride],crop[0,0,::stride,::stride],crop[0,1,::stride,::stride],color='white',scale=arrow_scale,width=.004)
        for k,label,c in [('A','Actual difference A','#ca8a04'),('U','Flow strength U','#0f766e'),('E','Relative difference E','#1d4ed8')]:
            a3.plot(times,[v[k]/first[k] for v in series],label=label,color=c,lw=1.8)
        a3.set(title='Why E changes',xlabel='Model time',ylabel='Fraction of initial value');a3.legend(fontsize=7);a3.grid(alpha=.15)
        a4.plot(times,[v['D']/first['D'] for v in series],color='#7c3aed',label='D / initial D')
        a4.plot(times,[v['B']/first['A'] for v in series],color='#be123c',label='Other shapes B / initial A')
        a4.set(title='Change in asymmetry pattern',xlabel='Model time',ylabel='Normalized measurement');a4.legend(fontsize=7);a4.grid(alpha=.15)
        cursors=[a3.axvline(0,color='black',lw=1),a4.axvline(0,color='black',lw=1)]
        title=fig.suptitle('',fontsize=12);frames=[]
        for i in np.unique(np.rint(np.linspace(0,r['steps'],33)).astype(int)):
            imgs[0].set_data(crop[i,2].T);imgs[1].set_data(crop[i,3].T)
            quiver.set_UVC(crop[i,0,::stride,::stride],crop[i,1,::stride,::stride])
            for c in cursors:c.set_xdata([times[i],times[i]])
            title.set_text(f"{name} | time {times[i]:.4f} | actual fluid slice z={r['plane_z']:.3f}\nFull-field measurements below; fixed colour scales")
            fig.canvas.draw();frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3]).convert('P',palette=Image.Palette.ADAPTIVE,colors=96))
            if i==r['steps']:fig.savefig(data/(name.replace(' ','-')+'-last.png'),dpi=100)
        frames[0].save(data/(name.replace(' ','-')+'.gif'),save_all=True,append_images=frames[1:],duration=160,loop=0,optimize=True)
        plt.close(fig)
    print(json.dumps({'verified_fields':verified,'largest_error':worst,'full_field_mirror_error':mirror_error,'cases':[{'case':r['case'],'changes':r['changes_percent'],'outside_fraction':r['final_squared_asymmetry_outside_template']} for r in rows]}),flush=True)

if __name__=='__main__':main()
