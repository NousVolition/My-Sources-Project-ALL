"""Phase portraits, basin map, and numerical controls from saved test results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from test_stability import SYSTEMS,competition

HERE=Path(__file__).resolve().parent;D=HERE/'data'


def field(ax,f,limit=2,density=1):
    grid=np.linspace(-limit,limit,101);X,Y=np.meshgrid(grid,grid)
    out=f(np.stack([X,Y],axis=-1))
    ax.streamplot(grid,grid,out[...,0],out[...,1],density=density,color='#9bafb8',linewidth=.65,arrowsize=.8)
    ax.axhline(0,color='#ddd',lw=.6);ax.axvline(0,color='#ddd',lw=.6)
    ax.set(xlim=(-limit,limit),ylim=(-limit,limit),xlabel='x',ylabel='y',aspect='equal')


def main():
    results=json.loads((D/'results.json').read_text());reverse=json.loads((D/'reversibility.json').read_text())
    data=np.load(D/'trajectories.npz');rd=np.load(D/'reversibility.npz')
    for archive in [data,rd]:
        for name in archive.files:assert np.isfinite(archive[name]).all(),name
    assert np.all(np.linalg.norm(data['grid_final'][data['grid_labels']==1]-[3,0],axis=1)<1e-5)
    assert np.all(np.linalg.norm(data['grid_final'][data['grid_labels']==2]-[0,2],axis=1)<1e-5)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    fig,axes=plt.subplots(4,2,figsize=(10,16),layout='constrained')
    fig.suptitle('Eight linear calibration systems',fontsize=19,fontweight='bold')
    short=['Stable node','Unstable spiral','Unstable defective node','Stable node','Center','Saddle','Line of equilibria','Saddle']
    eig=['−1, −2','2 ± i','1, 1 (one eigenvector)','−1, −4','±3i','1, −1','0, −2','1, −2']
    for ax,(name,raw,_,_,vectors),title,values in zip(axes.flat,SYSTEMS,short,eig):
        A=np.array(raw)
        field(ax,lambda z:z@A.T)
        for vector in vectors:
            v=np.array(vector,dtype=float);v/=np.linalg.norm(v)
            ax.plot([-4*v[0],4*v[0]],[-4*v[1],4*v[1]],'--',color='#cc642c',lw=1.2)
        if name=='5.2.7':
            values_Q,basis=np.linalg.eigh(np.array([[17.,5],[5,2.]]))
            angles=np.linspace(0,2*np.pi,501)
            for level in [.15,.6,1.5]:
                ellipse=basis@(np.sqrt(level/values_Q)[:,None]*np.array([np.cos(angles),np.sin(angles)]))
                ax.plot(*ellipse,color='#087e83',lw=1.1)
        if name=='5.2.9':ax.plot([-1.5,1.5],[-2,2],color='#087e83',lw=2)
        ax.plot(0,0,'o',color='#263c48',ms=4)
        ax.set_title(f'{name} | {title}\nλ = {values}',fontsize=11)
    fig.savefig(HERE/'linear_portraits.png',dpi=150);fig.savefig(HERE/'linear_portraits.svg');plt.close(fig)

    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    fig.suptitle('Does a local linear model predict the nonlinear outcome?',fontsize=18,fontweight='bold')
    starts=data['grid_starts'];labels=data['grid_labels']
    ax[0,0].scatter(starts[:,0],starts[:,1],c=labels,cmap=ListedColormap(['#b6dedd','#f3c6a7']),s=18,marker='s',edgecolor='none')
    curve=data['manifold'];ax[0,0].plot(curve[:,0],curve[:,1],color='#222',lw=1.4,label='Nonlinear basin boundary')
    xx=np.linspace(0,3.6,100)
    ax[0,0].plot(xx,1+(xx-1)/np.sqrt(2),'--',color='#a43c28',lw=1.2,label='Saddle tangent prediction')
    ax[0,0].plot([0,3,1],[2,0,1],'ko',ms=4)
    ax[0,0].text(2.65,.4,'Rabbits',color='#07615f');ax[0,0].text(.25,3,'Sheep',color='#864011')
    ax[0,0].set(title='A | The basin boundary curves away from its tangent',xlabel='Initial rabbit population x',ylabel='Initial sheep population y',xlim=(0,3.6),ylim=(0,3.6),aspect='equal')
    ax[0,0].legend(frameon=False,fontsize=8,loc='upper left')
    time=data['near_boundary_time'];tr=data['near_boundary_trajectories']
    for i,color in [(0,'#087e83'),(1,'#cc642c')]:
        ax[0,1].plot(time,tr[:,i,0],color=color,label=f'Start {i+1}: rabbits')
        ax[0,1].plot(time,tr[:,i,1],'--',color=color,label=f'Start {i+1}: sheep')
    ax[0,1].set(title='B | A 0.002 change in initial sheep count flips the outcome',xlabel='Time',ylabel='Scaled population',xlim=(0,40))
    ax[0,1].legend(frameon=False,fontsize=8)
    for a,color in [(-1.,'#087e83'),(0.,'#4763ab'),(1.,'#cc642c')]:
        ax[1,0].plot(data[f'radial_a{a:g}_time'],data[f'radial_a{a:g}_radius'],color=color,label=f'a={a:g}')
    ax[1,0].set(title='C | Same Jacobian, three different nonlinear outcomes',xlabel='Time',ylabel='Distance to origin',ylim=(0,.25))
    ax[1,0].legend(frameon=False,fontsize=8)
    for point,color in [(0.,'#087e83'),(1.,'#cc642c')]:
        rows=[r for r in results['cubic_local_tests'] if r['fixed_point'][0]==point]
        rho=[r['rho'] for r in rows];error=[r['exact_linearization_error'] for r in rows]
        order=np.polyfit(np.log(rho),np.log(error),1)[0]
        ax[1,1].loglog(rho,error,'-o',color=color,label=f'Near ({point:g},0): slope {order:.2f}')
    ax[1,1].set(title='D | Local approximation improves with smaller perturbations',xlabel='Initial component perturbation ρ',ylabel='Exact nonlinear / linear endpoint difference')
    ax[1,1].legend(frameon=False,fontsize=8)
    for a in ax.flat:a.spines[['top','right']].set_visible(False)
    fig.savefig(HERE/'stability_results.png',dpi=170);fig.savefig(HERE/'stability_results.svg');plt.close(fig)

    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    fig.suptitle('Reversibility and numerical false positives',fontsize=18,fontweight='bold')
    for h in [.1,.05,.025]:
        ax[0,0].plot(data[f'center_h{h:g}_time'],data[f'center_h{h:g}_energy'],label=f'RK4 h={h:g}')
    ax[0,0].axhline(1,color='#222',ls=':',label='Exact invariant')
    ax[0,0].set(title='A | RK4 can make a true center appear to decay',xlabel='Time',ylabel='Quadratic invariant / initial value')
    ax[0,0].legend(frameon=False,fontsize=8)
    for gamma,method,color,style in [(0,'RK4','#087e83','-o'),(.2,'RK4','#cc642c','-o'),(0,'Verlet','#4763ab','--s')]:
        rows=[r for r in reverse['rows'] if r['gamma']==gamma and r['method']==method]
        defects=[max(r['round_trip_defect'],1e-16) for r in rows]
        ax[0,1].loglog([r['h'] for r in rows],defects,style,color=color,label=f'{method}, γ={gamma:g}')
    ax[0,1].set(title='B | Damping produces a return error that survives refinement',xlabel='Time step h',ylabel='Forward / reverse return distance')
    ax[0,1].legend(frameon=False,fontsize=8)
    def f661(z):
        x,y=z[...,0],z[...,1]
        return np.stack([y*(1-x*x),1-y*y],axis=-1)
    def f662(z):
        x,y=z[...,0],z[...,1]
        return np.stack([y,x*np.cos(y)],axis=-1)
    field(ax[1,0],f661,2,1.1)
    for z in [-1,1]:
        ax[1,0].axvline(z,color='#c77b4d',ls=':',lw=1)
        ax[1,0].axhline(z,color='#c77b4d',ls=':',lw=1)
    ax[1,0].plot([1,-1,1,-1],[1,1,-1,-1],'ko',ms=4)
    ax[1,0].set_title('C | 6.6.1: stable / unstable nodes and two saddles',fontsize=11)
    field(ax[1,1],f662,3,1.1)
    for z in [-np.pi/2,np.pi/2]:ax[1,1].axhline(z,color='#c77b4d',ls=':',lw=1)
    ax[1,1].plot(0,0,'ko',ms=4)
    ax[1,1].set_title('D | 6.6.2: a saddle and invariant horizontal lines',fontsize=11)
    for a in ax.flat:a.spines[['top','right']].set_visible(False)
    fig.savefig(HERE/'reversibility_results.png',dpi=170);fig.savefig(HERE/'reversibility_results.svg');plt.close(fig)
    (D/'audit.json').write_text(json.dumps({'finite_arrays':True,'basin_labels_verified':True,'figure_count':3,'passed':True},indent=2),encoding='utf-8')
    print('Three figure sets and saved-array audit completed.')


if __name__=='__main__':main()
