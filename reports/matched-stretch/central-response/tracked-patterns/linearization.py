"""Check the local Jacobian of the actual saved-velocity interpolant."""
import json
import numpy as np
from scipy import fft
from analyze import ROOT,STUDY,wrap,trilinear,save
def jacobian(u,p):
    n=u.shape[1];q=((p+3)%6)*n/6;a=np.floor(q).astype(int);b=q-a;out=np.zeros((len(p),3,3))
    for i in (0,1):
        for j in (0,1):
            for k in (0,1):
                shift=np.array([i,j,k]);ind=(a+shift)%n;values=u[:,ind[:,0],ind[:,1],ind[:,2]].T
                weights=np.where(shift,b,1-b)
                for axis in range(3):
                    w=np.prod(weights[:,np.arange(3)!=axis],axis=1)*(2*shift[axis]-1)*n/6
                    out[:,:,axis]+=values*w[:,None]
    return out
def main():
    d=json.loads((ROOT/'measurements.json').read_text());ij=np.array([[e['i'],e['j']] for e in d['edges']]);types=np.array([e['kind'] for e in d['edges']]);out={}
    for rid,run in d['runs'].items():
        rows=[];h0=None
        for frame in run['frames']:
            t=frame['t'];p=np.array(frame['positions']);path=STUDY/'runs'/rid/f'field-t{t:.2f}.npz'
            with np.load(path,allow_pickle=False) as z:h=z['h']
            if h0 is None:h0=h.copy()
            weights=np.ones(h.shape[-1]);weights[1:-1]=2
            flow_change=float(np.sqrt(np.sum(weights*abs(h-h0)**2)/np.sum(weights*abs(h0)**2)))
            n=run['job']['n'];u=fft.irfftn(h,s=(n,)*3,axes=(-3,-2,-1),workers=1)
            dr=wrap(p[ij[:,1]]-p[ij[:,0]]);middle=p[ij[:,0]]+dr/2;G=jacobian(u,middle)
            probes={}
            for factor in (1.,.5,.25):
                measured=trilinear(u,middle+factor*dr/2)-trilinear(u,middle-factor*dr/2)
                predicted=np.einsum('nij,nj->ni',G,dr*factor)
                residual=np.linalg.norm(measured-predicted,axis=1);size=np.linalg.norm(measured,axis=1)
                probes[str(factor)]={kind:{'relative_rms_error':float(np.linalg.norm(residual[types==kind])/max(np.linalg.norm(size[types==kind]),1e-15)),'max_absolute_error':float(residual[types==kind].max())} for kind in ('central','outer','across')}
            dd=np.linalg.norm(wrap(p[:,None,:]-p[None,:,:]),axis=2);np.fill_diagonal(dd,np.inf);closest=np.unravel_index(dd.argmin(),dd.shape)
            rows.append({'t':t,'velocity_field_relative_change_from_start':flow_change,'minimum_label_separation':float(dd[closest]),'closest_labels':[int(v) for v in closest],'probe_comparisons':probes})
        out[rid]=rows;print('Local matrix '+rid+' complete',flush=True)
    save(ROOT/'linearization.json',{'scope':'Derivative of periodic trilinear saved-velocity interpolation; not a fixed-point stability classification','runs':out,'notes':['The interpolation Jacobian is piecewise defined; cell faces may change its value.','Half and quarter separations are synthetic local probes of the saved field, not newly evolved particles.','Error includes variation of the interpolation Jacobian along a segment.','No fixed autonomous two-coordinate model was introduced.']})
if __name__=='__main__':main()
