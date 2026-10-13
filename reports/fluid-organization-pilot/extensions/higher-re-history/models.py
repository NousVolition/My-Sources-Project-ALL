"""Training-only scaling; the same linear/RBF ridge validation grid for every input family."""
import numpy as np

def distance(a,b):return np.maximum((a*a).sum(1)[:,None]+(b*b).sum(1)[None,:]-2*a@b.T,0)
def kernel(a,b,kind,bw):return a@b.T/a.shape[1] if kind=='linear' else np.exp(-distance(a,b)/(2*bw*bw))
def predict(m,x):
 a=(np.asarray(x)-m['mean'])/m['scale'];return m['ymean']+m['yscale']*kernel(a,np.asarray(m['train']),m['kind'],m['bandwidth'])@np.asarray(m['alpha'])
def run_errors(pred,y,rows):return np.sqrt(((np.asarray(pred)-y)**2).reshape(-1,rows).mean(1))
def fit(x,y,xv,yv,grid,rows):
 mean=x.mean(0);scale=x.std(0);scale[scale<1e-12]=1.;a=(x-mean)/scale;v=(xv-mean)/scale
 d=distance(a,a);median=max(float(np.sqrt(np.median(d[np.triu_indices(len(d),1)]))),1e-8)
 ym=float(y.mean());ys=max(float(y.std()),1e-12);yy=(y-ym)/ys;best=None
 for kind,factors in [('linear',[1.]),('rbf',grid['bandwidth_factors'])]:
  for factor in factors:
   bw=median*factor;k=kernel(a,a,kind,bw);kv=kernel(v,a,kind,bw)
   for ridge in grid['ridge']:
    alpha=np.linalg.solve(k+ridge*np.eye(len(k)),yy);pred=ym+ys*kv@alpha;score=float(run_errors(pred,yv,rows).mean())
    if best is None or score<best[0]:
     m=dict(mean=mean.tolist(),scale=scale.tolist(),train=a.tolist(),ymean=ym,yscale=ys,alpha=alpha.tolist(),kind=kind,bandwidth=bw,factor=factor,ridge=ridge,validation_mean_run_rmse=score)
     best=(score,m)
 return best[1]
