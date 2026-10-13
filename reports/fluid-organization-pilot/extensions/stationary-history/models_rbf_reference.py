"""RBF kernel ridge; identical validation grid, training-only scaling for all families."""
import numpy as np
def distance(a,b):
    return np.maximum((a*a).sum(1)[:,None]+(b*b).sum(1)[None,:]-2*a@b.T,0)
def predict(model,x):
    a=(np.asarray(x)-model['mean'])/model['scale'];b=np.array(model['train'])
    return model['ymean']+np.exp(-distance(a,b)/(2*model['bandwidth']**2))@np.array(model['alpha'])*model['yscale']
def fit(x,y,xv,yv,grid):
    x=np.asarray(x);y=np.asarray(y);mean=x.mean(0);scale=x.std(0);scale[scale<1e-12]=1
    a=(x-mean)/scale;av=(xv-mean)/scale;d=distance(a,a);dv=distance(av,a)
    median=np.sqrt(np.median(d[np.triu_indices(len(d),1)]));median=max(median,1e-8)
    ym=float(y.mean());ys=max(float(y.std()),1e-12);yy=(y-ym)/ys;best=None
    for factor in grid['bandwidth_factors']:
        bw=median*factor;k=np.exp(-d/(2*bw*bw));kv=np.exp(-dv/(2*bw*bw))
        for regularization in grid['ridge']:
            alpha=np.linalg.solve(k+regularization*np.eye(len(k)),yy)
            pred=ym+ys*kv@alpha;rmse=float(np.sqrt(np.mean((pred-yv)**2)))
            if best is None or rmse<best[0]:
                model={'mean':mean.tolist(),'scale':scale.tolist(),'train':a.tolist(),'alpha':alpha.tolist(),'ymean':ym,'yscale':ys,'bandwidth':bw,'factor':factor,'ridge':regularization,'validation_rmse':rmse}
                best=(rmse,model)
    return best[1]
