"""Causal, label-invariant local geometric summaries and disjoint future targets."""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from simulate import L

PRESENT_NAMES = ([f"u_{i}" for i in range(3)]+[f"J_{i}{j}" for i in range(3) for j in range(3)]
                 +[f"strain_eigen_{i}" for i in range(3)] + ["strain_norm","rotation_norm","speed"]
                 +["distance_mean","distance_std","distance_min","distance_max"]
                 +[f"shape_{i}{j}" for i,j in [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]]
                 +["strain_alignment","angle_cos_mean","angle_cos_std"])
HISTORY_NAMES = ["distance_log_change_mean","distance_log_change_std","distance_log_change_max",
                 "distance_log_path","distance_rate_std","distance_midpoint_curvature",
                 "angle_endpoint_change","angle_path","angle_path_std",
                 "neighbor_endpoint_exchange","neighbor_exchange_mean","neighbor_exchange_max",
                 "strain_alignment_change","velocity_alignment_change","strain_shape_change",
                 "past_log_singular_0","past_log_singular_1","past_log_singular_2",
                 "past_log_volume","past_log_condition",
                 "shape_eigen_change_0","shape_eigen_change_1","shape_eigen_change_2"]


def wrap(x): return (x+L/2)%L-L/2


def neighbors(x,k):
    if not 3 <= k < len(x): raise ValueError("Need 3 <= neighbors < markers")
    # Tree avoids O(N^2) history scans. Random positions have no exact distance ties.
    _,idx=cKDTree((x+L/2)%L,boxsize=L).query((x+L/2)%L,k=k+1)
    return idx[:,1:]


def shape(r):
    d=np.linalg.norm(r,axis=-1)
    direction=r/np.maximum(d[...,None],1e-12)
    C=np.einsum("...ki,...kj->...ij",direction,direction)/r.shape[-2]
    return d,direction,C


def present(x,u,J,k):
    idx=neighbors(x,k)
    r=wrap(x[idx]-x[:,None])
    d,v,C=shape(r)
    S=(J+J.transpose(0,2,1))/2
    eig,axes=np.linalg.eigh(S)
    sn=np.linalg.norm(S,axis=(1,2))
    rn=np.linalg.norm((J-J.transpose(0,2,1))/2,axis=(1,2))
    alignment=np.mean(np.einsum("nki,ni->nk",v,axes[:,:,-1])**2,axis=1)
    pair=np.einsum("nki,nli->nkl",v,v)[:,np.triu_indices(k,1)[0],np.triu_indices(k,1)[1]]
    features=np.column_stack([u,J.reshape(len(x),9),eig,sn,rn,np.linalg.norm(u,axis=1),
                              d.mean(1),d.std(1),d.min(1),d.max(1),
                              C[:,0,0],C[:,0,1],C[:,0,2],C[:,1,1],C[:,1,2],C[:,2,2],
                              alignment,pair.mean(1),pair.std(1)])
    assert features.shape[1]==len(PRESENT_NAMES)
    return features,idx


def history(past_positions, current_u, current_J, k, duration, reverse=False):
    """Only accepts the past window; cannot inspect a later record by construction.

    Select the neighbors at the prediction instant, then follow their known past.
    Neighbor exchange is separately measured between consecutive past kNN sets.
    Reverse changes temporal ordering of this same measured past, not future data.
    """
    if len(past_positions)<3 or duration <= 0: raise ValueError("History needs >=3 snapshots")
    idx=neighbors(past_positions[-1],k)
    pos=past_positions[::-1] if reverse else past_positions
    r=wrap(pos[:,idx]-pos[:,:,None,:])
    d,v,C=shape(r)
    logd=np.log(np.maximum(d,1e-12))
    change=(logd[-1]-logd[0])/duration
    steps=np.diff(logd,axis=0)
    angular_step=1-np.sum(v[1:]*v[:-1],axis=-1)
    angular_end=1-np.sum(v[-1]*v[0],axis=-1)
    sets=[neighbors(p,k) for p in pos]
    def exchange(a,b):
        return 1-(a[:,:,None]==b[:,None,:]).any(axis=2).mean(axis=1)
    turnover=np.stack([exchange(a,b) for a,b in zip(sets[:-1],sets[1:])])
    S=(current_J+current_J.transpose(0,2,1))/2
    _,axes=np.linalg.eigh(S)
    direction=axes[:,:,-1]
    align=np.mean(np.einsum("tnki,ni->tnk",v,direction)**2,axis=-1)
    ud=current_u/np.maximum(np.linalg.norm(current_u,axis=1)[:,None],1e-12)
    ua=np.mean(np.einsum("tnki,ni->tnk",v,ud)**2,axis=-1)
    # Least-squares finite-neighborhood deformation, not true infinitesimal FTLE.
    # For row vectors r_final ~= r_initial @ B, F = B.T.
    gram=np.einsum("nki,nkj->nij",r[0],r[0])
    cross=np.einsum("nki,nkj->nij",r[0],r[-1])
    reg=1e-10*np.maximum(np.trace(gram,axis1=1,axis2=2),1e-12)[:,None,None]
    B=np.linalg.solve(gram+reg*np.eye(3),cross)
    sv=np.linalg.svd(B,compute_uv=False)
    lsv=np.log(np.maximum(sv,1e-12))/duration
    ce=np.linalg.eigvalsh(C)
    out=np.column_stack([change.mean(1),change.std(1),change.max(1),
                         abs(steps).sum(0).mean(1)/duration,
                         (steps/(duration/(len(pos)-1))).std(axis=(0,2)),
                         (logd[len(pos)//2]-(logd[0]+logd[-1])/2).mean(1),
                         angular_end.mean(1),angular_step.sum(0).mean(1),angular_step.sum(0).std(1),
                         exchange(sets[0],sets[-1]),turnover.mean(0),turnover.max(0),
                         (align[-1]-align[0])/duration,(ua[-1]-ua[0])/duration,
                         np.einsum("nij,nij->n",C[-1]-C[0],S)/duration,
                         lsv,lsv.sum(1),lsv[:,0]-lsv[:,-1],(ce[-1]-ce[0])/duration])
    assert out.shape[1]==len(HISTORY_NAMES) and np.isfinite(out).all()
    return out


def future_target(tangent_start,tangent_end,horizon):
    # M_end @ inv(M_start), avoiding explicit inversion.
    F=np.linalg.solve(tangent_start.transpose(0,2,1),tangent_end.transpose(0,2,1)).transpose(0,2,1)
    return np.log(np.linalg.svd(F,compute_uv=False)[:,0])/horizon


def index_at(time,t):
    i=int(np.argmin(abs(time-t)))
    if not np.isclose(time[i],t,atol=1e-10): raise ValueError(f"Missing exact time {t}")
    return i


def dataset(paths, protocol, lookback=.2, horizon=.2, k=8, count=512, stride=1):
    result={key:[] for key in ["present","history","reverse","irrelevant","target","group","anchor","particle","position"]}
    for path in paths:
        import json
        metadata=json.loads(Path(path).with_suffix(".json").read_text())
        seed=metadata["seed"]
        with np.load(path,allow_pickle=False) as data:
            time=data["time"]
            positions=data["positions"][:,:count]
            velocity=data["velocity"][:,:count]
            gradient=data["gradient"][:,:count]
            tangent=data["tangent"][:,:count]
        for t in protocol["anchors"]:
            i=index_at(time,t)
            b=index_at(time,t-lookback)
            start=index_at(time,t+protocol["gap"])
            end=index_at(time,t+protocol["gap"]+horizon)
            assert time[b]>=t-lookback-1e-10 and time[i]<time[start]<time[end]
            # Stride coarsens history observation, preserving present and target endpoints.
            window=positions[b:i+1:stride]
            if not np.array_equal(window[-1],positions[i]): raise ValueError("Stride drops present")
            A,_=present(positions[i],velocity[i],gradient[i],k)
            H=history(window,velocity[i],gradient[i],k,lookback)
            R=history(window,velocity[i],gradient[i],k,lookback,reverse=True)
            # An irrelevant marker is chosen as the spatially farthest current marker.
            dd=np.linalg.norm(wrap(positions[i,:,None]-positions[i,None,:]),axis=2)
            donor=dd.argmax(axis=1)
            y=future_target(tangent[start],tangent[end],horizon)
            for key,value in {"present":A,"history":H,"reverse":R,"irrelevant":H[donor],"target":y,
                              "group":np.full(count,seed),"anchor":np.full(count,t),"particle":np.arange(count),
                              "position":positions[i]}.items(): result[key].append(value)
    return {key:np.concatenate(values) for key,values in result.items()}


def shuffle_blocks(values,group,anchor,seed):
    rng=np.random.default_rng(seed)
    out=values.copy()
    for g in np.unique(group):
        for t in np.unique(anchor):
            ix=np.flatnonzero((group==g)&(anchor==t))
            out[ix]=values[rng.permutation(ix)]
    return out
