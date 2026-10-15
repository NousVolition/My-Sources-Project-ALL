"""Small reporting-verifier guard tests; no study data or fitting."""
import copy
import numpy as np
from verify_evidence import reconcile_boundary
from evaluate import summarize

def summary(prob):
    return summarize(y,event,pred,prob,group,threshold)
def refuses(fn):
    try:fn()
    except AssertionError:return
    raise AssertionError('Expected refusal')
event=np.array([1,0,1,0]);y=np.array([.2,.3,.6,.8]);pred=y+.01
group=np.array([1,1,2,2]);threshold=.4
reference=np.array([.4,.6,.2,.8]);live=reference.copy();live[0]=np.nextafter(.4,0.)
a=summary(live);b=summary(reference)
checked,record=reconcile_boundary(a,b,live,reference,event,group,threshold)
assert record['boundary_rows']==1 and record['row_indices']==[0]
assert checked['pooled']['precision']==b['pooled']['precision']
assert checked['per_run']['1']['recall']==b['per_run']['1']['recall']
assert checked['pooled']['AUPRC_AP']==a['pooled']['AUPRC_AP']
assert a==summary(live)
assert reconcile_boundary(b,b,reference,reference,event,group,threshold,True)[1] is None
refuses(lambda:reconcile_boundary(a,b,live,reference,event,group,threshold,True))
bad=live.copy();bad[0]=.4-2e-12
refuses(lambda:reconcile_boundary(summary(bad),b,bad,reference,event,group,threshold))
wrong=copy.deepcopy(a);wrong['per_run']['1']['recall']=.25
refuses(lambda:reconcile_boundary(wrong,b,live,reference,event,group,threshold))
wrong_ref=copy.deepcopy(b);wrong_ref['macro']['precision']=.123
refuses(lambda:reconcile_boundary(a,wrong_ref,live,reference,event,group,threshold))
far=reference.copy();far[0]=.4-5e-12;near=far.copy();near[0]-=1e-13
assert reconcile_boundary(summary(near),summary(far),near,far,event,group,threshold)[1] is None
assert np.max(abs(live-reference))==5.551115123125783e-17
print('Reporting roundoff guards: ULP boundary accepted; exact-mode, excessive drift, corrupted host/reference metrics rejected; nonbinary metrics preserved.')
