"""Fused double-precision loops; no fast-math approximations or parallel reductions."""
import numpy as np
from numba import njit

@njit(cache=True)
def project(h,k,inv,keep):
 out=np.empty_like(h)
 for q in range(h.shape[0]):
  for j in range(h.shape[2]):
   if keep[j]:
    dot=(k[0,j]*h[q,0,j]+k[1,j]*h[q,1,j]+k[2,j]*h[q,2,j])*inv[j]
    for c in range(3):out[q,c,j]=h[q,c,j]-k[c,j]*dot
   else:
    for c in range(3):out[q,c,j]=0
 return out
@njit(cache=True)
def curl(h,k):
 out=np.empty_like(h)
 for q in range(h.shape[0]):
  for j in range(h.shape[2]):
   out[q,0,j]=1j*(k[1,j]*h[q,2,j]-k[2,j]*h[q,1,j])
   out[q,1,j]=1j*(k[2,j]*h[q,0,j]-k[0,j]*h[q,2,j])
   out[q,2,j]=1j*(k[0,j]*h[q,1,j]-k[1,j]*h[q,0,j])
 return out
@njit(cache=True)
def inner(a,b,w,normalization):
 out=np.zeros(a.shape[0])
 for q in range(a.shape[0]):
  for c in range(3):
   for j in range(a.shape[2]):out[q]+=(a[q,c,j].real*b[q,c,j].real+a[q,c,j].imag*b[q,c,j].imag)*w[j]
 return out/normalization
@njit(cache=True)
def cross(u,w):
 out=np.empty_like(u);speed=0.
 for q in range(u.shape[0]):
  for j in range(u.shape[2]):
   a,b,c=u[q,0,j],u[q,1,j],u[q,2,j];x,y,z=w[q,0,j],w[q,1,j],w[q,2,j]
   out[q,0,j]=b*z-c*y;out[q,1,j]=c*x-a*z;out[q,2,j]=a*y-b*x
   speed=max(speed,a*a+b*b+c*c)
 return out,np.sqrt(speed)
@njit(cache=True)
def rhs(h,adv,k,inv,keep,k2,low,w,nu,power,normalization):
 bsize=h.shape[0];norm=np.zeros(bsize);diss=np.zeros(bsize)
 for q in range(bsize):
  for c in range(3):
   for j in range(h.shape[2]):
    e=(h[q,c,j].real**2+h[q,c,j].imag**2)*w[j]
    diss[q]+=nu*k2[j]*e
    if low[j]:norm[q]+=e
 coeff=np.zeros(bsize)
 if power!=0:
  for q in range(bsize):
   if norm[q]<1e-15*normalization:raise ValueError('Empty forcing band')
   coeff[q]=power*normalization/norm[q]
 out=np.empty_like(h);rates=np.empty((bsize,2))
 for q in range(bsize):
  rates[q,0]=power;rates[q,1]=diss[q]/normalization
  for j in range(h.shape[2]):
   dot=(k[0,j]*adv[q,0,j]+k[1,j]*adv[q,1,j]+k[2,j]*adv[q,2,j])*inv[j]
   factor=-nu*k2[j]+coeff[q]*low[j]
   for c in range(3):out[q,c,j]=(adv[q,c,j]-k[c,j]*dot)*keep[j]+factor*h[q,c,j]
 pdq=np.zeros((bsize-1,3))
 for q in range(1,bsize):
  net=0.;dissdiff=0.;forcing=0.
  for c in range(3):
   for j in range(h.shape[2]):
    d=h[q,c,j]-h[0,c,j];r=out[q,c,j]-out[0,c,j]
    net+=(d.real*r.real+d.imag*r.imag)*w[j]
    dissdiff+=nu*k2[j]*(d.real*d.real+d.imag*d.imag)*w[j]
    if low[j]:
     df=coeff[q]*h[q,c,j]-coeff[0]*h[0,c,j];forcing+=(d.real*df.real+d.imag*df.imag)*w[j]
  pdq[q-1,0]=(net+dissdiff-forcing)/normalization;pdq[q-1,1]=dissdiff/normalization;pdq[q-1,2]=forcing/normalization
 return out,rates,pdq
@njit(cache=True)
def axpy(h,k,factor):return h+factor*k
@njit(cache=True)
def rk4(h,a,b,c,d,factor):return h+factor*(a+2*b+2*c+d)
