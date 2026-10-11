"""Direct Fourier observation with bounded CUDA reductions.
Same retained polynomial and analytic derivatives as reference_solver.sample;
no interpolation, approximation, single precision or fluid RHS changes.
"""
import numpy as np
from reference_solver import Solver

CUDA=r'''
extern "C" __global__ void observe(
 const double2* h,const double* x,const long long* indices,const int* k,
 double* velocity,double* gradient,int modes,int cells,int points,int centers,double scale){
 int branch=blockIdx.y, point=blockIdx.x, tid=threadIdx.x;
 double sums[12]={0};
 const double* pos=x+(branch*points+point)*3;
 for(int m=tid;m<modes;m+=blockDim.x){
   const int* wave=k+m*3;
   double sn,cs; sincos(pos[0]*wave[0]+pos[1]*wave[1]+pos[2]*wave[2],&sn,&cs);
   double weight=(wave[2]==0?1.0:2.0)*scale;
   for(int c=0;c<3;++c){
     double2 a=h[(branch*3+c)*(long long)cells+indices[m]];
     sums[c]+=(a.x*cs-a.y*sn)*weight;
     double derivative=-(a.x*sn+a.y*cs)*weight;
     for(int j=0;j<3;++j)sums[3+c*3+j]+=derivative*wave[j];
   }
 }
 __shared__ double partial[12][128];
 for(int j=0;j<12;++j)partial[j][tid]=sums[j];
 __syncthreads();
 for(int stride=64;stride>0;stride/=2){
   if(tid<stride)for(int j=0;j<12;++j)partial[j][tid]+=partial[j][tid+stride];
   __syncthreads();
 }
 if(tid==0){
   for(int c=0;c<3;++c)velocity[(branch*points+point)*3+c]=partial[c][0];
   if(point<centers)for(int c=0;c<9;++c)gradient[(branch*centers+point)*9+c]=partial[3+c][0];
 }
}
'''

class FastSolver(Solver):
    def __init__(self,n,nu,backend='gpu'):
        super().__init__(n,nu,backend)
        if backend=='gpu':
            self.observation_kernel=self.xp.RawKernel(CUDA,'observe',options=('--std=c++11',))
            self.observation_indexes=self.xp.asarray(np.flatnonzero(self.keep_cpu.ravel()),dtype=np.int64)
            self.observation_modes=self.xp.asarray(self.modes_cpu,dtype=np.int32)
    def sample(self,h,positions,centers):
        if self.backend=='cpu':return super().sample(h,positions,centers)
        xp=self.xp;b=len(h);points=positions.shape[1]
        if points==0:return xp.zeros((b,0,3)),xp.zeros((b,centers,3,3))
        h=xp.ascontiguousarray(h);positions=xp.ascontiguousarray(positions)
        velocity=xp.empty((b,points,3));gradient=xp.empty((b,centers,3,3))
        self.observation_kernel((points,b),(128,),
          (h,positions,self.observation_indexes,self.observation_modes,velocity,gradient,
           np.int32(len(self.modes_cpu)),np.int32(np.prod(self.shape[1:])),np.int32(points),np.int32(centers),np.float64(1/self.n**3)))
        return velocity,gradient
