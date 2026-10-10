"""q-TIP4P/F, Habershon et al. JCP 131,024501 (2009), Eqs 1-4/Table I.
Units: nm, ps, dalton, kJ/mol, elementary charge. No isotope-dependent fitting.
"""
import numpy as np
import openmm as mm
from openmm import app,unit
from scipy.spatial.transform import Rotation

PARAMS=dict(epsilon=0.1852*4.184,sigma=0.31589,q=1.1128,gamma=0.73612,
            D=116.09*4.184,alpha=22.87,r0=.09419,ktheta=87.85*4.184,theta0=np.deg2rad(107.4))
MASSES={'H2O':[15.99491461957,1.00782503223,1.00782503223,0.],
        'D2O':[15.99491461957,2.01410177812,2.01410177812,0.]}

def make_system(n,isotope,side=None,cutoff=.9):
    p=PARAMS; system=mm.System()
    bond=mm.CustomBondForce('D*(a^2*(r-r0)^2-a^3*(r-r0)^3+7/12*a^4*(r-r0)^4)')
    for k,v in [('D',p['D']),('a',p['alpha']),('r0',p['r0'])]: bond.addGlobalParameter(k,v)
    angle=mm.HarmonicAngleForce()
    nb=mm.NonbondedForce()
    nb.setNonbondedMethod(nb.PME if side else nb.NoCutoff)
    if side:
        system.setDefaultPeriodicBoxVectors(*[mm.Vec3(*(side*np.eye(3)[k])) for k in range(3)])
        nb.setCutoffDistance(cutoff)
        nb.setEwaldErrorTolerance(1e-5)
        nb.setUseDispersionCorrection(True)
    for i in range(n):
        o=4*i
        for mass in MASSES[isotope]: system.addParticle(mass)
        system.setVirtualSite(o+3,mm.ThreeParticleAverageSite(o,o+1,o+2,p['gamma'],(1-p['gamma'])/2,(1-p['gamma'])/2))
        for h in (o+1,o+2): bond.addBond(o,h,[])
        angle.addAngle(o+1,o,o+2,p['theta0'],p['ktheta'])
        for q,s,e in [(0,p['sigma'],p['epsilon']),(p['q']/2,1,0),(p['q']/2,1,0),(-p['q'],1,0)]:
            nb.addParticle(q,s,e)
        for j in range(4):
            for k in range(j): nb.addException(o+j,o+k,0,1,0)
    system.addForce(bond);system.addForce(angle);system.addForce(nb)
    return system

def initial_positions(n,side,seed):
    model=app.Modeller(app.Topology(),[])
    model.addSolvent(app.ForceField('tip3p.xml'),numAdded=n)
    raw=np.asarray(model.positions.value_in_unit(unit.nanometer)).reshape(n,3,3)
    old=np.asarray(model.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer))[0,0]
    oxy=raw[:,0]*side/old
    rng=np.random.default_rng(seed)
    theta=PARAMS['theta0']/2;r=PARAMS['r0']
    ref=np.array([[0,0,0],[r*np.sin(theta),0,r*np.cos(theta)],[-r*np.sin(theta),0,r*np.cos(theta)]])
    rotations=Rotation.random(n,random_state=rng).as_matrix()
    xyz=np.zeros((n,4,3));xyz[:,:3]=np.einsum('nij,aj->nai',rotations,ref)+oxy[:,None,:]
    xyz[:,3]=PARAMS['gamma']*xyz[:,0]+(1-PARAMS['gamma'])*(xyz[:,1]+xyz[:,2])/2
    return xyz.reshape(-1,3)

def numpy_energy(xyz):
    """Independent nonperiodic reference for two or more molecules."""
    p=PARAMS;x=np.asarray(xyz).reshape(-1,4,3).copy()
    x[:,3]=p['gamma']*x[:,0]+(1-p['gamma'])*(x[:,1]+x[:,2])/2
    b=x[:,1:3]-x[:,0,None];r=np.linalg.norm(b,axis=-1);d=r-p['r0']
    energy=np.sum(p['D']*(p['alpha']**2*d*d-p['alpha']**3*d**3+7/12*p['alpha']**4*d**4))
    theta=np.arccos(np.sum(b[:,0]*b[:,1],axis=-1)/np.prod(r,axis=-1))
    energy+=np.sum(.5*p['ktheta']*(theta-p['theta0'])**2)
    charge=np.array([0,p['q']/2,p['q']/2,-p['q']])
    for i in range(len(x)):
        for j in range(i):
            d=np.linalg.norm(x[i,:,None,:]-x[j,None,:,:],axis=-1)
            energy+=np.sum(138.935457644382*np.outer(charge,charge)/d)
            energy+=4*p['epsilon']*((p['sigma']/d[0,0])**12-(p['sigma']/d[0,0])**6)
    return float(energy)
