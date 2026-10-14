"""Focused predeclared stress and accuracy controls; every attempt is retained."""
from experiment import run,ROOT,Parameters
from dataclasses import replace,asdict
import json,sys
base=Parameters(angles=24,layers=3,alpha=16,rho=160,viscosity=.15,preload=.35)
cases={}
def add(name,p,**kw):cases[name]=dict(p=p,**kw)
for n,l in [(24,3),(48,5),(72,7)]:
    add('mesh'+str(n)+'-two',replace(base,angles=n,layers=l),end=.2,max_seconds=1200)
for alpha,rho,visc,end in [(16,160,.15,1),(32,256,.075,1),(64,300,.0375,1),(128,500,.01875,.5)]:
    for fb in [False,True]:add('stress'+str(alpha)+'-'+('two' if fb else 'one'),replace(base,alpha=alpha,rho=rho,viscosity=visc,feedback=fb),end=end,max_seconds=1200)
for fb in [False,True]:
    add('accurate16-'+('two' if fb else 'one'),replace(base,feedback=fb),end=1,method='Radau',rtol=2e-10,atol=2e-12,max_seconds=1200)
    add('quad16-'+('two' if fb else 'one'),replace(base,quadrature=16,feedback=fb),end=.2,max_seconds=900)
    add('zero-'+('two' if fb else 'one'),replace(base,alpha=0,rho=28,viscosity=1.5,feedback=fb),end=2)
for i,seed in enumerate([(1,1,1),(-8,8,27),(1.000001,1,1)]):
    for fb in [False,True]:add('window'+str(i)+'-'+('two' if fb else 'one'),replace(base,alpha=8,rho=80,viscosity=.75,feedback=fb),end=6,initial=seed,max_seconds=1200)
def main():
    protocol={name:dict(parameters=asdict(c['p']),**{k:v for k,v in c.items() if k!='p'}) for name,c in cases.items()}
    if sys.argv[1:]==['--protocol']:
        (ROOT/'data/protocol.json').write_text(json.dumps(dict(cases=protocol,defaults=dict(rtol=2e-9,atol=2e-11,max_step=.01,record_step=.01,method='BDF'),targets=dict(short_state_scaled=1e-4,energy_balance=1e-3,mesh_gap_relative=.05),scope='Finite runs; spatial claims restricted to explicitly compared times'),indent=2)+'\n');return
    for name in sys.argv[1:]:
        if (ROOT/'data'/(name+'.json')).exists():raise ValueError('Use a new attempt name; no silent cache: '+name)
        run(name,**cases[name])
if __name__=='__main__':main()
