"""Optional early execution of the main protocol's forming/release controls."""
from dataclasses import replace
from run_experiments import run,Parameters
def main():
    for plastic in [False,True]:
        p=Parameters(preload=.35,alpha=0,rho=28,viscosity=1.5,plastic=plastic)
        _,state=run('form-'+str(plastic),p,end=6,group='memory',max_step=.04)
        run('release-'+str(plastic),replace(p,preload=0),end=12,state=state,group='memory',max_step=.04)
if __name__=='__main__':main()
