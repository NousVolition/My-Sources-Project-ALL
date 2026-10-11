"""Run independent numerical fields in isolated Python processes (not shared RNGs)."""
import argparse,concurrent.futures,json,os,subprocess,sys
from common import ROOT,read_protocol,save
def one(args):
    env=os.environ.copy();env.update(PYTHONPATH=os.pathsep.join(sys.path),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONIOENCODING='utf-8')
    cmd=[sys.executable,str(ROOT/'simulate.py')]+args
    r=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',env=env)
    log=ROOT/'results/logs'/('_'.join(args).replace('--','').replace('.','p')+'.txt');log.parent.mkdir(parents=True,exist_ok=True);log.write_text(r.stdout+r.stderr,encoding='utf-8')
    print(r.stdout or r.stderr,flush=True)
    if r.returncode:raise RuntimeError(' '.join(args)+'\n'+r.stderr)
    return args
def main(stage):
    p=read_protocol(stage.startswith('screen'));jobs=[]
    if stage in ('test','controls'):
        from forecast import verify_freeze
        verify_freeze()
    if stage=='screen_cubic':
        for regime in p['regimes']:
            for seed in p['marker_revalidation_seeds']:
                base=['--screen','--seed',str(seed),'--regime',regime]
                jobs += [base+['--label','cubic','--dt','.005' if regime=='active' else '.01'],base+['--label','cubic_fine','--n',str(p['pilot_fine_grid'][regime]),'--dt','.005'],base+['--label','cubic_half','--dt','.0025' if regime=='active' else '.005']]
    elif stage=='screen':
        for regime in p['regimes']:
            for seed in p['pilot_seeds']:
                jobs += [['--screen','--seed',str(seed),'--regime',regime,'--label','revised','--dt','.005' if regime=='active' else '.01'],['--screen','--seed',str(seed),'--regime',regime,'--label','fine','--n',str(p['pilot_fine_grid'][regime]),'--dt','.005']]
            jobs += [['--screen','--seed','90000','--regime',regime,'--label','half','--dt','.0025' if regime=='active' else '.005']]
    else:
        for regime in p['regimes']:
            seeds=(p['splits']['train']+p['splits']['validation'] if stage=='calibration' else p['splits'][stage]) if stage!='controls' else p['refinement_seeds']
            for seed in seeds:
                base=['--seed',str(seed),'--regime',regime]
                if stage=='controls':
                    jobs += [base+['--label','half','--dt',str(p['regimes'][regime]['dt']/2)],base+['--label','fine','--n',str(p['fine_grid'][regime]),'--dt',str(min(.005,p['regimes'][regime]['dt']))]]
                else:jobs.append(base)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(one,jobs))
    save('results/batch_'+stage+'.json',{'stage':stage,'jobs':jobs,'completed':len(jobs)})
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['screen','screen_cubic','train','validation','calibration','test','controls']);main(ap.parse_args().stage)
