"""Run named declared cases without duplicating another worker's assignments."""
import argparse,json
from run_experiments import DATA,run,Parameters
def main():
    ap=argparse.ArgumentParser();ap.add_argument('names',nargs='+');args=ap.parse_args()
    configs={c['name']:c for c in json.loads((DATA/'protocol.json').read_text())['configurations']}
    for name in args.names:
        c=dict(configs[name]);p=Parameters(**c.pop('parameters'));run(p=p,**c)
if __name__=='__main__':main()
