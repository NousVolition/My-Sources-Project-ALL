"""Optional worker for distant finer controls while main workers run stress."""
import json
from run_experiments import DATA,run,Parameters
def main():
    configs={c['name']:c for c in json.loads((DATA/'protocol.json').read_text())['configurations']}
    for name in ['steady-mesh72','prior-extreme-new-radau-fine-one','prior-extreme-new-bdf-fine-one','harder64-bdf-fine-one']:
        c=dict(configs[name]);p=Parameters(**c.pop('parameters'));run(p=p,**c)
if __name__=='__main__':main()
