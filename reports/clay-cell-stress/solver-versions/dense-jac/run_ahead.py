"""Optional worker for inexpensive future controls from the frozen protocol.

Run while the main workers are still in the initial extreme/sensitivity cases.
Do not use concurrently after main workers reach the same named configurations.
It changes scheduling only, and named results are reused by the main runner.
"""
import json
from run_experiments import DATA,run,Parameters
def main():
    configs=json.loads((DATA/'protocol.json').read_text())['configurations']
    for c in configs:
        name=c['name']
        if not (('-fine-two' in name) or name in ['steady-mesh24','steady-mesh48','steady-mesh72','old-form']):continue
        c=dict(c);p=Parameters(**c.pop('parameters'));run(p=p,**c)
if __name__=='__main__':main()
