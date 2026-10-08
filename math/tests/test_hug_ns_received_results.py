import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest

FOLDER=Path(__file__).resolve().parents[1]/'imports'/'hug-ns'


def load(name,filename):
    spec=importlib.util.spec_from_file_location(name,FOLDER/filename)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('name', sorted(p.name for p in (FOLDER/'results').glob('*.json')))
def test_received_result_bytes_match_source_archive_manifest(name):
    manifest=json.loads((FOLDER/'provenance.json').read_text(encoding='utf-8'))
    entry=next(r for r in manifest['source_files']+manifest.get('supplemental_files',[])
               if r['published_path']=='results/'+name)
    assert hashlib.sha256((FOLDER/'results'/name).read_bytes()).hexdigest()==entry['sha256']


def test_supplied_result_arithmetic_and_times():
    previous=sys.modules.get('solver')
    sys.modules['solver']=load('hug_ns_results_solver','solver.py')
    try:
        checks=load('hug_ns_received_check','check_received_results.py')
    finally:
        if previous is None:
            del sys.modules['solver']
        else:
            sys.modules['solver']=previous
    actual=checks.check_results()
    saved=json.loads((FOLDER/'received-results-check.json').read_text(encoding='utf-8'))
    assert actual['all_numeric_values_finite']
    assert actual['sha256']==saved['sha256']
    assert actual['series']==saved['series']
