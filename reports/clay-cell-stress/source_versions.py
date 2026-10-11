"""Resolve per-run source bytes for both documented solver storage versions."""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parent
def resolve_source(name,digest):
    for folder in [ROOT,ROOT/'solver-versions/dense-jac']:
        p=folder/name
        if p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==digest:return p
    raise AssertionError('Unresolved source version: '+name+' '+digest)
def same_equations(recorded):
    for name,digest in recorded.items():resolve_source(name,digest)
    return all(hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==recorded[n] for n in ['clay_hug.py','contact_hug.py','cell_hug.py'])
