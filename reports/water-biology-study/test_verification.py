"""Regression checks for incomplete reproduction and nested integrity records."""
from pathlib import Path
import json
import numpy as np
import pytest
import verify_package as v


@pytest.fixture
def package(tmp_path,monkeypatch):
    root=tmp_path/'package';actual=tmp_path/'actual'
    (root/'data').mkdir(parents=True);(actual/'data').mkdir(parents=True)
    manifest={}
    for i in range(48):
        name=f'row-{i}.json';manifest['data/'+name]='unused'
        for folder in [root,actual]:(folder/'data'/name).write_text(str(i))
    for name in v.NON_RUN_DATA:manifest['data/'+name]='unused'
    monkeypatch.setattr(v,'ROOT',root)
    return root,actual,manifest


def test_all_documented_reproductions_match(package):
    _,actual,manifest=package
    assert v.compare_reproduction(actual,manifest)==48


@pytest.mark.parametrize('change',['missing','unexpected','empty','missing_directory'])
def test_reject_incomplete_or_unexpected_reproduction(package,change):
    _,actual,manifest=package
    if change=='missing':(actual/'data/row-0.json').unlink()
    elif change=='unexpected':(actual/'data/extra.json').write_text('extra')
    else:
        for p in (actual/'data').iterdir():p.unlink()
        if change=='missing_directory':(actual/'data').rmdir()
    with pytest.raises(ValueError):v.compare_reproduction(actual,manifest)


def test_reject_equal_count_wrong_filename(package):
    _,actual,manifest=package
    (actual/'data/row-0.json').rename(actual/'data/replacement.json')
    with pytest.raises(ValueError,match='missing=.*row-0'):v.compare_reproduction(actual,manifest)


def test_reject_manifest_with_missing_expected_output(package):
    _,actual,manifest=package;manifest.pop('data/row-0.json')
    with pytest.raises(ValueError,match='Expected 48'):v.compare_reproduction(actual,manifest)


@pytest.mark.parametrize('change',['shape','nan','keys'])
def test_reject_invalid_reproduced_array(package,change):
    root,actual,manifest=package
    manifest.pop('data/row-0.json');manifest['data/row-0.npz']='unused'
    for folder in [root,actual]:(folder/'data/row-0.json').unlink()
    np.savez(root/'data/row-0.npz',x=np.ones((2,2)))
    if change=='shape':np.savez(actual/'data/row-0.npz',x=np.ones((1,2)))
    elif change=='nan':np.savez(actual/'data/row-0.npz',x=np.full((2,2),np.nan))
    else:np.savez(actual/'data/row-0.npz',other=np.ones((2,2)))
    with pytest.raises((ValueError,AssertionError)):v.compare_reproduction(actual,manifest)


def test_nested_manifests_are_hashed(tmp_path,monkeypatch):
    monkeypatch.setattr(v,'ROOT',tmp_path)
    for name in ['manifest.json','package-verification.json','coupled-feedback/manifest.json',
                 'coupled-feedback/package-verification.json','coupled-feedback/data.json']:
        p=tmp_path/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('{}')
    assert {p.relative_to(tmp_path).as_posix() for p in v.included()}=={
        'coupled-feedback/manifest.json','coupled-feedback/package-verification.json','coupled-feedback/data.json'}


def test_manifest_coverage_rejects_unlisted_nested_manifest(tmp_path,monkeypatch):
    monkeypatch.setattr(v,'ROOT',tmp_path)
    (tmp_path/'manifest.json').write_text('{}')
    (tmp_path/'coupled-feedback').mkdir()
    (tmp_path/'coupled-feedback/manifest.json').write_text('{}')
    monkeypatch.setattr('sys.argv',['verify_package.py'])
    with pytest.raises(ValueError,match='unlisted=.*coupled-feedback/manifest.json'):v.main()
