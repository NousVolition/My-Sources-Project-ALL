"""Saved diagnostic gates must reject unsafe, omitted and stale CFL records."""
import pytest
from analyze_coupled import saved_cfl_check
from verify_coupled import verify_saved_cfl


def rows(value):return [{'rows':[{'baseline':{'cfl':.1},'disturbed':{'cfl':value}}]}]


@pytest.mark.parametrize('value',[.5,.6,float('nan'),float('inf'),-.1])
def test_saved_cfl_rejects_invalid_values(value):
    assert saved_cfl_check(rows(value),.5)[1] is False


def test_cfl_passes_below_limit_and_uses_largest_arm():
    assert saved_cfl_check(rows(.3),.5)==(.3,True)


def test_empty_cfl_record_fails():
    assert saved_cfl_check([],.5)[1] is False


@pytest.mark.parametrize('analysis',[
    {'gates':{},'max_saved_advective_cfl':.3},
    {'gates':{'cfl':True}},
    {'gates':{'cfl':True},'max_saved_advective_cfl':.2},
    {'gates':{'cfl':True},'max_saved_advective_cfl':float('nan')},
])
def test_verifier_rejects_missing_or_stale_cfl_gate(analysis):
    with pytest.raises(ValueError):verify_saved_cfl(rows(.3),analysis,{'gates':{'max_advective_cfl':.5}})


def test_verifier_recomputes_saved_values_instead_of_trusting_pass_flag():
    with pytest.raises(ValueError):
        verify_saved_cfl(rows(.6),{'gates':{'cfl':True},'max_saved_advective_cfl':.6},
                         {'gates':{'max_advective_cfl':.5}})
