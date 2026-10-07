import itertools

from survival import check_aerospace_state, check_survival, check_system_state, simulate_balance_metaphor


def test_nominal_survives():
    assert simulate_balance_metaphor(1, 1, 1, 1, 1) == 1


def test_one_zero_collapses():
    assert simulate_balance_metaphor(1, 1, 1, 0, 1) == 0


def test_any_zero_collapses():
    for bits in itertools.product([0, 1], repeat=5):
        if 0 in bits:
            assert simulate_balance_metaphor(*bits) == 0


def test_check_survival_strings():
    assert check_survival(True, True, True, True, True).endswith("(1)")
    assert check_survival(True, True, False, True, True).endswith("(0)")


def test_ways_to_failure_names_every_off_bit():
    value, log = check_system_state(0, 1, 0, 1, 0)
    assert value == 0
    assert "Disembodied (E=0)" in log
    assert "Asymmetric (S=0)" in log
    assert "Obsolete (U=0)" in log
    assert "Literalized" not in log


def test_ways_only_all_ones_survive():
    survivors = []
    for combo in itertools.product([0, 1], repeat=5):
        value, log = check_system_state(*combo)
        if value == 1:
            survivors.append(combo)
            assert log == "Nominal: Metaphor survives."
        else:
            assert log.startswith("Failure due to:")
    assert survivors == [(1, 1, 1, 1, 1)]


def test_aerospace_only_all_ones_survive():
    survivors = []
    for combo in itertools.product([0, 1], repeat=5):
        value, log = check_aerospace_state(*combo)
        if value == 1:
            survivors.append(combo)
            assert "Stable Flight" in log
        else:
            assert log.startswith("System Failure due to:")
    assert survivors == [(1, 1, 1, 1, 1)]


def test_aerospace_fault_names():
    value, log = check_aerospace_state(0, 1, 1, 0, 1)
    assert value == 0
    assert "Airframe Compromise (B=0)" in log
    assert "Thrust Stall (T=0)" in log
    assert "Avionics Dark" not in log
