from regime import regime_at, zone_status


def test_regime_is_nuanced_below_threshold():
    omega, regime, intensity = regime_at(18)
    assert round(omega, 2) == 14.40
    assert regime == "Nuanced (Stealth)"
    assert round(intensity, 2) == 2.88


def test_regime_switches_at_time_19():
    omega, regime, intensity = regime_at(19)
    assert round(omega, 2) == 15.20
    assert regime == "Aggressive (Nuance Lost)"
    assert round(intensity, 2) == 10.94


def test_intensity_jumps_when_formula_changes():
    _, _, before = regime_at(18)
    _, _, after = regime_at(19)
    assert after > before * 3


def test_zone_1_is_stealth_and_zone_2_is_forced():
    _, _, first = zone_status(1)
    _, _, second = zone_status(2)
    assert first == "Nuanced (Stealth)"
    assert second == "Forced Aggression"


def test_zones_crash_from_6_onward():
    for zone in range(6, 11):
        _, intensity, status = zone_status(zone)
        assert status == "CRASHED"
        assert intensity > 100
