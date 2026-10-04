from navier import PurePythonNavierStokes3D


def test_five_steps_match_script_output():
    sim = PurePythonNavierStokes3D(N=8)
    means = []
    for _ in range(5):
        sim.step(dt=0.1, P_U=3.0)
        means.append(round(sim.mean_S(), 5))
        assert round(sim.max_w(), 5) == 0.02000
    assert means == [0.00059, 0.00117, 0.00176, 0.00234, 0.00293]


def test_scalar_grows_when_source_is_positive():
    sim = PurePythonNavierStokes3D(N=4)
    before = sim.mean_S()
    sim.step(dt=0.1, P_U=3.0)
    assert sim.mean_S() > before


def test_zero_source_does_not_raise_scalar():
    sim = PurePythonNavierStokes3D(N=4)
    sim.step(dt=0.1, P_U=0.0)
    assert sim.mean_S() == 0.0


def test_grid_shape_stays_cubic():
    sim = PurePythonNavierStokes3D(N=4)
    sim.step(dt=0.1, P_U=1.0)
    assert len(sim.S) == 4
    assert len(sim.S[0]) == 4
    assert len(sim.S[0][0]) == 4
