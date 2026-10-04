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


def test_divergence_at_is_zero_for_uniform_field():
    sim = PurePythonNavierStokes3D(N=8)
    for i in range(8):
        for j in range(8):
            for k in range(8):
                assert sim.divergence_at(i, j, k) == 0.0


def test_divergence_at_wraps_periodically():
    sim = PurePythonNavierStokes3D(N=4)
    sim.u[1][0][0] = 1.0
    sim.u[3][0][0] = -1.0
    assert sim.divergence_at(0, 0, 0) == (1.0 - (-1.0)) / 2.0


def test_project_leaves_uniform_field_unchanged():
    sim = PurePythonNavierStokes3D(N=4)
    sim.project(iterations=10)
    assert sim.u[0][0][0] == 0.05
    assert sim.v[0][0][0] == 0.01
    assert sim.w[0][0][0] == 0.02
    assert sim.max_w() == 0.02


def test_project_returns_none_and_keeps_grid_shape():
    sim = PurePythonNavierStokes3D(N=4)
    result = sim.project(iterations=5)
    assert result is None
    assert len(sim.u) == 4
    assert len(sim.u[0]) == 4
    assert len(sim.u[0][0]) == 4


def test_step_with_pressure_advances_scalar():
    sim = PurePythonNavierStokes3D(N=4)
    sim.step_with_pressure(dt=0.1, P_U=3.0)
    assert sim.mean_S() > 0.0


def test_step_with_pressure_with_zero_source_keeps_uniform_field():
    sim = PurePythonNavierStokes3D(N=4)
    sim.step_with_pressure(dt=0.1, P_U=0.0)
    assert sim.mean_S() == 0.0
    assert round(sim.u[0][0][0], 5) == 0.05
    assert round(sim.v[0][0][0], 5) == 0.01
    assert round(sim.w[0][0][0], 5) == 0.02
