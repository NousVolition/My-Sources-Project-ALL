"""Independent physical-space measurements of saved paired fields; no evolution."""
import numpy as np

from diagnostics import ROOT, compare, read, save, sha
from run_pair import Flow, verify_sources


def main():
    p, r = read(ROOT / "protocol.json"), read(ROOT / "result.json")
    verify_sources(p)
    assert r["protocol_sha256"] == sha(ROOT / "protocol.json")
    f = Flow(p["grid"], p["viscosity"], workers=1)
    initial_reference = np.load(ROOT / "initial-reference.npy", allow_pickle=False)
    initial_disturbed = np.load(ROOT / "initial-disturbed.npy", allow_pickle=False)
    assert np.array_equal(initial_disturbed, f.perturb(initial_reference, p["amplitude"], p["seed"]))
    assert abs(np.sqrt(f.inner(initial_reference, initial_reference))/p["reference_initial_l2"]-1) < 1e-12
    checks, max_norm_error, max_metric_error = [], 0., 0.
    for i, row in enumerate(r["series"]):
        fields = {}
        for name in ("reference", "disturbed"):
            path = ROOT / f"fields/{name}/field-{i:03d}.npy"
            assert sha(path) == row["field_hashes"][name]
            h = np.load(path, allow_pickle=False)
            assert np.isfinite(h).all()
            fields[name] = h
            u = np.fft.irfftn(h, s=(p["grid"],)*3, axes=(-3, -2, -1))
            w = np.fft.irfftn(f.curl(h), s=(p["grid"],)*3, axes=(-3, -2, -1))
            observed = row if name == "disturbed" else row["reference_diagnostics"]
            actual = np.array([np.sqrt(np.sum(w*w, axis=0)).max(), .5*np.sum(u*u)*f.dx**3, .5*np.sum(w*w)*f.dx**3])
            recorded = np.array([observed[k] for k in ("Wmax", "energy", "enstrophy")])
            error = float(np.max(abs(actual-recorded)/np.maximum(abs(recorded), 1e-30)))
            assert error < 1e-10
            max_metric_error = max(max_metric_error, error)
            assert np.max(abs(u.mean(axis=(1, 2, 3))-p["mean_velocity"])) < 1e-10
            divh = 1j*sum(k*c for k, c in zip(f.k, h))
            div = np.fft.irfftn(divh, s=(p["grid"],)*3, axes=(-3, -2, -1))
            assert np.max(abs(div))*f.dx/max(observed["peak_speed"], 1e-12) < 1e-10
            assert abs(observed["energy_budget_relative_residual"]) < .005
        measured = compare(f, fields["disturbed"], fields["reference"], p["reference_initial_l2"], p["amplitude"])
        error = abs(measured["physical_D"]-row["D_relative_initial_reference_l2"])/max(measured["physical_D"], 1e-30)
        assert error < 1e-10
        max_norm_error = max(max_norm_error, error)
        for key in ("difference_energy", "difference_energy_net_rate", "difference_energy_viscous_loss_rate", "difference_energy_transfer_from_reference"):
            assert abs(measured[key]-row[key]) < 1e-9*max(abs(row[key]), 1e-12)
        assert abs(row["difference_energy_budget_relative_residual"]) < .005
        checks.append({"t": row["t"], "field_hashes": row["field_hashes"]})
    assert r["status"] != "complete" or (len(checks) == 41 and checks[-1]["t"] == .4)
    save(ROOT / "saved-field-verification.json", {"status": "passed", "through": checks[-1]["t"],
         "field_pairs": len(checks), "physical_norm_max_relative_error": max_norm_error,
         "physical_W_energy_enstrophy_max_relative_error": max_metric_error,
         "initial_field_prescribed": True, "fields": checks, "new_simulations": 0})
    print(f"Verified {len(checks)} paired fields through {checks[-1]['t']:.2f}.")


if __name__ == "__main__":
    main()
