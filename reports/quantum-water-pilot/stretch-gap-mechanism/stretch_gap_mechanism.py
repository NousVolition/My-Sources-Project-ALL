"""Separate stretch-gap mechanism diagnostic from saved GFN2-xTB matrices.

Run from the repository root:
python reports/quantum-water-pilot/stretch-gap-mechanism/stretch_gap_mechanism.py
The default source is the unchanged parent study directory.
Requires NumPy. Source files are read only; JSON is written to stdout.
No electronic-structure calculation or alternative-method comparison is run.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

# CODATA 2018 SI values; c is the defined speed of light.
BOHR_A = 0.529177210903
CM1_FACTOR = np.sqrt(4.3597447222071e-18 / (BOHR_A * 1e-10)**2
                     / 1.66053906660e-27) / (2 * np.pi * 299792458 * 100)
Q = np.array([[1., 1., 0.], [1., -1., 0.], [0., 0., np.sqrt(2)]]) / np.sqrt(2)
ORDER = ["bend", "symmetric stretch", "antisymmetric stretch"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def internals(x):
    b = x[1:] - x[0]
    r = np.linalg.norm(b, axis=1)
    return np.array([*r, np.arccos(np.clip(np.dot(b[0], b[1]) / np.prod(r), -1, 1))])


def coordinate_jacobian(x):
    b = x[1:] - x[0]
    r = np.linalg.norm(b, axis=1)
    u = b / r[:, None]
    co = np.dot(u[0], u[1])
    si = np.sqrt(1 - co**2)
    B = np.zeros((3, 9))
    B[0, :3], B[0, 3:6] = -u[0], u[0]
    B[1, :3], B[1, 6:] = -u[1], u[1]
    B[2, 3:6] = -(u[1] - co * u[0]) / (r[0] * si)
    B[2, 6:] = -(u[0] - co * u[1]) / (r[1] * si)
    B[2, :3] = -B[2, 3:6] - B[2, 6:]
    return B, u


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=Path(__file__).resolve().parent.parent)
    args = ap.parse_args()
    root = args.source
    original_manifest = json.loads((root / "checksums.json").read_text(encoding="utf-8"))
    gap_manifest = json.loads((root / "stretch-gap-checksums.json").read_text(encoding="utf-8"))
    cfg = json.loads((root / "protocol.json").read_text(encoding="utf-8"))
    saved = json.loads((root / "vibration-final/runs.json").read_text(encoding="utf-8"))
    records = [r for r in saved if r["method"] == "GFN2-xTB"]
    assert len(records) == 10 and len({r["tag"] for r in records}) == 10
    for iso in ["H2O", "D2O"]:
        subset = [r for r in records if r["isotope"] == iso]
        assert len(subset) == 5
        assert {(r["start"], r["hessian_step_bohr"]) for r in subset} == {
            (0, .005), (1, .005), (2, .005), (0, .0025), (0, .00125)}
    input_names = ["protocol.json", "vibration_comparison.csv", "vibration-final/runs.json",
                   "normal_modes.py", "run_vibrations.py"]
    for rec in records:
        input_names += [f"vibration-final/{rec['tag']}/{name}"
                        for name in ["hessian", "xtbopt.xyz", "settings.inp"]]
    input_hashes = {name: digest(root / name) for name in input_names}
    assert all(original_manifest[name] == sha for name, sha in input_hashes.items())
    assert digest(root / "checksums.json") == "565d465fcf67c2ad34c630a213be2746c4d86361a847a44b3104539e892e6ace"
    assert digest(root / "stretch-gap-checksums.json") == "73ecada483d083e9311962969b3f4fc7ccaed222a6ca51dafacd95b400a63084"
    for name in ["vibration_comparison.csv", "stretch-gap-results.json"]:
        assert digest(root / name) == gap_manifest[name]
    gap = json.loads((root / "stretch-gap-results.json").read_text(encoding="utf-8"))
    with (root / "vibration_comparison.csv").open(encoding="utf-8", newline="") as f:
        position_rows = {(r["isotope"], r["mode"]): float(r["predicted_cm1"])
                         for r in csv.DictReader(f)}
    results = []
    primaries = {}
    primary_hessians = {}
    for rec in records:
        folder = root / "vibration-final" / rec["tag"]
        vals = [float(v.replace("D", "E"))
                for line in (folder / "hessian").read_text(encoding="utf-8").splitlines()
                if not line.strip().startswith("$") for v in line.split()]
        hraw = np.array(vals).reshape(9, 9)
        assert np.max(np.abs(hraw - hraw.T)) < 1e-7
        H = (hraw + hraw.T) / 2
        lines = (folder / "xtbopt.xyz").read_text(encoding="utf-8").splitlines()
        assert int(lines[0]) == 3
        assert [line.split()[0] for line in lines[2:5]] == ["O", "H", "H"]
        xyz_A = np.array([[float(v) for v in line.split()[1:4]] for line in lines[2:5]])
        x = xyz_A / BOHR_A
        B, units = coordinate_jacobian(x)
        # Check the analytic coordinate derivative without altering any force matrix.
        B_fd = np.empty_like(B)
        for j in range(9):
            step = np.zeros((3, 3))
            step.flat[j] = 1e-5
            B_fd[:, j] = (internals(x + step) - internals(x - step)) / 2e-5
        jac_error = float(np.max(np.abs(B - B_fd)))
        iso = rec["isotope"]
        mH = cfg["mass_u"]["H1" if iso == "H2O" else "D2"]
        masses = np.repeat([cfg["mass_u"]["O16"], mH, mH], 3)
        inv_mass = 1 / masses
        weighted = H * np.sqrt(np.outer(inv_mass, inv_mass))
        eig, vectors = np.linalg.eigh(weighted)
        assert np.all(eig[-3:] > 0)
        assignment = {"bend": 6}
        changes_by_mode = {}
        for i in [7, 8]:
            displacement = (vectors[:, i] / np.sqrt(masses)).reshape(3, 3)
            dr = np.sum((displacement[1:] - displacement[0]) * units, axis=1)
            name = "symmetric stretch" if np.prod(dr) > 0 else "antisymmetric stretch"
            assignment[name] = i
            changes_by_mode[name] = dr.tolist()
        assert set(assignment) == set(ORDER)
        freq = {name: float(CM1_FACTOR * np.sqrt(eig[i]))
                for name, i in assignment.items()}
        replay_error = max(abs(freq[name] - rec["frequency_cm1"][i])
                           for i, name in enumerate(ORDER))
        P = np.linalg.pinv(B)
        F = P.T @ H @ P
        G = (B * inv_mass) @ B.T
        recon_error = float(np.linalg.norm(H - B.T @ F @ B) / np.linalg.norm(H))
        fg_eig = np.linalg.eigvals(F @ G)
        assert np.max(abs(fg_eig.imag)) < 1e-12 and np.min(fg_eig.real) > 0
        internal_freq = np.sort(np.sqrt(fg_eig.real) * CM1_FACTOR)
        internal_error = float(np.max(abs(internal_freq - sorted(freq.values()))))
        assert jac_error < 1e-8
        assert replay_error < 1e-5 and internal_error < 1e-5 and recon_error < 1e-8
        row = {
            "tag": rec["tag"], "isotope": iso, "start": rec["start"],
            "hessian_step_bohr": rec["hessian_step_bohr"],
            "frequency_cm1": freq,
            "gap_cm1": freq["antisymmetric stretch"] - freq["symmetric stretch"],
            "saved_frequency_reproduction_max_cm1": replay_error,
            "internal_coordinate_frequency_reproduction_max_cm1": internal_error,
            "relative_hessian_reconstruction_error": recon_error,
            "coordinate_jacobian_check_max": jac_error,
            "bond_changes_for_mode_assignment": changes_by_mode,
        }
        results.append(row)
        if rec["start"] == 0 and rec["hessian_step_bohr"] == .00125:
            Fs, Gs = Q @ F @ Q.T, Q @ G @ Q.T
            k, coupling = (F[0, 0] + F[1, 1]) / 2, F[0, 1]
            g, shared_O = (G[0, 0] + G[1, 1]) / 2, G[0, 1]
            la = eig[assignment["antisymmetric stretch"]]
            ls = eig[assignment["symmetric stretch"]]
            square_terms = {
                "mass_geometry": -2 * k * shared_O,
                "bond_force_coupling": -2 * coupling * g,
            }
            square_terms["bend_and_small_numerical_remainder"] = la - ls - sum(square_terms.values())
            scale = CM1_FACTOR**2 / (freq["antisymmetric stretch"] + freq["symmetric stretch"])
            terms = {name: float(value * scale) for name, value in square_terms.items()}
            assert abs(sum(terms.values()) - row["gap_cm1"]) < 1e-8
            # Isolate symmetric-stretch/bend mixing in this coordinate representation.
            block = np.ix_([0, 2], [0, 2])
            symmetric_root = max(np.linalg.eigvals(Fs[block] @ Gs[block]).real)
            bend_term = float((Fs[0, 0] * Gs[0, 0] - symmetric_root) * scale)
            remainder = terms["bend_and_small_numerical_remainder"] - bend_term
            assert abs(remainder) < 1e-5
            csv_error = max(abs(freq[name] - position_rows[(iso, name)]) for name in ORDER)
            saved_gap = next(r for r in gap["results"] if r["isotope"] == iso)
            assert csv_error < 1e-5 and abs(row["gap_cm1"] - float(saved_gap["predicted_gap_cm1"])) < 1e-5
            primaries[iso] = {
                "tag": rec["tag"], "internal_coordinate_order": ["r1_bohr", "r2_bohr", "angle_radians"],
                "F_energy_units": "Eh with coordinates in bohr and radians",
                "F": F.tolist(), "G": G.tolist(),
                "F_symmetric_antisymmetric_angle": Fs.tolist(),
                "G_symmetric_antisymmetric_angle": Gs.tolist(),
                "angle_deg": float(np.degrees(internals(x)[2])),
                "bond_lengths_A": (internals(x)[:2] * BOHR_A).tolist(),
                "gap_components_cm1": terms,
                "bend_only_component_cm1": bend_term,
                "asymmetry_and_reconstruction_remainder_cm1": float(remainder),
                "stretch_force_constant_ratio_symmetric_over_antisymmetric": float(Fs[0, 0] / Fs[1, 1]),
                "stretch_kinetic_factor_ratio_symmetric_over_antisymmetric": float(Gs[0, 0] / Gs[1, 1]),
                "cancelled_fraction_of_mass_geometry_stretch_term": float(-square_terms["bond_force_coupling"] / square_terms["mass_geometry"]),
                "published_csv_reproduction_max_cm1": csv_error,
            }
            primary_hessians[iso] = H
    isotope_hessian_difference = float(np.max(abs(primary_hessians["H2O"] - primary_hessians["D2O"])))
    assert isotope_hessian_difference < 1e-12
    assert all(digest(root / name) == sha for name, sha in input_hashes.items())
    output = {
        "label": "Separate diagnostic addition: source of the small stretch-band gap",
        "source_repository": "NousVolition/My-Sources-Project-ALL",
        "source_commit": "f571de5d8ba4caa1b6393b3954ad7cb04dd4f97a",
        "original_study_commit": "3469351c1ef9e0fdd0e121179a06d5b7fb0db639",
        "scope": "Post-processing saved matrices; no new xTB or Lean run; no parameter changes or alternative-method comparison.",
        "gap_acceptance_criterion": None,
        "gap_classification": "diagnostic only",
        "component_convention": "Algebraic decomposition of squared-frequency difference in r1, r2, angle coordinates, then divide by the full frequency sum. Components are not separate simulations or unique causal attributions.",
        "constants": {"basis": "CODATA 2018", "bohr_angstrom": BOHR_A, "frequency_conversion_cm1": float(CM1_FACTOR)},
        "input_sha256": input_hashes,
        "published_gap_result_sha256": digest(root / "stretch-gap-results.json"),
        "published_gap_manifest_sha256": digest(root / "stretch-gap-checksums.json"),
        "numerical_consistency_checks_passed": True,
        "checks_note": "Checks above concern data reproduction and algebra, not a new scientific gap pass mark.",
        "primary_H_D_hessian_max_difference": isotope_hessian_difference,
        "saved_run_gap_ranges_cm1": {
            iso: {"minimum": min(r["gap_cm1"] for r in results if r["isotope"] == iso),
                  "maximum": max(r["gap_cm1"] for r in results if r["isotope"] == iso)}
            for iso in ["H2O", "D2O"]
        },
        "primary_mechanism": primaries,
        "saved_run_reanalysis": results,
    }
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
