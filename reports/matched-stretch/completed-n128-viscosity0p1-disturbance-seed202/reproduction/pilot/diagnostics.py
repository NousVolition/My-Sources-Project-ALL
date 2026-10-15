"""Saved-field comparison only; never integrate the reference trajectory."""
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + "." + str(os.getpid()) + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for attempt in range(20):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            if attempt == 19:
                raise
            time.sleep(.1)


def compare(flow, disturbed, reference, norm0, initial_D):
    delta = disturbed - reference
    wh = flow.curl(delta)
    difference_energy = .5 * flow.inner(delta, delta)
    D = np.sqrt(2 * difference_energy) / norm0
    # Independently check the norm with a physical-space inverse FFT.
    du = np.fft.irfftn(delta, s=(flow.n,) * 3, axes=(-3, -2, -1))
    physical_D = np.sqrt(np.sum(du * du) * flow.dx ** 3) / norm0
    error = abs(D - physical_D) / max(D, 1e-30)
    assert error < 1e-10, "Fourier/physical separation norms disagree"
    # Instantaneous energy of the difference: net rate = transfer - viscous loss.
    rp, _, _ = flow.rhs(disturbed)
    rb, _, _ = flow.rhs(reference)
    net = flow.inner(delta, rp - rb)
    loss = flow.nu * flow.inner(wh, wh)
    transfer = net + loss
    # Independent periodic incompressible identity: -integral delta_i delta_j d_j U_i.
    direct_transfer = 0.
    for i in range(3):
        for j in range(3):
            grad = flow.real(1j * flow.k[j] * reference[i])
            direct_transfer -= float(np.sum(du[i] * du[j] * grad) * flow.dx ** 3)
    transfer_error = abs(transfer - direct_transfer) / max(abs(transfer), abs(loss), 1e-20)
    assert transfer_error < 1e-8, "Difference-energy identity does not close"
    return {"D_relative_initial_reference_l2": float(D),
            "amplification": float(D / initial_D),
            "difference_energy": float(difference_energy),
            "difference_energy_net_rate": net,
            "difference_energy_transfer_from_reference": transfer,
            "difference_energy_viscous_loss_rate": loss,
            "instantaneous_log_norm_rate": net / max(2 * difference_energy, 1e-300),
            "physical_D": float(physical_D), "norm_relative_error": float(error),
            "difference_energy_identity_relative_error": float(transfer_error),
            "difference_mean_max": float(np.max(abs(delta[:, 0, 0, 0].real / flow.n ** 3)))}
