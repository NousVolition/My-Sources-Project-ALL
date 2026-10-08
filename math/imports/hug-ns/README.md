> **Reviewed import.** Start with [REVIEW.md](REVIEW.md). The gate and fluid evolution are retained. The pressure diagnostic was corrected and connected to the runner; [the exact changes](review-changes.patch) are recorded. Three numerical defects remain documented. Longer-run numbers below are supplied claims whose logs were not included. See also [the bound note and its review correction](BOUND.md).

<!-- END IMPORT REVIEW BANNER -->

# Hug start, then the fluid equation

The circle is the thing. The cube is only the grid the calculation sits in. Opposite sides of the cube are the same side, so the flow that leaves one side comes back on the other.

The hug is the smooth start. A flat gate is zero near the sides of the cube and one in the middle. It multiplies the stream function, not the velocity. The velocity is then made divergence-free. After time zero that gate does not act. It is not a wall, a pressure limit, or a memory inside the run.

The step after that is unforced periodic Navier-Stokes:

    d_t u + (u·grad) u = -grad p + nu * laplacian u
    div u = 0

No body force.

## Run

```
pip install -r requirements.txt
python run.py hug --n 32 --time 0.4
python run.py sharper --n 64 --time 2
python run.py tube --n 48 --time 2
python run.py tubes --n 48 --time 2
python run.py compare --n 48 --time 2 --out compare.json
```

`compare` is the pulled tube, and it also prints steepening against smoothing.

## What the runs showed

- Hug ring, nu 0.01: speed falls. On 48³ it was still falling at time 40, peak about 0.06.
- Tighter ring, 64³, nu 0.01: peak 5.47 to 0.86 by time 2. Slope 66 to 3.7.
- Pulled tube, nu 0.002: peak rises, 3.74 to 5.22, then falls to 1.46 by time 4. Energy falls the whole time.
- Same tube, 64³: the climb is larger, 3.60 to 6.16, then below the start by time 1.60.
- Two opposite tubes: peak never beats the start, 2.98 to 2.19 by time 2.
- On the pulled tube, the biggest steepening stays larger than the biggest smoothing at every sample. The speed still fades. Those two maxima may not be in the same spot, and pressure is not in that comparison.

These runs do not prove that every smooth start stays smooth for all time.
