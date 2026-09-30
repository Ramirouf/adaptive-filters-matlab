"""Generate data/h.csv, the room impulse response of Sec. 4 and Fig. 2 of the paper.

The paper computes h with the software of its ref. [32], rir-generator, and states the
geometry but not the speed of sound or whether the generator's high-pass filter is on.
c = 340 m/s with the high-pass filter on is the setting that reproduces Fig. 2 (peak
0.0501 at tap 33, dip at tap 32, undershoot at tap 35). h is NOT normalized: the paper's
Fig. 2 is on a scale of 0.05, and the knobs of its Tables 1-2 only hold at that scale.

Run from the repo root (needs numpy and rir-generator):
    python data/make_h.py
"""

from pathlib import Path

import numpy as np
import rir_generator

h = rir_generator.generate(
    c=340.0,                      # m/s, identified from Fig. 2
    fs=8000.0,                    # Hz
    r=[1.0, 1.5, 1.0],            # receiver position, m
    s=[1.0, 2.5, 2.0],            # source position, m
    L=[5.0, 10.0, 6.0],           # room dimensions, m
    reverberation_time=0.2,       # RT60, s
    nsample=128,                  # M
    hp_filter=True,               # identified from Fig. 2
)[:, 0]

out = Path(__file__).resolve().parent / "h.csv"
np.savetxt(out, h, fmt="%.17g",
           header="h, M = 128 taps, rir-generator c=340 hp_filter=True (data/make_h.py)",
           comments="% ")
print(f"wrote {out}: ||h|| = {np.linalg.norm(h):.5f}, "
      f"peak {h.max():.4f} at tap {h.argmax()} (0-based)")
