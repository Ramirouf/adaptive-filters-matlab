"""Write verify/reference/, the Python side of the cross-validation.

The reference is af.algorithms.run() from the Python reproduction of the paper
(github.com/Ramirouf/bayes-robust-adaptive-filtering). For each of 14 configurations this
script records the misalignment at every sample and the weights every W_EVERY samples,
and it writes the exact input and noise run() drew, so that verify/compare_python.m can
feed the same signals to af.run_filter and compare.

run() does not return its signals, but it draws them with af.signals.signal_chunks from
np.random.default_rng(seed) and uses that generator for nothing else, so calling
signal_chunks the same way reproduces them exactly. If that ever stops being true, the
comparison fails loudly rather than passing by accident.

Run from this repo's root (needs the Python repo on the path, with run(w_every=...)):
    PYTHONPATH=/path/to/bayes-robust-adaptive-filtering python verify/export_reference.py
"""

from pathlib import Path

import numpy as np

from af import signals
from af.algorithms import run
from af.experiments import TABLE1, TABLE2

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "verify" / "reference"

T, N, SEED, W_EVERY = 4000, 2, 1, 100      # 0.5 s at 8 kHz, 2 realizations
A, V_U, SNR_DB, BETA_STAR = 0.9, 1.0, 5.0, 0.2
VAR_X = V_U / (1.0 - A**2)

# 5 variants x beta in {2, 1} with I = 0, plus the iterated robust filters
CONFIGS = [(v, beta, 0) for beta in (2.0, 1.0) for v in ("SG", "fKF", "sKF", "vKF", "KF")]
CONFIGS += [(v, 1.0, 1) for v in ("fKF", "sKF", "vKF", "KF")]


def knobs(variant, beta, n_iter):
    """Paper Tables 1-2 at -20 dB, as transcribed in the Python repo. vKF takes the sKF value."""
    table = TABLE1[beta] if n_iter == 0 else TABLE2[(-20, n_iter)]
    return table["sKF" if variant == "vKF" else variant]


def save(name, arr, header):
    np.savetxt(OUT / name, np.atleast_2d(arr), fmt="%.17g", delimiter=",",
               header=header, comments="% ")


h = np.loadtxt(ROOT / "data" / "h.csv", comments="%")
v_eta = signals.noise_variance_for_snr(h, SNR_DB, a=A, var_x=VAR_X)

rng = np.random.default_rng(SEED)
blocks = list(signals.signal_chunks(T, N, rng, a=A, var_x=VAR_X, sigma_eta=np.sqrt(v_eta),
                                    beta_true=BETA_STAR))
x = np.concatenate([b[0] for b in blocks], axis=1)       # (N, T)
eta = np.concatenate([b[1] for b in blocks], axis=1)

OUT.mkdir(parents=True, exist_ok=True)
save("setup.csv", [T, N, W_EVERY, A, V_U, SNR_DB, BETA_STAR, v_eta],
     "T, N, w_every, a, v_u, snr_db, beta_star, v_eta")
save("h.csv", h[:, None], "h (column)")
save("x.csv", x.T, "input x, T x N")
save("eta.csv", eta.T, "noise eta, T x N")

for variant, beta, n_iter in CONFIGS:
    r = run(h, variant, beta, T, N=N, seed=SEED, a=A, var_x=VAR_X, snr_db=SNR_DB,
            beta_true=BETA_STAR, n_iter=n_iter, w_every=W_EVERY,
            **knobs(variant, beta, n_iter))
    name = f"{variant}_beta{int(beta)}_I{n_iter}"
    save(f"{name}_mis.csv", r["misalignment"][:, None],
         "misalignment (linear), every sample")
    # (K, N, M) -> K rows of N*M, tap index fastest; matches reshape(w_traj, M*N, K)' in MATLAB
    save(f"{name}_w.csv", r["w_traj"].reshape(r["w_traj"].shape[0], -1),
         f"weights at samples 0, {W_EVERY}, 2*{W_EVERY}, ...; row = snapshot, N*M columns")
    print(f"{name:14s} final misalignment {10 * np.log10(r['misalignment'][-1]):7.2f} dB")

print(f"wrote {OUT}")
