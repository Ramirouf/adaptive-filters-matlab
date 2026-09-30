# Decisions

## 2026-09-29: Reproduce the paper's own experiment, not a synthetic benchmark

- **Decided:** reproduce Sec. 4 of the paper (M = 128 room impulse response, AR(1) input, generalized Gaussian noise with beta* = 0.2, SNR 5 dB, N = 100) and its Figs. 3 and 4, with the knobs of Tables 1 and 2 used verbatim.
- **Why:** every result can be checked against a published figure, and the paper's knobs need no recalibration.
- **Rejected:** a textbook benchmark (M = 16 random taps, white input, Gaussian / Laplacian / impulsive noise). No published results to check against, every filter would need a new calibration, and it reproduces nothing.

## 2026-09-29: The paper is the ground truth; the Python notebook is the cross-check

- **Decided:** the code is written from the paper's equations, with the equation number on each update. The Python reproduction (`bayesian_adaptive_filters.ipynb`, package `af/`) is used only to cross-validate, on identical signals.
- **Why:** the notebook is itself an attempt to reproduce the paper. Agreement with it shows two codes compute the same thing, not that either matches the paper; the match with the paper is checked separately (Figs. 3-4, closed-form posterior in `verify/check_model.m`).

## 2026-09-29: h from rir-generator, the software the paper cites

- **Decided:** `data/h.csv` comes from `rir-generator` (Python port of ref. [32]) with c = 340 m/s and the high-pass filter on, not normalized (`||h|| = 0.0599`).
- **Why:** the paper names the generator but not these two settings; this pair reproduces Fig. 2 (peak 0.0501 at tap 33). The knobs of Tables 1-2 only hold at the paper's scale of `h`.
- **Rejected:** the notebook's own image-source `h`. It is not what the paper used. The notebook found its own `h` gives slightly better agreement with Tables 1-2 when recalibrating (rms log ratio 0.108 against 0.144), but with the published knobs the rir-generator `h` already lands every curve within about 1 dB of its target (`results/summary.csv`).

## 2026-09-29: Values the paper leaves unstated

| Unstated | Used | Why |
|---|---|---|
| `v_u` | 1, so `var(x) = 5.263` | Unit-variance driving noise; the notebook uses 1 too. With `v_u = 1` the published knobs land every curve near its target. |
| `v_0` | `\|\|h\|\|^2 / M` | The notebook's choice (0 dB initial misalignment). The paper, p. 3, calls it a free parameter. |
| `tau` at `beta = 1` | `sigma / sqrt(2)`, general Eq. (69) | Eq. (80) prints `sigma sqrt(2)`: that is (69) with `kappa(1) = sqrt(2)`, but Eq. (61) gives `1/sqrt(2)`. Eq. (77) has the same error. |
| vKF knob | the sKF `eps` | Not in the tables; the paper omits vKF from its figures because it matched sKF. |
| Run lengths | 12, 1.2, 1.0, 1.2, 2.0 s | Read off the axes of Figs. 3a-4c. |
| Seed | 1 for every panel | Not stated. All curves of a panel share one input and noise realization. |

## 2026-09-29: One recursion file; the KF runs one realization at a time

- **Decided:** `+af/run_filter.m` holds all five variants, as in the Python `run()`. SG, fKF, sKF and vKF process all N realizations as matrix columns; the full-covariance KF loops over realizations.
- **Why:** the paper's point is that the variants are one recursion, and a single file that mirrors the Python makes cross-validation mismatches easy to locate. Measured in Octave 11.3 (M = 128, N = 100): 24 µs per realization-sample for one realization at a time with a 2-D covariance, against 38 µs looping inside a 3-D array and 48 µs vectorized over a 3-D array. Octave has no `pagemtimes`.
- **Rejected:** one file per algorithm (`lms.m`, `nlms.m`, ...), as first planned.

## 2026-09-29: Cross-validation through a trajectory option in the Python code

- **Decided:** added `w_every` to `af.algorithms.run()` in the Python repo (it stores `w_t` every k samples). `verify/export_reference.py` regenerates the signals `run()` drew and writes them with the Python results; `verify/compare_python.m` compares weights and misalignment. 14 configurations, tolerance 1e-9.
- **Why:** comparing weight vectors is a stronger check than comparing averaged curves. The INRS `verify.py` suite still passes after the change.
- **Checked:** six planted bugs (e.g. `e_{t,i}` in the numerator of Eq. (38), `tau` from Eq. (80), a wrong Table 2 entry) each fail exactly the configurations they affect.

## 2026-09-29: Octave first, MATLAB-compatible subset only

- **Decided:** develop in Octave 11.3. No toolbox functions (own gamma sampler instead of `randg`/`gamrnd`), no Octave-only syntax, one function per file, `QUICK = true` mode (N = 10) for MATLAB Online's 15-minute compute limit.
- **Why:** Octave is the only local install; MATLAB Online Basic is the only MATLAB access.
- **Checked (2026-09-29):** in MATLAB R2026b with no toolboxes installed, Code Analyzer reports no issues and the full `run_all` passes every check (4.3 min on an M1, against 17 min in Octave). The only MATLAB-specific change needed: R2025a and later follow a dark system theme when saving figures, so `sim/make_figure.m` sets `theme(fig, 'light')` inside a `try` (Octave has no `theme`).

## 2026-09-29: Repo outside Google Drive

- **Decided:** `~/Documents/adaptive-filters-matlab`.
- **Why:** a git repo inside a synced Drive folder risks conflicted copies of `.git` files.

## 2026-09-29: Check against the published curves, not only against Python

- **Decided:** the markers on the fKF, sKF and KF curves of Figs. 3a and 3b (27 points; SG has none) are read off a 400 dpi render of the paper by `data/digitize_fig3.py` and stored in `data/paper_fig3_markers.csv`. `verify/compare_paper.m` compares our curves with them (tolerance 1 dB); `run_all.m` runs it.
- **Why:** agreement with the Python code shows two implementations compute the same thing; only a comparison with the published figures checks the reproduction itself. Result at N = 100: largest gap 0.52 dB, within 0.2 dB while the curves are still falling.
- **Also answered:** an apparent early-transient difference in the Fig. 3a sKF, seen on a low-resolution render of the page, was a misreading. At 400 dpi the paper's sKF matches ours (-7.06 dB against -7.20 dB at 1 s).

## 2026-09-29: No golden vectors in this repo

- **Decided:** the FIR golden vectors for the Verilog project live in that project's repo, next to the RTL they check.
- **Why:** they have nothing to do with the paper, and the golden model must change in the same commit as the RTL's rounding and saturation rule.
- **Rejected:** a `golden/` folder here, copied across repos.
