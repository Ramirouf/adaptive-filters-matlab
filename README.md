# Bayesian adaptive filters in MATLAB/Octave

An independent reimplementation, for study purposes, of the numerical examples of

> L. Szczecinski, J. Benesty, E. V. Kuhn, "A unified Bayesian perspective for conventional and
> robust adaptive filters", *Signal Processing* 238 (2026) 110189.

The paper derives LMS, regularized NLMS, the broadband Kalman filter and the Kalman filter, plus
their robust (Laplacian-noise) counterparts, from one Bayesian recursion. This repository
implements that recursion in MATLAB/Octave, reproduces Figs. 3 and 4 of the paper with its
published parameters, and cross-validates the result against an earlier Python reproduction
([bayes-robust-adaptive-filtering](https://github.com/Ramirouf/bayes-robust-adaptive-filtering))
on identical input and noise.

**Status: work in progress.** Results, figures and the full write-up are still to come.

## Run

```
octave --no-gui run_all.m          % full run, N = 100 as in the paper
```

In MATLAB (including MATLAB Online, which limits continuous compute to 15 minutes):

```
QUICK = true; run_all              % N = 10
```

No toolboxes are needed.

## Layout

```
run_all.m                 checks, cross-validation, then Figs. 3-4
+af/run_filter.m          the recursion: SG, fKF, sKF, vKF, KF, for any beta and I
+af/                      noise model, AR(1) input, SNR, the paper's Tables 1-2
sim/                      one panel of Fig. 3 or 4: run, plot, save
verify/check_model.m      checks against the paper and closed-form results
verify/compare_python.m   cross-validation against the Python reproduction
verify/export_reference.py  writes verify/reference/ from the Python code
data/h.csv                the room impulse response of Fig. 2 (data/make_h.py)
DECISIONS.md              what was decided, including values the paper leaves unstated
```
