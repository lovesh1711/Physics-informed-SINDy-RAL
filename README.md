# Physics-Informed Sparse-Identification MPC for Agile Quadrotors: Beyond Lifted-Linear Koopman Predictors

Code to reproduce the simulation results of the paper. It compares a
**physics-informed SINDy** model, identified directly in the state space and
driven by a linear-time-varying MPC, against a **lifted-linear Koopman/EDMD**
predictor on aggressive quadrotor trajectory tracking.

The paper's central point is a *structural* limitation of the lifted-linear
predictor. Quadrotor thrust acts along the body axis, so its effect on the
state depends on attitude, and a **constant** Koopman input matrix can only
hold the *average* thrust direction over the training data. No richer
dictionary removes this defect, and it can be read off a single identified
coefficient, β_z, before flight. The SINDy model carries the bilinear thrust
term explicitly and tracks where the Koopman baselines diverge.

## Requirements

- MATLAB R2021b or newer
- Optimization Toolbox (`quadprog`)
- Statistics and Machine Learning Toolbox (`mvnrnd`, `RandStream`)

No external MPC solver (e.g. CVX) is required.

## Quick start

```matlab
>> main
```

`main.m` runs the full pipeline from scratch (training-data generation, model
identification, and closed-loop MPC) and writes every table's data and every
figure to `results/`. All random seeds are fixed, so every number is
reproduced exactly, except the online compute times, which depend on the
machine.

## What `main.m` produces

| Stage | Script | Paper item | Output |
|-------|--------|------------|--------|
| 1 | `run_diagnostics.m` | Tables I, II | predicted vs fitted β_z on datasets D1–D4; error floor on D4; `results/diagnostics_results.mat` |
| 2 | `run_comparison.m` | Table III | open-loop nRMSE, closed-loop tracking RMSE, online compute time; `results/comparison_results.mat` |
| 3 | `make_paper_figures.m` | Figs. 3, 4 | `results/states.png`, `results/control_inputs.png` |
| 4 | `run_figure8.m` | Sec. IV-D | figure-eight with SINDy, K-narrow (β_z < 0), K-broad (β_z > 0); `results/figure8_results.mat` |
| 5 | `make_figure8_plots.m` | Fig. 5 | `results/figure8_traj.png`, `results/figure8_err.png` |
| 6 | `run_bilinear_ablation.m` | Table IV | EDMD with and without the bilinear thrust atom, trained on D2 |

Stages communicate only through the `.mat` files, so any stage can also be run
on its own once its input `.mat` exists.

## Repository layout

```
main.m                    runs the whole pipeline
run_diagnostics.m         Tables I-II: beta_z diagnostic and error floor
run_comparison.m          Table III: Koopman-EDMD vs SINDy, closed-loop tracking
make_paper_figures.m      Figs. 3-4: states.png and control_inputs.png
run_figure8.m             figure-eight experiment (three controllers)
make_figure8_plots.m      Fig. 5: the two figure-eight panels
run_bilinear_ablation.m   Table IV: bilinear thrust atom ablation
utils/                    quadrotor dynamics, RK4, quaternion/SO(3) helpers,
                          differentially-flat figure-eight reference
koopman/                  EDMD observables (basic and enriched) and
                          reference linearization
sindy/                    candidate library, STLSQ, SINDy RHS and Jacobian
results/                  generated figures and data (created by `main`)
```

## Model and controller, in brief

- **Plant.** Rigid-body quadrotor, state `x = [p; v; q; ω]` (13-D), input
  `u = [f_t; M]`, NED frame, integrated with fixed-step RK4 (1 ms sub-steps,
  10 ms control period).
- **Koopman-EDMD.** State lifted by the basic rigid-body dictionary
  (`N_z = 40`); a constant `(A, B)` fitted by least squares; MPC with an
  integral channel for offset-free tracking.
- **SINDy.** Sparse identification (STLSQ) over a 48-feature physics-informed
  library that includes the bilinear thrust atom `f_t Rᵀe₃`, the
  quaternion-kinematics products and the gyroscopic products. Derivatives are
  obtained by finite differencing of the simulated states. A deviation-form
  LTV-MPC linearizes the model along the reference, so each step solves only
  a quadratic program.

## Citation

If you use this code, please cite the paper. *(Citation details to be added.)*
