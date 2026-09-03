# Physics-Informed SINDy-LTV-MPC for Agile Quadrotor Trajectory Tracking

Code to reproduce the simulation results of the paper. It compares a
**physics-informed SINDy** model, identified directly in the state space and
driven by a linear-time-varying MPC, against a **lifted-linear Koopman/EDMD**
predictor on aggressive quadrotor trajectory tracking.

The paper's central point is a *structural* limitation of the lifted-linear
predictor: because quadrotor thrust acts along the body axis, its effect on the
state depends on attitude, and a **constant** Koopman input matrix can only hold
the *average* thrust direction over the training data — a defect no richer
dictionary can remove, and one that can be read off a single identified
coefficient before flight. The SINDy model carries the bilinear thrust term
explicitly and tracks where the Koopman baselines diverge.

## Requirements

- MATLAB R2021b or newer
- Optimization Toolbox (`quadprog`)
- Statistics and Machine Learning Toolbox (`mvnrnd`, `RandStream`)

No external MPC solver (e.g. CVX) is required.

## Quick start

```matlab
>> main
```

`main.m` runs the full pipeline from scratch — training-data generation, model
identification, and closed-loop MPC — and writes every figure to `results/`.
A complete run takes a few minutes.

## What `main.m` produces

The stages run in the order the figures appear in the paper:

| Stage | Script | Output |
|-------|--------|--------|
| 1 | `run_comparison.m` | `results/comparison_results.mat`; prints open-loop nRMSE, closed-loop tracking RMSE, and per-step compute time |
| 2 | `make_paper_figures.m` | `results/states.png` (state tracking), `results/control_inputs.png` (inputs) |
| 3 | `run_figure8.m` | `results/figure8_results.mat`; figure-eight with SINDy, K-narrow (β_z < 0), K-broad (β_z > 0) |
| 4 | `make_figure8_plots.m` | `results/figure8_traj.png` (trajectory), `results/figure8_err.png` (error vs. time) |

Stages communicate only through the `.mat` files, so any stage can also be run
on its own once its input `.mat` exists.

## Repository layout

```
main.m                 orchestrator: runs the whole pipeline in paper order
run_comparison.m       train Koopman-EDMD + SINDy, closed-loop tracking, metrics
make_paper_figures.m   render states.png and control_inputs.png
run_figure8.m          figure-eight experiment (three controllers)
make_figure8_plots.m   render the two figure-eight panels
utils/                 quadrotor dynamics, RK4, quaternion/SO(3) helpers,
                       differentially-flat figure-eight reference
koopman/               EDMD observables and reference linearization
sindy/                 candidate library, STLSQ, SINDy RHS and Jacobian
results/               generated figures and data (created by `main`)
```

## Model and controller, in brief

- **Plant.** Rigid-body quadrotor, state `x = [p; v; q; ω]` (13-D), input
  `u = [f_t; M]`; integrated with fixed-step RK4.
- **Koopman-EDMD.** State lifted by a rigid-body dictionary; a constant
  `(A, B)` fitted by least squares; TV-MPC with an integral channel for
  offset-free tracking.
- **SINDy.** Sparse identification (STLSQ) over a physics-informed library that
  includes the bilinear thrust atom `f_t Rᵀe₃` and the quaternion-kinematics
  terms; a deviation-form LTV-MPC re-linearizes the model along the reference,
  so each step solves only a quadratic program.

## Citation

If you use this code, please cite the paper. *(Citation details to be added.)*
