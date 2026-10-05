%% ========================================================================
%  main.m -- reproduce the simulation results of
%  "Physics-Informed Sparse-Identification MPC for Agile Quadrotors:
%   Beyond Lifted-Linear Koopman Predictors"
%  ------------------------------------------------------------------------
%  Runs the full pipeline from scratch -- training-data generation, model
%  identification (Koopman-EDMD and physics-informed SINDy), and closed-loop
%  MPC -- and regenerates every table and data figure in the paper:
%
%    [1] run_diagnostics       Tables I and II (Section III-E): predicted vs
%                              fitted beta_z on datasets D1-D4, and the
%                              irreducible error floor on D4
%    [2] run_comparison        Table III: trains Koopman-EDMD and SINDy, runs
%                              both closed-loop controllers on the aggressive
%                              randomized reference, prints open-loop nRMSE,
%                              tracking RMSE and online compute time, and
%                              saves results/comparison_results.mat
%    [3] make_paper_figures    Fig. 3  results/states.png
%                              Fig. 4  results/control_inputs.png
%    [4] run_figure8           Figure-eight with three controllers -- SINDy,
%                              K-narrow (beta_z < 0), K-broad (beta_z > 0) --
%                              saves results/figure8_results.mat
%    [5] make_figure8_plots    Fig. 5a results/figure8_traj.png
%                              Fig. 5b results/figure8_err.png
%    [6] run_bilinear_ablation Table IV: EDMD with and without the bilinear
%                              thrust atom, trained on D2
%
%  Requirements: MATLAB (R2021b or newer)
%                Optimization Toolbox                 (quadprog)
%                Statistics and Machine Learning Tbx  (mvnrnd, RandStream)
%
%  Usage:   >> main
%  All outputs are written to ./results/. Online compute times depend on
%  the machine; every other number is deterministic (fixed random seeds).
%
%  Note: each stage is a self-contained script that clears the workspace,
%  so stages communicate only through the .mat files in results/.
%  ========================================================================
clc; clear; close all;
cd(fileparts(mfilename('fullpath')));        % run from the repository root

fprintf('\n===== [1/6] Tables I-II: beta_z diagnostic and error floor =====\n');
run_diagnostics

fprintf('\n===== [2/6] Table III: SINDy-LTV-MPC vs Koopman-MPC =====\n');
run_comparison

fprintf('\n===== [3/6] Figs. 3-4: state and input figures =====\n');
make_paper_figures

fprintf('\n===== [4/6] Figure-eight tracking with three controllers =====\n');
run_figure8

fprintf('\n===== [5/6] Fig. 5: figure-eight panels =====\n');
make_figure8_plots

fprintf('\n===== [6/6] Table IV: bilinear thrust atom ablation =====\n');
run_bilinear_ablation

fprintf('\n===== DONE. All results written to %s%sresults =====\n', pwd, filesep);
