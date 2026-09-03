%% ========================================================================
%  main.m -- reproduce the simulation figures of
%  "Physics-Informed SINDy-LTV-MPC for Agile Quadrotor Trajectory Tracking"
%  ------------------------------------------------------------------------
%  Runs the full pipeline from scratch -- training-data generation, model
%  identification (Koopman-EDMD and physics-informed SINDy), and closed-loop
%  MPC -- and regenerates every data figure in the paper, in the order the
%  figures appear:
%
%    [1] run_comparison      Trains the Koopman-EDMD and SINDy models, runs
%                            both closed-loop controllers on the aggressive
%                            randomized reference, prints the open-loop /
%                            tracking / timing numbers, and saves
%                            results/comparison_results.mat
%    [2] make_paper_figures  Renders results/states.png          (Fig. 3)
%                            and     results/control_inputs.png  (Fig. 4)
%    [3] run_figure8         Runs the figure-eight experiment with three
%                            controllers -- SINDy, K-narrow (beta_z < 0) and
%                            K-broad (beta_z > 0) -- and saves
%                            results/figure8_results.mat
%    [4] make_figure8_plots  Renders results/figure8_traj.png    (Fig. 5a)
%                            and     results/figure8_err.png     (Fig. 5b)
%
%  Requirements: MATLAB (R2021b or newer)
%                Optimization Toolbox                 (quadprog)
%                Statistics and Machine Learning Tbx  (mvnrnd, RandStream)
%
%  Usage:   >> main
%  All outputs are written to ./results/.
%
%  Note: each stage below is a self-contained script that clears the
%  workspace and saves its results to disk, so the stages communicate only
%  through the .mat files in results/. A full run takes a few minutes.
%  ========================================================================
clc; clear; close all;
cd(fileparts(mfilename('fullpath')));        % run from the repository root

fprintf('\n===== [1/4] Closed-loop tracking: SINDy-LTV-MPC vs Koopman-MPC =====\n');
run_comparison

fprintf('\n===== [2/4] Rendering state and input figures =====\n');
make_paper_figures

fprintf('\n===== [3/4] Figure-eight tracking with three controllers =====\n');
run_figure8

fprintf('\n===== [4/4] Rendering figure-eight panels =====\n');
make_figure8_plots

fprintf('\n===== DONE. All figures written to %s%sresults =====\n', pwd, filesep);
