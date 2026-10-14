
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
