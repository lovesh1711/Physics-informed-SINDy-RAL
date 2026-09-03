%% make_figure8_plots.m --- figure-8 figures from figure8_results.mat
% Emits TWO separate panels (each with its own top legend) so they can be
% dropped into LaTeX as subfigures:
%   results/figure8_traj.png  -- trajectory in the North-East plane
%   results/figure8_err.png   -- position error vs time (log)
clc; clear; close all;
thisDir=fileparts(mfilename('fullpath'));
L=load(fullfile(thisDir,'results','figure8_results.mat'));
N=L.Nsim; t=L.tvec(1:N);
% match states.png (make_paper_figures.m): reference grey solid,
% SINDy vermillion dashed, Koopman blue dash-dot. K-narrow keeps the states
% Koopman blue; K-broad gets a distinct purple. Line styles matched below.
cRef=[0.35 0.35 0.35]; cS=[0.835 0.369 0.000];
cN=[0.000 0.447 0.698]; cB=[0.494 0.184 0.556];

pe=@(X) sqrt(sum((X(1:3,1:N)-L.Xref(1:3,1:N)).^2,1));

%% ---- Panel (a): xy trajectory, zoomed to the figure-8 ------------------
figA=figure('Color','w','Position',[80 80 500 500]);
axA=axes(figA); hold(axA,'on');
h0=plot(L.Xref(2,1:N),L.Xref(1,1:N),'-', 'Color',cRef,'LineWidth',2.4);
h1=plot(L.Xs(2,1:N), L.Xs(1,1:N), '--','Color',cS,'LineWidth',2.2);
h2=plot(L.Xk_n(2,1:N),L.Xk_n(1,1:N),'-.','Color',cN,'LineWidth',2.0);
h3=plot(L.Xk_b(2,1:N),L.Xk_b(1,1:N),'-.','Color',cB,'LineWidth',2.0);
axis(axA,'equal'); xlim(axA,[-3 3]); ylim(axA,[-3 3]); grid(axA,'on'); box(axA,'on');
xlabel(axA,'East  y [m]'); ylabel(axA,'North  x [m]');
set(axA,'FontSize',12,'LineWidth',1.1,'GridAlpha',0.15);
lgA=legend(axA,[h0 h1 h2 h3], ...
    {'reference','SINDy','K-narrow (\beta_z<0)','K-broad (\beta_z>0)'}, ...
    'Orientation','horizontal','NumColumns',2,'Box','off','FontSize',11, ...
    'Location','northoutside');
exportgraphics(figA,fullfile(thisDir,'results','figure8_traj.png'),'Resolution',240);

%% ---- Panel (b): position error vs time (log) --------------------------
figB=figure('Color','w','Position',[80 80 560 380]);
axB=axes(figB); hold(axB,'on');
g1=semilogy(axB,t,pe(L.Xs),'--','Color',cS,'LineWidth',2.2);
g2=semilogy(axB,t,pe(L.Xk_n),'-.','Color',cN,'LineWidth',2.0);
g3=semilogy(axB,t,pe(L.Xk_b),'-.','Color',cB,'LineWidth',2.0);
set(axB,'YScale','log');
grid(axB,'on'); box(axB,'on');
xlabel(axB,'time [s]'); ylabel(axB,'position error  |p - p_r|  [m]');
ylim(axB,[1e-2 3e2]);
set(axB,'FontSize',12,'LineWidth',1.1,'GridAlpha',0.15);
lgB=legend(axB,[g1 g2 g3], ...
    {'SINDy','K-narrow (\beta_z<0)','K-broad (\beta_z>0)'}, ...
    'Orientation','horizontal','NumColumns',3,'Box','off','FontSize',11, ...
    'Location','northoutside');
exportgraphics(figB,fullfile(thisDir,'results','figure8_err.png'),'Resolution',240);

fprintf('saved figure8_traj.png and figure8_err.png\n');

% also report when each Koopman error exceeds 1 m (divergence onset)
en=pe(L.Xk_n); eb=pe(L.Xk_b);
in=find(en>1,1); ib=find(eb>1,1);
fprintf('divergence onset (err>1m): K-narrow t=%.2fs, K-broad t=%.2fs\n', t(in), t(ib));
fprintf('altitude drift at end: K-narrow z=%.1f, K-broad z=%.1f, SINDy z=%.2f (ref z=0)\n', ...
        L.Xk_n(3,N), L.Xk_b(3,N), L.Xs(3,N));
fprintf('DONE\n');
