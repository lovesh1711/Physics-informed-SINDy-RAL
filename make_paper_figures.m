%% ========================================================================
%  make_paper_figures.m
%  ------------------------------------------------------------------------
%  Renders the two figures used in the RA-L manuscript in the visual style
%  of the SE(3) Koopman-MPC paper (Narayanan et al., 2023):
%      * axis labels placed INSIDE each subplot to save space
%      * compact tile spacing
%      * tick labels visible on every panel
%      * a single legend in the top-left tile only
%
%  Run order:
%      >> run_comparison        % generates results/comparison_results.mat
%      >> make_paper_figures    % re-renders the two PNGs in paper style
%
%  Output:
%      results/states.png          (6 x 2  : states tracking)
%      results/control_inputs.png  (2 x 2  : f_t, M_1, M_2, M_3)
%  ========================================================================

clc; clear; close all;

% ---------- load saved comparison data ------------------------------------
thisDir    = fileparts(mfilename('fullpath'));
resultsDir = fullfile(thisDir, 'results');
S = load(fullfile(resultsDir, 'comparison_results.mat'));

tvec = S.tvec;                     % input times: u_k applied from t=(k-1)*ts
Nsim = numel(tvec);
ts_  = tvec(2) - tvec(1);
tst  = (0:Nsim) * ts_;             % state times: x(t_0) ... x(t_Nsim)
% logged state k is at t = k*ts; prepend the common initial state x0 so
% every state trace and the reference are drawn on the same time grid
withx0 = @(X, rows) [S.x0(rows), X];

% ---------- continuous (unwrapped) attitude rotation vectors --------------
%   The principal log vee(log R) wraps whenever the rotation angle crosses
%   pi, producing 2*pi*n jumps in the per-axis traces. We instead accumulate
%   the small relative-rotation logs between consecutive samples, which is
%   jump-free by construction. Recomputed from the saved quaternion series.
addpath(fullfile(thisDir, 'utils'), '-begin');

theta_ref   = rotvec_unwrap(S.q_ref(:, 1:Nsim+1));
theta_sindy = rotvec_unwrap(withx0(S.q_sindy, 7:10));
theta_koop  = rotvec_unwrap(withx0(S.q_koop,  7:10));

% ---------- style ---------------------------------------------------------
%  Three-color publication palette (Okabe-Ito-inspired, colour-blind safe):
%    reference: dark grey   -- sits visually 'behind' the data
%    SINDy    : vermillion  -- our proposed method (dashed)
%    Koopman  : deep blue   -- baseline (dash-dot)
c_ref   = [0.35  0.35  0.35];
c_sindy = [0.835 0.369 0.000];
c_koop  = [0.000 0.447 0.698];

lw       = 2.0;
lw_axis  = 0.1;

fs_tick  = 10;
fs_label = 13;
fs_leg   = 10;

% Common helper to place an unboxed, bold in-axis label
inLabel = @(ax, str) text(ax, 0.05, 0.78, str, ...
    'Units','normalized', 'Interpreter','tex', ...
    'FontSize', fs_label, 'FontWeight','bold', ...
    'BackgroundColor', [1 1 1], 'Margin', 1, ...
    'EdgeColor', 'none');

%% ===================== FIGURE 1: STATES (6 x 2) =========================
fig1 = figure('Name','States (paper style)','Color',[1 1 1], ...
              'Units','centimeters','Position',[2 2 14.3 21.2]);
tl1 = tiledlayout(fig1, 6, 2, ...
                  'TileSpacing','compact','Padding','compact');

% Left column data  : p_x, p_y, p_z, v_x, v_y, v_z
left_ref   = [S.p_ref_ctrl(:,1:Nsim+1); S.v_ref_ctrl(:,1:Nsim+1)];
left_sindy = [withx0(S.p_sindy,1:3);    withx0(S.v_sindy,4:6)];
left_koop  = [withx0(S.p_koop,1:3);     withx0(S.v_koop,4:6)];
left_lbl   = {'x','y','z','v_x','v_y','v_z'};

% Right column data : theta_x, theta_y, theta_z, omega_x, omega_y, omega_z
right_ref   = [theta_ref;   S.omega_ref_ctrl(:,1:Nsim+1)];
right_sindy = [theta_sindy; withx0(S.omega_sindy,11:13)];
right_koop  = [theta_koop;  withx0(S.omega_koop,11:13)];
right_lbl   = {'\theta_x','\theta_y','\theta_z','\omega_x','\omega_y','\omega_z'};

for r = 1:6
    % --- Left column tile ---
    ax = nexttile(tl1, (r-1)*2 + 1);
    plot(ax, tst, left_ref(r,:),   '-',   'Color', c_ref,   'LineWidth', lw); hold(ax,'on');
    plot(ax, tst, left_sindy(r,:), '--',  'Color', c_sindy, 'LineWidth', lw);
    plot(ax, tst, left_koop(r,:),  '-.',  'Color', c_koop,  'LineWidth', lw);
    grid(ax,'on'); box(ax,'on');
    set(ax,'FontSize',fs_tick,'LineWidth',lw_axis,'FontWeight','bold');
    if r < 6
        set(ax,'XTickLabel',[]);
    else
        xlabel(ax,'t (s)','FontSize',fs_label,'FontWeight','bold');
    end
    inLabel(ax, left_lbl{r});
    if r == 1
        lg = legend(ax, {'Reference','SINDy','Koopman'}, ...
                    'Location','southwest','FontSize',fs_leg, ...
                    'Box','on');
        lg.ItemTokenSize = [14 8];
    end

    % --- Right column tile ---
    ax = nexttile(tl1, (r-1)*2 + 2);
    plot(ax, tst, right_ref(r,:),   '-',   'Color', c_ref,   'LineWidth', lw); hold(ax,'on');
    plot(ax, tst, right_sindy(r,:), '--',  'Color', c_sindy, 'LineWidth', lw);
    plot(ax, tst, right_koop(r,:),  '-.',  'Color', c_koop,  'LineWidth', lw);
    grid(ax,'on'); box(ax,'on');
    set(ax,'FontSize',fs_tick,'LineWidth',lw_axis,'FontWeight','bold');
    if r < 6
        set(ax,'XTickLabel',[]);
    else
        xlabel(ax,'t (s)','FontSize',fs_label,'FontWeight','bold');
    end
    inLabel(ax, right_lbl{r});
end

exportgraphics(fig1, fullfile(resultsDir,'states.png'), ...
               'Resolution', 400);

%% ===================== FIGURE 2: INPUTS (2 x 2) =========================
fig2 = figure('Name','Inputs (paper style)','Color',[1 1 1], ...
              'Units','centimeters','Position',[2 2 20.3 14.6]);
tl2 = tiledlayout(fig2, 2, 2, ...
                  'TileSpacing','compact','Padding','compact');

%  Layout: (f_t, M_1) on top, (M_2, M_3) below.
in_sindy = {S.u_sindy(1,:), S.u_sindy(2,:); ...
            S.u_sindy(3,:), S.u_sindy(4,:)};
in_koop  = {S.u_koop(1,:),  S.u_koop(2,:); ...
            S.u_koop(3,:),  S.u_koop(4,:)};
in_lbl   = {'f_t','M_1';
            'M_2','M_3'};

for r = 1:2
    for c = 1:2
        ax = nexttile(tl2);
        plot(ax, tvec, in_sindy{r,c}, '-', 'Color', c_sindy, 'LineWidth', lw); hold(ax,'on');
        plot(ax, tvec, in_koop{r,c},  '-', 'Color', c_koop,  'LineWidth', lw);
        grid(ax,'on'); box(ax,'on');
        set(ax,'FontSize',fs_tick,'LineWidth',lw_axis,'FontWeight','bold');
        if r < 2
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'t (s)','FontSize',fs_label,'FontWeight','bold');
        end
        inLabel(ax, in_lbl{r,c});
        if r == 1 && c == 1
            lg = legend(ax, {'SINDy','Koopman'}, ...
                        'Location','southeast','FontSize',fs_leg, ...
                        'Box','on');
            lg.ItemTokenSize = [14 8];
        end
    end
end

exportgraphics(fig2, fullfile(resultsDir,'control_inputs.png'), ...
               'Resolution', 400);

fprintf('\nPaper figures written to:\n  %s\n  %s\n', ...
        fullfile(resultsDir,'states.png'), ...
        fullfile(resultsDir,'control_inputs.png'));


%% ====================== LOCAL FUNCTIONS ==================================
function th = rotvec_unwrap(Q)
%ROTVEC_UNWRAP  Continuous rotation-vector trace from a quaternion series.
%   Q  : 4 x N quaternion time series (each column [qw;qx;qy;qz]).
%   th : 3 x N continuous axis-angle vector, accumulated from the small
%        relative-rotation logs between consecutive samples (jump-free).
    N  = size(Q, 2);
    th = zeros(3, N);
    Rprev    = QuatToRot(normalize_quat(Q(:,1)));
    th(:,1)  = vee(so3_log(Rprev));
    for k = 2:N
        Rk      = QuatToRot(normalize_quat(Q(:,k)));
        dR      = Rk * Rprev.';                 % relative rotation, small
        th(:,k) = th(:,k-1) + vee(so3_log(dR)); % accumulate the increment
        Rprev   = Rk;
    end
end
