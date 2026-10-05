%% ========================================================================
%  run_bilinear_ablation.m  --  Table IV of the paper (Section IV-E)
%  ------------------------------------------------------------------------
%  Adds the bilinear thrust atom f_t * R(q)' * e3 to the EDMD input and
%  compares open-loop prediction accuracy against plain EDMD:
%
%      (1) EDMD            : u     = [f_t; M_1; M_2; M_3]          (4-dim)
%      (2) EDMD + atom     : u_aug = [u; f_t * R(q)' e3]            (7-dim)
%
%  Both use the basic rigid-body dictionary (Nz = 40) and are trained on
%  the broad-attitude dataset D2 (same recipe and seed as run_diagnostics.m)
%  and validated on a separate broad-attitude set (seed 12, 50 x 0.5 s).
%  Also prints the fitted gain of the bilinear channel against -Ts/m.
%  ========================================================================

clc; clear; close all;

% ---------- setup ---------------------------------------------------------
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'utils'),   '-begin');
addpath(fullfile(thisDir, 'koopman'), '-begin');

sTrain = RandStream('mt19937ar','Seed',11);   % D2 seed
sVal   = RandStream('mt19937ar','Seed',12);   % separate validation seed

params.m = 4.34;
params.J = diag([0.0820; 0.0845; 0.1377]);
params.g = 9.81;

u_min = [10;  -30; -30; -30];
u_max = [300;  30;  30;  30];

ts    = 0.01;
dt    = 0.001;
steps = round(ts/dt);
n     = 13;
e3    = [0;0;1];

% ---------- training data: broad-attitude regime (D2 of Section III-E) ----
%   600 rollouts of 1 s from randomly tilted attitudes (tilt up to pi/4),
%   inputs ~ N(0, 15 I) + mg on thrust, held 3 or 7 sample periods
%   (alternating), seed 11 -- identical to D2 and to K-broad (run_figure8.m).
Sigma_u = diag([15; 15; 15; 15]);
[X_list, Xp_list, U_list] = gen_broad(sTrain, 600, 1.0, Sigma_u, params, ts, dt, steps, u_min, u_max);
X    = X_list;
Xp   = Xp_list;
Uall = U_list;
K    = size(X, 2);
fprintf('Training data: %d snapshots\n', K);

% ---------- shared lifting ------------------------------------------------
p_lift = 3;
Nz     = n + 9*p_lift;   % = 40

Z  = zeros(Nz, K);
Zp = zeros(Nz, K);
for k = 1:K
    xk  = X(:,k);   xk(7:10)  = normalize_quat(xk(7:10));
    xkp = Xp(:,k);  xkp(7:10) = normalize_quat(xkp(7:10));
    Z(:,k)  = observable_phi(xk,  p_lift);
    Zp(:,k) = observable_phi(xkp, p_lift);
end

% ---------- augmented input: u_aug = [u; f_t * R^T e_3] -------------------
%   At every training snapshot, compute the bilinear atom from the
%   contemporary state and input.
Uaug = zeros(7, K);
for k = 1:K
    qk    = normalize_quat(X(7:10, k));
    Rk    = QuatToRot(qk);
    b3_k  = Rk' * e3;
    f_t_k = Uall(1, k);
    Uaug(:, k) = [Uall(:, k); f_t_k * b3_k];
end

% ---------- train both EDMD models ----------------------------------------
lambda_edmd = 1e-6;

% (1) Baseline
ZU       = [Z; Uall];
G_base   = Zp * ZU' / (ZU * ZU' + lambda_edmd * eye(size(ZU,1)));
A_base   = G_base(:, 1:Nz);
B_base   = G_base(:, Nz+1:end);

% (2) Augmented input
ZUa      = [Z; Uaug];
G_aug    = Zp * ZUa' / (ZUa * ZUa' + lambda_edmd * eye(size(ZUa,1)));
A_aug    = G_aug(:, 1:Nz);
B_aug    = G_aug(:, Nz+1:end);

fprintf('Baseline EDMD : A %dx%d, B %dx%d\n',  size(A_base,1),size(A_base,2),size(B_base,1),size(B_base,2));
fprintf('Augmented EDMD: A %dx%d, B %dx%d  (input dim = %d)\n', ...
        size(A_aug,1),size(A_aug,2),size(B_aug,1),size(B_aug,2), size(B_aug,2));

% ---------- open-loop validation ------------------------------------------
Nval_traj  = 50;
Nsteps_val = round(0.5 / ts);
RandStream.setGlobalStream(sVal);

C_state = [eye(n), zeros(n, 9*p_lift)];   % project lifted -> state

num_p_b = 0; den_p_b = 0;  num_v_b = 0; den_v_b = 0;
num_th_b = 0; den_th_b = 0; num_w_b = 0; den_w_b = 0;
num_p_a = 0; den_p_a = 0;  num_v_a = 0; den_v_a = 0;
num_th_a = 0; den_th_a = 0; num_w_a = 0; den_w_a = 0;

for tr = 1:Nval_traj
    x_true = broad_ic();                       % same broad initial conditions
    z_b = observable_phi(x_true, p_lift);
    z_a = observable_phi(x_true, p_lift);
    u   = [params.m*params.g; 0; 0; 0];
    if mod(tr,2)==0, Nu_hold = 3; else, Nu_hold = 7; end

    for k = 1:Nsteps_val
        if (k == 1) || (mod(k-1, Nu_hold) == 0)
            u = mvnrnd(zeros(4,1), Sigma_u).';
            u(1) = u(1) + params.m*params.g;
            u = min(max(u, u_min), u_max);
        end

        % --- True next state
        for j = 1:steps, x_true = rk4_step(x_true, u, params, dt); end
        x_true(7:10) = normalize_quat(x_true(7:10));

        % --- Baseline Koopman next state
        z_b = A_base * z_b + B_base * u;
        x_b = C_state * z_b;

        % --- Augmented Koopman next state (uses x_a current to build aug input)
        x_a_curr = C_state * z_a;
        x_a_curr(7:10) = normalize_quat(x_a_curr(7:10));
        R_a = QuatToRot(x_a_curr(7:10));
        u_aug = [u; u(1) * (R_a' * e3)];
        z_a = A_aug * z_a + B_aug * u_aug;
        x_a = C_state * z_a;

        % --- Accumulate squared errors per block (1-indexed: p, v, th, w)
        % position
        num_p_b = num_p_b + sum((x_b(1:3)  - x_true(1:3)).^2);
        num_p_a = num_p_a + sum((x_a(1:3)  - x_true(1:3)).^2);
        den_p_b = den_p_b + sum(x_true(1:3).^2);
        den_p_a = den_p_b;
        % velocity
        num_v_b = num_v_b + sum((x_b(4:6)  - x_true(4:6)).^2);
        num_v_a = num_v_a + sum((x_a(4:6)  - x_true(4:6)).^2);
        den_v_b = den_v_b + sum(x_true(4:6).^2);
        den_v_a = den_v_b;
        % attitude (raw quaternion error norm -- crude but consistent for both)
        num_th_b = num_th_b + sum((x_b(7:10) - x_true(7:10)).^2);
        num_th_a = num_th_a + sum((x_a(7:10) - x_true(7:10)).^2);
        den_th_b = den_th_b + sum(x_true(7:10).^2);
        den_th_a = den_th_b;
        % angular velocity
        num_w_b = num_w_b + sum((x_b(11:13) - x_true(11:13)).^2);
        num_w_a = num_w_a + sum((x_a(11:13) - x_true(11:13)).^2);
        den_w_b = den_w_b + sum(x_true(11:13).^2);
        den_w_a = den_w_b;
    end
end

eps_d = 1e-6;
nRMSE_b_p  = 100 * sqrt(num_p_b  / max(den_p_b,  eps_d));
nRMSE_b_v  = 100 * sqrt(num_v_b  / max(den_v_b,  eps_d));
nRMSE_b_th = 100 * sqrt(num_th_b / max(den_th_b, eps_d));
nRMSE_b_w  = 100 * sqrt(num_w_b  / max(den_w_b,  eps_d));
nRMSE_a_p  = 100 * sqrt(num_p_a  / max(den_p_a,  eps_d));
nRMSE_a_v  = 100 * sqrt(num_v_a  / max(den_v_a,  eps_d));
nRMSE_a_th = 100 * sqrt(num_th_a / max(den_th_a, eps_d));
nRMSE_a_w  = 100 * sqrt(num_w_a  / max(den_w_a,  eps_d));

fprintf('\n================ Open-loop nRMSE comparison ================\n');
fprintf('%-14s %12s %12s %12s\n', 'State', 'Baseline', 'Augmented', 'Change');
fprintf('%-14s %11.2f%% %11.2f%% %+11.2f%%\n', ...
        'Position',          nRMSE_b_p,  nRMSE_a_p,  nRMSE_a_p  - nRMSE_b_p);
fprintf('%-14s %11.2f%% %11.2f%% %+11.2f%%\n', ...
        'Velocity',          nRMSE_b_v,  nRMSE_a_v,  nRMSE_a_v  - nRMSE_b_v);
fprintf('%-14s %11.2f%% %11.2f%% %+11.2f%%\n', ...
        'Attitude (quat)',   nRMSE_b_th, nRMSE_a_th, nRMSE_a_th - nRMSE_b_th);
fprintf('%-14s %11.2f%% %11.2f%% %+11.2f%%\n', ...
        'Angular vel',       nRMSE_b_w,  nRMSE_a_w,  nRMSE_a_w  - nRMSE_b_w);
fprintf('============================================================\n');

% ---------- inspect the new B columns (sanity / interpretation) ------------
%   In the augmented model the last 3 input columns map (f_t * b3) into the
%   lifted dynamics.  Print the row of B_aug that hits velocity v_z to see if
%   EDMD found the analytic -t_s/m sign.
row_vz = 6;   % v_z is state index 6, also lifted index 6
fprintf('\nAnalytic Euler v_z thrust gain: -t_s/m = %+.3e\n', -ts/params.m);
fprintf('Baseline  B(v_z, f_t)               = %+.3e\n', B_base(row_vz, 1));
fprintf('Augmented B(v_z, f_t)               = %+.3e   (direct channel)\n', ...
        B_aug(row_vz, 1));
fprintf('Augmented B(v_z, f_t*b3_z) column   = %+.3e   (bilinear channel)\n', ...
        B_aug(row_vz, 7));


%% ---------- local functions -------------------------------------------
function x = broad_ic()
    % random tilted initial state (D2 recipe of Section III-E)
    p0 = 1.0*(2*rand(3,1)-1);  v0 = 2.0*(2*rand(3,1)-1);
    ax = randn(3,1); ax = ax/max(norm(ax),1e-9); ang = (pi/4)*rand;
    q0 = normalize_quat([cos(ang/2); sin(ang/2)*ax]);
    w0 = 3.0*(2*rand(3,1)-1);
    x  = [p0; v0; q0; w0];
end

function [X, Xp, U] = gen_broad(stream, Ntraj, T, Sigma_u, params, ts, dt, steps, u_min, u_max)
    RandStream.setGlobalStream(stream);
    n = 13; Ns = round(T/ts); X = zeros(n, Ntraj*Ns); Xp = X; U = zeros(4, Ntraj*Ns); kk = 0;
    mg = params.m*params.g;
    for tr = 1:Ntraj
        x = broad_ic();
        if mod(tr,2)==0, hl = 3; else, hl = 7; end
        u = [mg; 0; 0; 0];
        for k = 1:Ns
            if (k==1) || (mod(k-1,hl)==0)
                u = mvnrnd(zeros(4,1), Sigma_u).'; u(1) = u(1) + mg;
                u = min(max(u, u_min), u_max);
            end
            xk = x;
            for j = 1:steps, x = rk4_step(x, u, params, dt); end
            kk = kk+1; X(:,kk) = xk; Xp(:,kk) = x; U(:,kk) = u;
        end
    end
end
