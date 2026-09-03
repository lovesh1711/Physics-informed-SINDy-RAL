%% ========================================================================
%  run_comparison.m  —  Koopman-EDMD vs SINDy-C for Quadrotor MPC
%  ========================================================================
%  Self-contained comparison script for the paper.
%  Run section-by-section (Ctrl+Enter in MATLAB) or all at once.
%
%  Sections:
%    1. Setup (paths, params, seeds)
%    2. Shared training data generation
%    3. Koopman EDMD training + open-loop validation
%    4. SINDy-C training + open-loop validation
%    5. Reference trajectory generation
%    6. Koopman TV-MPC (with integral action)
%    7. SINDy LTV-MPC
%    8. Comparison plots + metrics
%    9. Save results
%  ========================================================================

clc; clear; close all;

%% ======================== 1. SETUP ====================================
% Add subfolders to path (use -begin so paper/ versions take priority)
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'utils'),   '-begin');
addpath(fullfile(thisDir, 'koopman'), '-begin');
addpath(fullfile(thisDir, 'sindy'),   '-begin');

% Reproducible streams (identical data for both methods)
sTrain = RandStream('mt19937ar', 'Seed', 1);
sVal   = RandStream('mt19937ar', 'Seed', 2);
sRef   = RandStream('mt19937ar', 'Seed', 3);

% Quadrotor parameters
params.m = 4.34;
params.J = diag([0.0820; 0.0845; 0.1377]);
params.g = 9.81;

% Input bounds
u_min = [10;  -30; -30; -30];
u_max = [300;  30;  30;  30];

% Timing
ts    = 0.01;           % EDMD / MPC sample period (s)
dt    = 0.001;          % RK4 integration step (s)
steps = round(ts/dt);   % sub-steps per sample period

% State dimensions
n  = 13;                % [p(3); v(3); q(4); omega(3)]
nu = 4;                 % [ft; M1; M2; M3]

fprintf('=== Setup complete ===\n');

%% ======================== 2. TRAINING DATA ============================
T_train = 0.5;          % trajectory length (s)
Ntraj   = 100;          % number of training trajectories
Nsteps_train = round(T_train / ts);

Nu_hold = 5;            % piecewise-constant input hold (steps)
mu_u    = zeros(4, 1);
Sigma_u = diag([10; 10; 10; 10]);

X_list    = {};
Xp_list   = {};
U_list    = {};
Xdot_list = {};          % also store derivatives (needed by SINDy)

RandStream.setGlobalStream(sTrain);

for tr = 1:Ntraj
    x = [zeros(3,1); zeros(3,1); [1;0;0;0]; zeros(3,1)];
    x(7:10) = normalize_quat(x(7:10));
    u = [params.m*params.g; 0; 0; 0];

    for k = 1:Nsteps_train
        if (k == 1) || (mod(k-1, Nu_hold) == 0)
            u = mvnrnd(mu_u, Sigma_u).';
            u(1) = u(1) + params.m*params.g;
            u = min(max(u, u_min), u_max);
        end

        x_k    = x;
        xdot_k = quad_dynamics(0, x_k, u, params);   % continuous derivative

        for j = 1:steps
            x = rk4_step(x, u, params, dt);
        end
        x_kp1 = x;

        X_list{end+1}    = x_k;
        Xp_list{end+1}   = x_kp1;
        U_list{end+1}    = u;
        Xdot_list{end+1} = xdot_k;
    end
end

X    = cell2mat(X_list);       % 13 x K
Xp   = cell2mat(Xp_list);     % 13 x K
Uall = cell2mat(U_list);      %  4 x K
Xdot = cell2mat(Xdot_list);   % 13 x K

K_data = size(X, 2);
fprintf('Training data: %d trajectories, %d snapshots\n', Ntraj, K_data);

%% ======================== 3. KOOPMAN EDMD =============================
fprintf('\n--- Koopman EDMD training ---\n');

p_lift = 3;                         % observable order
Nz     = n + 9*p_lift;              % lifted dimension = 40

Z  = zeros(Nz, K_data);
Zp = zeros(Nz, K_data);
for k = 1:K_data
    xk  = X(:,k);   xk(7:10)  = normalize_quat(xk(7:10));
    xkp = Xp(:,k);  xkp(7:10) = normalize_quat(xkp(7:10));
    Z(:,k)  = observable_phi(xk,  p_lift);
    Zp(:,k) = observable_phi(xkp, p_lift);
end

ZU     = [Z; Uall];
lambda_edmd = 1e-6;
G_edmd = Zp * ZU' / (ZU*ZU' + lambda_edmd*eye(size(ZU,1)));

A_koop = G_edmd(:, 1:Nz);
B_koop = G_edmd(:, Nz+1:end);
C_koop = [eye(n), zeros(n, 9*p_lift)];

fprintf('Koopman: A (%dx%d), B (%dx%d), lifted dim = %d\n', ...
    size(A_koop,1), size(A_koop,2), size(B_koop,1), size(B_koop,2), Nz);

% --- Koopman open-loop validation ---
Nval_traj  = 50;
T_val      = 0.5;
Nsteps_val = round(T_val / ts);
mu_val     = zeros(4, 1);
Sigma_val  = diag([10; 10; 10; 10]);

RandStream.setGlobalStream(sVal);

koop_num_p = 0; koop_den_p = 0;
koop_num_v = 0; koop_den_v = 0;
koop_num_th = 0; koop_den_th = 0;
koop_num_w = 0; koop_den_w = 0;
eps_den = 1e-6;

for tr = 1:Nval_traj
    x_true = [zeros(3,1); zeros(3,1); [1;0;0;0]; zeros(3,1)];
    x_true(7:10) = normalize_quat(x_true(7:10));
    z = observable_phi(x_true, p_lift);
    u = [params.m*params.g; 0; 0; 0];

    for k = 1:Nsteps_val
        if (k == 1) || (mod(k-1, Nu_hold) == 0)
            u = mvnrnd(mu_val, Sigma_val).';
            u(1) = u(1) + params.m*params.g;
            u = min(max(u, u_min), u_max);
        end

        x = x_true;
        for j = 1:steps
            x = rk4_step(x, u, params, dt);
        end
        x_true = x;

        z = A_koop*z + B_koop*u;
        x_pred = C_koop*z;
        x_pred(7:10) = normalize_quat(x_pred(7:10));

        p_err = norm(x_pred(1:3)   - x_true(1:3));
        v_err = norm(x_pred(4:6)   - x_true(4:6));
        w_err = norm(x_pred(11:13) - x_true(11:13));

        R_true = QuatToRot(normalize_quat(x_true(7:10)));
        R_hat  = QuatToRot(normalize_quat(x_pred(7:10)));
        th_err = norm(vee(so3_log(R_hat' * R_true)));
        th_tru = norm(vee(so3_log(R_true)));

        koop_num_p = koop_num_p + p_err^2;
        koop_den_p = koop_den_p + max(norm(x_true(1:3))^2, eps_den);
        koop_num_v = koop_num_v + v_err^2;
        koop_den_v = koop_den_v + max(norm(x_true(4:6))^2, eps_den);
        koop_num_th = koop_num_th + th_err^2;
        koop_den_th = koop_den_th + max(th_tru^2, eps_den);
        koop_num_w = koop_num_w + w_err^2;
        koop_den_w = koop_den_w + max(norm(x_true(11:13))^2, eps_den);
    end
end

koop_nRMSE_p  = 100*sqrt(koop_num_p  / koop_den_p);
koop_nRMSE_v  = 100*sqrt(koop_num_v  / koop_den_v);
koop_nRMSE_th = 100*sqrt(koop_num_th / koop_den_th);
koop_nRMSE_w  = 100*sqrt(koop_num_w  / koop_den_w);

fprintf('Koopman nRMSE:  pos=%.2f%%  vel=%.2f%%  theta=%.2f%%  omega=%.2f%%\n', ...
    koop_nRMSE_p, koop_nRMSE_v, koop_nRMSE_th, koop_nRMSE_w);

%% ======================== 4. SINDy-C TRAINING ==========================
fprintf('\n--- SINDy-C training ---\n');

[Theta_lib, feature_names] = sindy_library(X, Uall);
Nfeat = size(Theta_lib, 1);

lambda_sindy = 0.10;
n_iter_sindy = 10;
Xi = stlsq(Theta_lib, Xdot, lambda_sindy, n_iter_sindy);

fprintf('SINDy: n=%d, Nfeat=%d, nnz(Xi)=%d\n', n, Nfeat, nnz(Xi));

state_names = {'p1','p2','p3','v1','v2','v3','qw','qx','qy','qz','w1','w2','w3'};
print_sindy_model(Xi, feature_names, state_names, 1e-12);

% --- SINDy open-loop validation ---
RandStream.setGlobalStream(sVal);   % SAME seed as Koopman validation

sindy_num_p = 0; sindy_den_p = 0;
sindy_num_v = 0; sindy_den_v = 0;
sindy_num_th = 0; sindy_den_th = 0;
sindy_num_w = 0; sindy_den_w = 0;

for tr = 1:Nval_traj
    x_true = [zeros(3,1); zeros(3,1); [1;0;0;0]; zeros(3,1)];
    x_true(7:10) = normalize_quat(x_true(7:10));
    x_hat = x_true;
    u = [params.m*params.g; 0; 0; 0];

    for k = 1:Nsteps_val
        if (k == 1) || (mod(k-1, Nu_hold) == 0)
            u = mvnrnd(mu_val, Sigma_val).';
            u(1) = u(1) + params.m*params.g;
            u = min(max(u, u_min), u_max);
        end

        x = x_true;
        for j = 1:steps
            x = rk4_step(x, u, params, dt);
        end
        x_true = x;
        x_true(7:10) = normalize_quat(x_true(7:10));

        x_hat = rk4_step_sindy(x_hat, u, Xi, dt, steps);
        x_hat(7:10) = normalize_quat(x_hat(7:10));

        p_err = norm(x_hat(1:3)   - x_true(1:3));
        v_err = norm(x_hat(4:6)   - x_true(4:6));
        w_err = norm(x_hat(11:13) - x_true(11:13));

        R_true = QuatToRot(normalize_quat(x_true(7:10)));
        R_hat  = QuatToRot(normalize_quat(x_hat(7:10)));
        th_err = norm(vee(so3_log(R_hat' * R_true)));
        th_tru = norm(vee(so3_log(R_true)));

        sindy_num_p = sindy_num_p + p_err^2;
        sindy_den_p = sindy_den_p + max(norm(x_true(1:3))^2, eps_den);
        sindy_num_v = sindy_num_v + v_err^2;
        sindy_den_v = sindy_den_v + max(norm(x_true(4:6))^2, eps_den);
        sindy_num_th = sindy_num_th + th_err^2;
        sindy_den_th = sindy_den_th + max(th_tru^2, eps_den);
        sindy_num_w = sindy_num_w + w_err^2;
        sindy_den_w = sindy_den_w + max(norm(x_true(11:13))^2, eps_den);
    end
end

sindy_nRMSE_p  = 100*sqrt(sindy_num_p  / sindy_den_p);
sindy_nRMSE_v  = 100*sqrt(sindy_num_v  / sindy_den_v);
sindy_nRMSE_th = 100*sqrt(sindy_num_th / sindy_den_th);
sindy_nRMSE_w  = 100*sqrt(sindy_num_w  / sindy_den_w);

fprintf('SINDy nRMSE:  pos=%.2f%%  vel=%.2f%%  theta=%.2f%%  omega=%.2f%%\n', ...
    sindy_nRMSE_p, sindy_nRMSE_v, sindy_nRMSE_th, sindy_nRMSE_w);

%% ======================== 5. REFERENCE TRAJECTORY =====================
fprintf('\n--- Generating reference trajectory ---\n');

Nh    = 7;              % prediction horizon
t_sim = 1.2;            % simulation time (s)
Nsim  = round(t_sim / ts);

Nref_total = Nsim * steps;
x_ref = zeros(n, Nref_total);

xr = [zeros(3,1); zeros(3,1); [1;0;0;0]; zeros(3,1)];
xr(7:10) = normalize_quat(xr(7:10));

RandStream.setGlobalStream(sRef);
mu_ref    = [2; 2; 2; 2];
Sigma_ref = diag([30, 30, 30, 30]);

u_ref_high = mvnrnd(mu_ref, Sigma_ref, Nref_total)';
u_ref_high(1,:) = u_ref_high(1,:) + params.m*params.g;
u_ref_high = min(max(u_ref_high, u_min), u_max);

for k = 1:Nref_total
    xr = rk4_step(xr, u_ref_high(:,k), params, dt);
    xr(7:10) = normalize_quat(xr(7:10));
    x_ref(:,k) = xr;
end

ref_idx    = 1:steps:Nref_total;
x_ref_ctrl = x_ref(:, ref_idx);
Nref_ctrl  = size(x_ref_ctrl, 2);

u_ref_ctrl = zeros(nu, Nref_ctrl);
for k = 1:Nref_ctrl
    u_ref_ctrl(:,k) = u_ref_high(:, steps*(k-1)+1);
end

% Extract reference signals for plotting
p_ref_ctrl     = x_ref_ctrl(1:3, :);
v_ref_ctrl     = x_ref_ctrl(4:6, :);
omega_ref_ctrl = x_ref_ctrl(11:13, :);
theta_ref_ctrl = zeros(3, Nref_ctrl);
for k = 1:Nref_ctrl
    Rk = QuatToRot(normalize_quat(x_ref_ctrl(7:10,k)));
    theta_ref_ctrl(:,k) = vee(so3_log(Rk));
end

% Lift reference for Koopman
z_ref_ctrl = zeros(Nz, Nref_ctrl);
for k = 1:Nref_ctrl
    xk = x_ref_ctrl(:,k);
    xk(7:10) = normalize_quat(xk(7:10));
    z_ref_ctrl(:,k) = observable_phi(xk, p_lift);
end

fprintf('Reference: %d control-rate steps over %.2f s\n', Nref_ctrl, t_sim);

%% ======================== 6. KOOPMAN TV-MPC ===========================
fprintf('\n--- Running Koopman TV-MPC ---\n');

% Augmented state: [z(Nz); eta(3); 1]  (integral + constant)
Nz_aug = Nz + 3 + 1;
C_p    = [eye(3), zeros(3, n-3)] * C_koop;    % position extraction from z

% Cost weights
Q_pos = diag([500 500 500]);
Q_vel = diag([50 50 50]);
Q_q   = diag([10 10 10 10]);
Q_w   = diag([100 100 100]);
w_I   = 1;
S_int = w_I * eye(3);
eta_max = 0.5;                          % anti-windup clamp

Qtilde_koop = zeros(Nz, Nz);
Qtilde_koop(1:13, 1:13) = blkdiag(Q_pos, Q_vel, Q_q, Q_w);
Qtilde_aug = blkdiag(Qtilde_koop, S_int, 0);
Qbar_koop  = kron(eye(Nh), Qtilde_aug);

w_ft = 0.01; w_M = 5;
Ru_koop  = diag([w_ft, w_M, w_M, w_M]);
Rbar_koop = kron(eye(Nh), Ru_koop);

% Tighter input bounds (300N on a 4.34kg quad = 7g → unstable flips)
u_min_koop = [10;  -10; -10; -10];
u_max_koop = [100;  10;  10;  10];
U_lb_koop = repmat(u_min_koop, Nh, 1);
U_ub_koop = repmat(u_max_koop, Nh, 1);

% Init
x_true_koop = [zeros(3,1); zeros(3,1); [1;0;0;0]; zeros(3,1)];
x_true_koop(7:10) = normalize_quat(x_true_koop(7:10));
z   = observable_phi(x_true_koop, p_lift);
eta = zeros(3, 1);
z_aug = [z; eta; 1];

% Storage
p_koop     = zeros(3, Nsim);
v_koop     = zeros(3, Nsim);
theta_koop = zeros(3, Nsim);
q_koop     = zeros(4, Nsim);         % quaternion log (for SO(3) geodesic error)
omega_koop = zeros(3, Nsim);
u_koop     = zeros(nu, Nsim);
t_koop_mpc = zeros(1, Nsim);         % MPC solve time

opts = optimoptions('quadprog', 'Display', 'none', 'Algorithm', 'interior-point-convex');
localLin = @(zr, ur) local_koopman_linearization(zr, ur, A_koop, B_koop);

for k = 1:Nsim
    tic_k = tic;

    ref_k = min(k, Nref_ctrl);

    % Build TV prediction matrices
    Aseq = cell(Nh, 1);
    Bseq = cell(Nh, 1);
    for i = 1:Nh
        idx = min(ref_k + i - 1, Nref_ctrl);
        uref_i = u_ref_high(:, steps*(idx-1)+1);

        [Ai, Bi] = localLin(z_ref_ctrl(:,idx), uref_i);
        pref = p_ref_ctrl(:, idx);

        Aseq{i} = [Ai,           zeros(Nz,3),  zeros(Nz,1);
                    ts*C_p,       eye(3),       -ts*pref;
                    zeros(1,Nz),  zeros(1,3),    1         ];
        Bseq{i} = [Bi; zeros(3,nu); zeros(1,nu)];
    end

    Aqp = zeros(Nz_aug*Nh, Nz_aug);
    Bqp = zeros(Nz_aug*Nh, nu*Nh);
    Ap = eye(Nz_aug);
    for i = 1:Nh
        Ap = Aseq{i} * Ap;
        Aqp((i-1)*Nz_aug+1:i*Nz_aug, :) = Ap;
        for j = 1:i
            Aprod = eye(Nz_aug);
            for s = j+1:i, Aprod = Aseq{s} * Aprod; end
            Bqp((i-1)*Nz_aug+1:i*Nz_aug, (j-1)*nu+1:j*nu) = Aprod * Bseq{j};
        end
    end

    Y = zeros(Nz_aug*Nh, 1);
    for i = 1:Nh
        idx  = min(ref_k + i - 1, Nref_ctrl);
        base = (i-1)*Nz_aug;
        Y(base+1:base+Nz)      = z_ref_ctrl(:, idx);
        Y(base+Nz+1:base+Nz+3) = zeros(3,1);
        Y(base+Nz+4)           = 1;
    end

    H = 2*(Bqp'*Qbar_koop*Bqp + Rbar_koop);
    f = 2*(Bqp'*Qbar_koop*(Aqp*z_aug - Y));
    H = 0.5*(H + H') + 1e-8*eye(size(H));

    [Uopt, ~, ef] = quadprog(H, f, [], [], [], [], U_lb_koop, U_ub_koop, [], opts);
    if ef <= 0 || any(isnan(Uopt))
        Uopt = -(H \ f);
        Uopt = min(max(Uopt, U_lb_koop), U_ub_koop);
    end

    u0 = min(max(Uopt(1:nu), u_min_koop), u_max_koop);
    u_koop(:,k) = u0;

    t_koop_mpc(k) = toc(tic_k);

    % Propagate plant
    x = x_true_koop;
    for j = 1:steps
        x = rk4_step(x, u0, params, dt);
        x(7:10) = normalize_quat(x(7:10));
    end
    x_true_koop = x;

    % Update lifted + integral (with anti-windup clamp)
    z   = observable_phi(x_true_koop, p_lift);
    ep  = x_true_koop(1:3) - p_ref_ctrl(:, ref_k);
    eta = eta + ep * ts;
    eta = max(min(eta, eta_max), -eta_max);   % anti-windup
    z_aug = [z; eta; 1];

    % Log
    p_koop(:,k)     = x_true_koop(1:3);
    v_koop(:,k)     = x_true_koop(4:6);
    q_koop(:,k)     = normalize_quat(x_true_koop(7:10));
    theta_koop(:,k) = vee(so3_log(QuatToRot(q_koop(:,k))));
    omega_koop(:,k) = x_true_koop(11:13);
end

fprintf('Koopman MPC: mean solve time = %.4f ms\n', 1000*mean(t_koop_mpc));

%% ====== 6b. Koopman precomputation timing (online-only) ===============
% Koopman is a linear predictor and linearises about the reference, so its
% per-instant matrices are precomputable, exactly as for SINDy. Online work
% is: lift the measured state, one matvec f = M*z_aug - f_c, and the QP.
fprintf('\n--- Koopman timing: precomputation analysis ---\n');
t0k = tic;
Hk_pre = cell(Nsim,1); Mk_pre = cell(Nsim,1); fck_pre = cell(Nsim,1);
for k = 1:Nsim
    ref_k = min(k, Nref_ctrl);
    Aseq = cell(Nh,1); Bseq = cell(Nh,1);
    for i = 1:Nh
        idx = min(ref_k+i-1, Nref_ctrl);
        uref_i = u_ref_high(:, steps*(idx-1)+1);
        [Ai,Bi] = localLin(z_ref_ctrl(:,idx), uref_i);
        pref = p_ref_ctrl(:, idx);
        Aseq{i} = [Ai, zeros(Nz,3), zeros(Nz,1); ts*C_p, eye(3), -ts*pref; zeros(1,Nz), zeros(1,3), 1];
        Bseq{i} = [Bi; zeros(3,nu); zeros(1,nu)];
    end
    Aqp = zeros(Nz_aug*Nh, Nz_aug); Bqp = zeros(Nz_aug*Nh, nu*Nh); Ap = eye(Nz_aug);
    for i = 1:Nh
        Ap = Aseq{i}*Ap; Aqp((i-1)*Nz_aug+1:i*Nz_aug,:) = Ap;
        for j = 1:i
            Aprod = eye(Nz_aug);
            for s = j+1:i, Aprod = Aseq{s}*Aprod; end
            Bqp((i-1)*Nz_aug+1:i*Nz_aug,(j-1)*nu+1:j*nu) = Aprod*Bseq{j};
        end
    end
    Yk = zeros(Nz_aug*Nh,1);
    for i = 1:Nh
        idx = min(ref_k+i-1, Nref_ctrl); base = (i-1)*Nz_aug;
        Yk(base+1:base+Nz) = z_ref_ctrl(:,idx);
        Yk(base+Nz+4)      = 1;
    end
    Hk = 2*(Bqp'*Qbar_koop*Bqp + Rbar_koop);
    Hk = 0.5*(Hk+Hk') + 1e-8*eye(size(Hk));
    Hk_pre{k}=Hk; Mk_pre{k}=2*(Bqp'*Qbar_koop*Aqp); fck_pre{k}=2*(Bqp'*Qbar_koop*Yk);
end
t_precomp_k = toc(t0k);

% online-only pass: lift + f + QP timed (plant propagation excluded)
x_kv = [zeros(6,1);1;0;0;0;zeros(3,1)]; x_kv(7:10)=normalize_quat(x_kv(7:10));
z_kv = observable_phi(x_kv, p_lift); eta_kv = zeros(3,1); zaug_kv = [z_kv; eta_kv; 1];
u_koop_v = zeros(nu,Nsim); t_online_k = zeros(1,Nsim);
for k = 1:Nsim
    ref_k = min(k,Nref_ctrl);
    t1 = tic;                                    % f + QP
    f = Mk_pre{k}*zaug_kv - fck_pre{k};
    [Uopt,~,ef] = quadprog(Hk_pre{k}, f, [], [], [], [], U_lb_koop, U_ub_koop, [], opts);
    if ef<=0 || any(isnan(Uopt)), Uopt = -(Hk_pre{k}\f); Uopt = min(max(Uopt,U_lb_koop),U_ub_koop); end
    u0 = min(max(Uopt(1:nu), u_min_koop), u_max_koop);
    ta = toc(t1);
    u_koop_v(:,k) = u0;
    x = x_kv;                                    % propagate (not timed)
    for j = 1:steps, x = rk4_step(x, u0, params, dt); x(7:10)=normalize_quat(x(7:10)); end
    x_kv = x;
    t2 = tic;                                    % lift + integral (online)
    z_kv   = observable_phi(x_kv, p_lift);
    eta_kv = max(min(eta_kv + (x_kv(1:3)-p_ref_ctrl(:,ref_k))*ts, eta_max), -eta_max);
    zaug_kv = [z_kv; eta_kv; 1];
    t_online_k(k) = ta + toc(t2);
end
du_max_k = max(abs(u_koop_v(:) - u_koop(:)));
fprintf('  verification  max|u_precomp - u_online| = %.2e  (should be ~0)\n', du_max_k);
fprintf('  Koopman online per-step        : %7.3f ms\n', 1000*mean(t_online_k));
fprintf('  offline precompute (total)     : %7.1f ms  (%.3f ms / instant)\n', ...
        1000*t_precomp_k, 1000*t_precomp_k/Nsim);

%% ======================== 7. SINDy LTV-MPC ============================
fprintf('\n--- Running SINDy LTV-MPC ---\n');

% Cost weights (on deviation states)
Q_pos_s = diag([300 200 1]);
Q_vel_s = diag([100 100 100]);
Q_q_s   = diag([1 1 1 1]);
Q_w_s   = diag([2500 2500 2500]);
Qx_sindy = blkdiag(Q_pos_s, Q_vel_s, Q_q_s, Q_w_s);
Qbar_sindy = kron(eye(Nh), Qx_sindy);

w_ft_s = 0.2; w_M_s = 4;
Ru_sindy  = diag([w_ft_s, w_M_s, w_M_s, w_M_s]);
Rbar_sindy = kron(eye(Nh), Ru_sindy);

% Init
x_true_sindy = [zeros(3,1); zeros(3,1); [1;0;0;0]; zeros(3,1)];
x_true_sindy(7:10) = normalize_quat(x_true_sindy(7:10));

p_sindy     = zeros(3, Nsim);
v_sindy     = zeros(3, Nsim);
theta_sindy = zeros(3, Nsim);
q_sindy     = zeros(4, Nsim);        % quaternion log (for SO(3) geodesic error)
omega_sindy = zeros(3, Nsim);
u_sindy     = zeros(nu, Nsim);
t_sindy_mpc = zeros(1, Nsim);
% split timers (before-optimisation breakdown; validates the Table III label)
t_jac_s = 0; t_cond_s = 0; t_qp_s = 0;

for k = 1:Nsim
    tic_k = tic;

    ref_k = min(k, Nref_ctrl);
    xref0 = x_ref_ctrl(:, ref_k);
    dx0   = x_true_sindy - xref0;
    dx0(7:10) = normalize_quat(dx0(7:10));

    Aseq_s = cell(Nh, 1);
    Bseq_s = cell(Nh, 1);
    DU_lb  = zeros(nu*Nh, 1);
    DU_ub  = zeros(nu*Nh, 1);

    tj = tic;                                    % --- Jacobian build
    for i = 1:Nh
        idx  = min(ref_k + i - 1, Nref_ctrl);
        xr_i = x_ref_ctrl(:, idx);
        ur_i = u_ref_ctrl(:, idx);

        [Ac, Bc] = sindy_jacobian_fd(xr_i, ur_i, Xi);
        Ad = eye(n) + ts*Ac;
        Bd = ts*Bc;

        Aseq_s{i} = Ad;
        Bseq_s{i} = Bd;

        base = (i-1)*nu;
        DU_lb(base+1:base+nu) = u_min - ur_i;
        DU_ub(base+1:base+nu) = u_max - ur_i;
    end
    t_jac_s = t_jac_s + toc(tj);

    tc = tic;                                    % --- condensation + Hessian
    Aqp_s = zeros(n*Nh, n);
    Bqp_s = zeros(n*Nh, nu*Nh);
    Ap = eye(n);
    for i = 1:Nh
        Ap = Aseq_s{i} * Ap;
        Aqp_s((i-1)*n+1:i*n, :) = Ap;
        for j = 1:i
            Aprod = eye(n);
            for s = j+1:i, Aprod = Aseq_s{s} * Aprod; end
            Bqp_s((i-1)*n+1:i*n, (j-1)*nu+1:j*nu) = Aprod * Bseq_s{j};
        end
    end

    Y_s = zeros(n*Nh, 1);   % deviation target = 0

    H = 2*(Bqp_s'*Qbar_sindy*Bqp_s + Rbar_sindy);
    f = 2*(Bqp_s'*Qbar_sindy*(Aqp_s*dx0 - Y_s));
    H = 0.5*(H + H') + 1e-9*eye(size(H));
    t_cond_s = t_cond_s + toc(tc);

    tq = tic;                                    % --- QP solve
    [DUopt, ~, ef] = quadprog(H, f, [], [], [], [], DU_lb, DU_ub, [], opts);
    if ef <= 0 || any(isnan(DUopt))
        DUopt = -(H \ f);
        DUopt = min(max(DUopt, DU_lb), DU_ub);
    end
    t_qp_s = t_qp_s + toc(tq);

    du0   = DUopt(1:nu);
    uref0 = u_ref_ctrl(:, ref_k);
    u0    = min(max(uref0 + du0, u_min), u_max);
    u_sindy(:,k) = u0;

    t_sindy_mpc(k) = toc(tic_k);

    % Propagate plant
    x = x_true_sindy;
    for j = 1:steps
        x = rk4_step(x, u0, params, dt);
        x(7:10) = normalize_quat(x(7:10));
    end
    x_true_sindy = x;

    % Log
    p_sindy(:,k)     = x_true_sindy(1:3);
    v_sindy(:,k)     = x_true_sindy(4:6);
    q_sindy(:,k)     = normalize_quat(x_true_sindy(7:10));
    theta_sindy(:,k) = vee(so3_log(QuatToRot(q_sindy(:,k))));
    omega_sindy(:,k) = x_true_sindy(11:13);
end

fprintf('SINDy MPC: mean solve time = %.4f ms\n', 1000*mean(t_sindy_mpc));

%% =============== 7b. Precomputation timing analysis ===================
% The LTV-MPC linearises about the REFERENCE (x_ref, u_ref), which is known
% ahead of time. Hence every per-instant matrix -- Ad,Bd, the condensed
% S_x,S_u, the Hessian H, the input bounds, and M = 2 S_u' Qbar S_x -- is
% independent of the state measurement and can be built offline. The only
% online work is one matvec f = M*dx0 and the box-constrained QP.
fprintf('\n--- SINDy timing: precomputation analysis ---\n');

% (a) offline: discrete Jacobians once at each reference point
t0 = tic;
Ad_ref = cell(Nref_ctrl,1); Bd_ref = cell(Nref_ctrl,1);
for j = 1:Nref_ctrl
    [Ac,Bc] = sindy_jacobian_fd(x_ref_ctrl(:,j), u_ref_ctrl(:,j), Xi);
    Ad_ref{j} = eye(n) + ts*Ac;
    Bd_ref{j} = ts*Bc;
end
% (b) offline: per-instant condensation, Hessian, bounds, and M matrix
Hpre  = cell(Nsim,1);
Mpre  = cell(Nsim,1);              % (nu*Nh x n): online gradient f = Mpre{k}*dx0
LBpre = zeros(nu*Nh, Nsim);
UBpre = zeros(nu*Nh, Nsim);
for k = 1:Nsim
    ref_k = min(k, Nref_ctrl);
    Aseq = cell(Nh,1); Bseq = cell(Nh,1);
    for i = 1:Nh
        idx = min(ref_k+i-1, Nref_ctrl);
        Aseq{i} = Ad_ref{idx};
        Bseq{i} = Bd_ref{idx};
        base = (i-1)*nu;
        LBpre(base+1:base+nu, k) = u_min - u_ref_ctrl(:, idx);
        UBpre(base+1:base+nu, k) = u_max - u_ref_ctrl(:, idx);
    end
    Aqp = zeros(n*Nh, n); Bqp = zeros(n*Nh, nu*Nh); Ap = eye(n);
    for i = 1:Nh
        Ap = Aseq{i}*Ap;
        Aqp((i-1)*n+1:i*n,:) = Ap;
        for j = 1:i
            Aprod = eye(n);
            for s = j+1:i, Aprod = Aseq{s}*Aprod; end
            Bqp((i-1)*n+1:i*n,(j-1)*nu+1:j*nu) = Aprod*Bseq{j};
        end
    end
    Hk = 2*(Bqp'*Qbar_sindy*Bqp + Rbar_sindy);
    Hk = 0.5*(Hk+Hk') + 1e-9*eye(size(Hk));
    Hpre{k} = Hk;
    Mpre{k} = 2*(Bqp'*Qbar_sindy*Aqp);
end
t_precomp = toc(t0);

% (c) online-only pass: re-simulate the plant, timing just f = M*dx0 + QP
x_true_v  = [zeros(6,1); 1;0;0;0; zeros(3,1)];
x_true_v(7:10) = normalize_quat(x_true_v(7:10));
u_sindy_v = zeros(nu, Nsim);
t_online  = zeros(1, Nsim);
for k = 1:Nsim
    ton   = tic;
    ref_k = min(k, Nref_ctrl);
    dx0   = x_true_v - x_ref_ctrl(:, ref_k);
    dx0(7:10) = normalize_quat(dx0(7:10));
    f = Mpre{k} * dx0;
    [DUopt, ~, ef] = quadprog(Hpre{k}, f, [], [], [], [], ...
                              LBpre(:,k), UBpre(:,k), [], opts);
    if ef <= 0 || any(isnan(DUopt))
        DUopt = -(Hpre{k} \ f);
        DUopt = min(max(DUopt, LBpre(:,k)), UBpre(:,k));
    end
    t_online(k) = toc(ton);
    u0 = min(max(u_ref_ctrl(:, ref_k) + DUopt(1:nu), u_min), u_max);
    u_sindy_v(:,k) = u0;
    x = x_true_v;
    for j = 1:steps, x = rk4_step(x, u0, params, dt); x(7:10)=normalize_quat(x(7:10)); end
    x_true_v = x;
end

% (d) verify the precomputed pass reproduces the official inputs
du_max = max(abs(u_sindy_v(:) - u_sindy(:)));

% (e) report (online-only; split retained in vars t_jac_s/t_cond_s/t_qp_s if needed)
fprintf('  verification  max|u_precomp - u_online| = %.2e  (should be ~0)\n', du_max);
fprintf('  SINDy online per-step          : %7.3f ms\n', 1000*mean(t_online));
fprintf('  offline precompute (total)     : %7.1f ms  (%.3f ms / instant)\n', ...
        1000*t_precomp, 1000*t_precomp/Nsim);

%% ======================== 8. COMPARISON PLOTS =========================
%  Two figures total (Koopman vs SINDy, both against the reference):
%    fig_states : 6 x 2 grid
%        Left  col : p_x, p_y, p_z, v_x, v_y, v_z
%        Right col : theta_x, theta_y, theta_z, omega_x, omega_y, omega_z
%        where theta = vee(log(R)) is the SO(3) log-map rotation vector.
%    fig_inputs : 4 x 1 stack (f_t, M_1, M_2, M_3)
%  ----------------------------------------------------------------------
fprintf('\n--- Generating comparison plots ---\n');

tvec = (0:Nsim-1) * ts;
resultsDir = fullfile(thisDir, 'results');

% --- Continuous (unwrapped) per-axis rotation vector ---------------------
%   The principal log vee(log R) wraps by 2*pi*n whenever the rotation angle
%   crosses pi. We instead accumulate the small relative-rotation logs
%   between consecutive samples, giving a jump-free axis-angle trace that
%   can be plotted component-by-component like p and v.
theta_ref_vec   = zeros(3, Nsim);
theta_sindy_vec = zeros(3, Nsim);
theta_koop_vec  = zeros(3, Nsim);
Rref_prev   = QuatToRot(normalize_quat(x_ref_ctrl(7:10, 1)));
Rsindy_prev = QuatToRot(normalize_quat(q_sindy(:, 1)));
Rkoop_prev  = QuatToRot(normalize_quat(q_koop(:, 1)));
theta_ref_vec(:, 1)   = vee(so3_log(Rref_prev));
theta_sindy_vec(:, 1) = vee(so3_log(Rsindy_prev));
theta_koop_vec(:, 1)  = vee(so3_log(Rkoop_prev));
for k = 2:Nsim
    R_ref_k   = QuatToRot(normalize_quat(x_ref_ctrl(7:10, k)));
    R_sindy_k = QuatToRot(normalize_quat(q_sindy(:, k)));
    R_koop_k  = QuatToRot(normalize_quat(q_koop(:, k)));
    theta_ref_vec(:, k)   = theta_ref_vec(:, k-1)   + vee(so3_log(R_ref_k   * Rref_prev.'));
    theta_sindy_vec(:, k) = theta_sindy_vec(:, k-1) + vee(so3_log(R_sindy_k * Rsindy_prev.'));
    theta_koop_vec(:, k)  = theta_koop_vec(:, k-1)  + vee(so3_log(R_koop_k  * Rkoop_prev.'));
    Rref_prev = R_ref_k;  Rsindy_prev = R_sindy_k;  Rkoop_prev = R_koop_k;
end

% Color / style conventions
c_ref   = [0 0 0];                  % black, dashed
c_sindy = [0    0.45 0.74];         % blue
c_koop  = [0.85 0.33 0.10];         % orange
lw_ref  = 1.2;
lw_data = 1.4;

% --- Figure 1: States grid (6 rows x 2 cols = 12 panels) -----------------
%   Left column (top->bot) : p_x p_y p_z v_x v_y v_z
%   Right column (top->bot): theta_x theta_y theta_z omega_x omega_y omega_z
fig_states = figure('Name','State Tracking Grid','Color',[1 1 1], ...
                    'Units','centimeters','Position',[2 2 9 18]);
tl = tiledlayout(fig_states, 6, 2, ...
                 'TileSpacing','tight','Padding','compact');

left_labels  = {'p_x (m)','p_y (m)','p_z (m)', ...
                'v_x (m/s)','v_y (m/s)','v_z (m/s)'};
right_labels = {'\theta_x (rad)','\theta_y (rad)','\theta_z (rad)', ...
                '\omega_x (rad/s)','\omega_y (rad/s)','\omega_z (rad/s)'};

left_ref   = [p_ref_ctrl(:,1:Nsim); v_ref_ctrl(:,1:Nsim)];
left_sindy = [p_sindy;              v_sindy];
left_koop  = [p_koop;               v_koop];
right_ref   = [theta_ref_vec;   omega_ref_ctrl(:,1:Nsim)];
right_sindy = [theta_sindy_vec; omega_sindy];
right_koop  = [theta_koop_vec;  omega_koop];

for r = 1:6
    % Left column
    ax = nexttile(tl, (r-1)*2 + 1);
    plot(ax, tvec, left_ref(r,:),   '--', 'Color', c_ref,   'LineWidth', lw_ref);  hold(ax,'on');
    plot(ax, tvec, left_sindy(r,:), '-',  'Color', c_sindy, 'LineWidth', lw_data);
    plot(ax, tvec, left_koop(r,:),  '-',  'Color', c_koop,  'LineWidth', lw_data);
    ylabel(ax, left_labels{r});
    grid(ax,'on'); box(ax,'on'); set(ax,'FontSize',8);
    if r < 6, set(ax,'XTickLabel',[]); end
    if r == 1
        legend(ax, {'Reference','SINDy-LTV-MPC','Koopman-MPC'}, ...
               'Location','best','FontSize',7);
    end

    % Right column
    ax = nexttile(tl, (r-1)*2 + 2);
    plot(ax, tvec, right_ref(r,:),   '--', 'Color', c_ref,   'LineWidth', lw_ref);  hold(ax,'on');
    plot(ax, tvec, right_sindy(r,:), '-',  'Color', c_sindy, 'LineWidth', lw_data);
    plot(ax, tvec, right_koop(r,:),  '-',  'Color', c_koop,  'LineWidth', lw_data);
    ylabel(ax, right_labels{r});
    grid(ax,'on'); box(ax,'on'); set(ax,'FontSize',8);
    if r < 6, set(ax,'XTickLabel',[]); end
end

xlabel(tl, 'Time (s)');
exportgraphics(fig_states, fullfile(resultsDir, 'states_grid.png'), ...
               'Resolution', 300);

% --- Figure 2: Control inputs (4 rows x 1 col) ---------------------------
fig_inputs = figure('Name','Control Inputs','Color',[1 1 1], ...
                    'Units','centimeters','Position',[2 2 9 14]);
ti = tiledlayout(fig_inputs, 4, 1, ...
                 'TileSpacing','tight','Padding','compact');
labels_u = {'f_t (N)','M_1 (N\cdotm)','M_2 (N\cdotm)','M_3 (N\cdotm)'};
for i = 1:4
    ax = nexttile(ti);
    plot(ax, tvec, u_sindy(i,:), '-', 'Color', c_sindy, 'LineWidth', lw_data); hold(ax,'on');
    plot(ax, tvec, u_koop(i,:),  '-', 'Color', c_koop,  'LineWidth', lw_data);
    ylabel(ax, labels_u{i});
    grid(ax,'on'); box(ax,'on');
    set(ax,'FontSize',9);
    if i == 1
        legend(ax, {'SINDy-LTV-MPC','Koopman-MPC'}, ...
               'Location','best','FontSize',8);
    end
    if i < 4
        set(ax,'XTickLabel',[]);
    end
end
xlabel(ti, 'Time (s)');
exportgraphics(fig_inputs, fullfile(resultsDir, 'inputs_grid.png'), ...
               'Resolution', 300);

%% ======================== 9. SUMMARY TABLE ============================
fprintf('\n');
fprintf('================================================================\n');
fprintf('                   COMPARISON SUMMARY\n');
fprintf('================================================================\n');

% Tracking RMSE (MPC closed-loop)
trk_koop_p = sqrt(mean(sum((p_koop - p_ref_ctrl(:,1:Nsim)).^2, 1)));
trk_koop_v = sqrt(mean(sum((v_koop - v_ref_ctrl(:,1:Nsim)).^2, 1)));
trk_koop_th = sqrt(mean(sum((theta_koop - theta_ref_ctrl(:,1:Nsim)).^2, 1)));
trk_koop_w = sqrt(mean(sum((omega_koop - omega_ref_ctrl(:,1:Nsim)).^2, 1)));

trk_sindy_p = sqrt(mean(sum((p_sindy - p_ref_ctrl(:,1:Nsim)).^2, 1)));
trk_sindy_v = sqrt(mean(sum((v_sindy - v_ref_ctrl(:,1:Nsim)).^2, 1)));
trk_sindy_th = sqrt(mean(sum((theta_sindy - theta_ref_ctrl(:,1:Nsim)).^2, 1)));
trk_sindy_w = sqrt(mean(sum((omega_sindy - omega_ref_ctrl(:,1:Nsim)).^2, 1)));

fprintf('\n--- Open-Loop Prediction nRMSE (%%) ---\n');
fprintf('%-12s %10s %10s\n', 'State', 'Koopman', 'SINDy');
fprintf('%-12s %9.2f%% %9.2f%%\n', 'Position',  koop_nRMSE_p,  sindy_nRMSE_p);
fprintf('%-12s %9.2f%% %9.2f%%\n', 'Velocity',  koop_nRMSE_v,  sindy_nRMSE_v);
fprintf('%-12s %9.2f%% %9.2f%%\n', 'Theta',     koop_nRMSE_th, sindy_nRMSE_th);
fprintf('%-12s %9.2f%% %9.2f%%\n', 'Omega',     koop_nRMSE_w,  sindy_nRMSE_w);

fprintf('\n--- Closed-Loop Tracking RMSE ---\n');
fprintf('%-12s %10s %10s\n', 'State', 'Koopman', 'SINDy');
fprintf('%-12s %10.4f %10.4f\n', 'Position',  trk_koop_p,  trk_sindy_p);
fprintf('%-12s %10.4f %10.4f\n', 'Velocity',  trk_koop_v,  trk_sindy_v);
fprintf('%-12s %10.4f %10.4f\n', 'Theta',     trk_koop_th, trk_sindy_th);
fprintf('%-12s %10.4f %10.4f\n', 'Omega',     trk_koop_w,  trk_sindy_w);

fprintf('\n--- MPC Computation Time ---\n');
fprintf('Koopman: mean=%.2f ms, max=%.2f ms\n', 1000*mean(t_koop_mpc), 1000*max(t_koop_mpc));
fprintf('SINDy:   mean=%.2f ms, max=%.2f ms\n', 1000*mean(t_sindy_mpc), 1000*max(t_sindy_mpc));
fprintf('================================================================\n');

%% ======================== 10. SAVE WORKSPACE ==========================
% Quaternion series (for continuous/unwrapped attitude plotting downstream)
q_ref = x_ref_ctrl(7:10, 1:Nsim);
save(fullfile(resultsDir, 'comparison_results.mat'), ...
    'koop_nRMSE_p','koop_nRMSE_v','koop_nRMSE_th','koop_nRMSE_w', ...
    'sindy_nRMSE_p','sindy_nRMSE_v','sindy_nRMSE_th','sindy_nRMSE_w', ...
    'trk_koop_p','trk_koop_v','trk_koop_th','trk_koop_w', ...
    'trk_sindy_p','trk_sindy_v','trk_sindy_th','trk_sindy_w', ...
    't_koop_mpc','t_sindy_mpc', ...
    'p_koop','v_koop','theta_koop','omega_koop','u_koop','q_koop', ...
    'p_sindy','v_sindy','theta_sindy','omega_sindy','u_sindy','q_sindy', ...
    'p_ref_ctrl','v_ref_ctrl','theta_ref_ctrl','omega_ref_ctrl','q_ref', ...
    'tvec','A_koop','B_koop','Xi','params');

fprintf('\nResults saved to %s\n', resultsDir);
fprintf('Simulation complete.\n');
