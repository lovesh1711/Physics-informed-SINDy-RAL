%% ========================================================================
%  run_diagnostics.m  --  Tables I and II of the paper (Section III-E)
%  ------------------------------------------------------------------------
%  Generates the four training datasets D1-D4 and fits EDMD models to test
%  the predictions of Section III:
%    Table I  : predicted vs fitted beta_z = B*(v_z, f_t) for each
%               dataset / dictionary pair (eq. betaz)
%    Table II : one-step velocity error on D4 for three dictionaries of
%               increasing size, against the irreducible floor (Cor. 1)
%
%  Datasets (thrust = mg + zero-mean Gaussian perturbation in all four):
%    D1  near hover, held     : 100 x 0.5 s from hover, Sigma = 10 I,
%                               each input held 5 sample periods
%    D2  broad, held          : 600 x 1 s, tilt up to 45 deg, Sigma = 15 I,
%                               inputs held 3 or 7 sample periods
%    D3  near hover, i.i.d.   : 100 x 0.5 s from hover,
%                               Sigma = diag(25,1,1,1), new input every step
%    D4  broad, i.i.d.        : 600 x 1 s, tilt up to 90 deg,
%                               Sigma = diag(25,40,40,40), new input every step
%
%  Dictionaries: basic (rigid-body lifting, Nz = 40), enriched (basic + R(:),
%  Nz = 49), and '+const' (one constant observable appended).
%
%  Output: printed tables, results/diagnostics_results.mat
% ========================================================================
clc; clear; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'utils'),   '-begin');
addpath(fullfile(thisDir, 'koopman'), '-begin');

params.m = 4.34;
params.J = diag([0.0820; 0.0845; 0.1377]);
params.g = 9.81;
u_min = [10;  -30; -30; -30];
u_max = [300;  30;  30;  30];
ts    = 0.01;  dt = 0.001;  steps = round(ts/dt);
mg    = params.m * params.g;
lam   = 1e-6;                 % EDMD ridge
e3    = [0; 0; 1];
i_vz  = 6;                    % v_z row of the lifted state
j_ft  = 1;                    % thrust column of B

%% ---------- datasets D1-D4 ---------------------------------------------
cfg = struct('seed',1,'Ntraj',100,'T',0.5,'hold',5,'Sigma',diag([10;10;10;10]), ...
             'ic','hover','max_ang',0,'v0',0,'w0',0,'p0',0);
D.D1 = gen_data(cfg, params, ts, dt, steps, mg, u_min, u_max);
cfg = struct('seed',11,'Ntraj',600,'T',1.0,'hold',-1,'Sigma',diag([15;15;15;15]), ...
             'ic','broad','max_ang',pi/4,'v0',2.0,'w0',3.0,'p0',1.0);
D.D2 = gen_data(cfg, params, ts, dt, steps, mg, u_min, u_max);
cfg = struct('seed',21,'Ntraj',100,'T',0.5,'hold',1,'Sigma',diag([25;1;1;1]), ...
             'ic','hover','max_ang',0,'v0',0,'w0',0,'p0',0);
D.D3 = gen_data(cfg, params, ts, dt, steps, mg, u_min, u_max);
cfg = struct('seed',31,'Ntraj',600,'T',1.0,'hold',1,'Sigma',diag([25;40;40;40]), ...
             'ic','broad','max_ang',pi/2,'v0',2.0,'w0',6.0,'p0',1.0);
D.D4 = gen_data(cfg, params, ts, dt, steps, mg, u_min, u_max);

% body-z axis R(q)' e3 at every snapshot: mean and spread per dataset
fn = fieldnames(D);
for i = 1:numel(fn)
    d = D.(fn{i});
    b3 = zeros(3, d.K);
    for k = 1:d.K
        b3(:,k) = QuatToRot(normalize_quat(d.X(7:10,k)))' * e3;
    end
    D.(fn{i}).b3      = b3;
    D.(fn{i}).b3_mean = mean(b3, 2);
end

%% ---------- Table I: predicted vs fitted beta_z -------------------------
fits = { ...
  'D1', 'basic'      , @(x) observable_phi(x,3)               , 40; ...
  'D2', 'basic'      , @(x) observable_phi(x,3)               , 40; ...
  'D2', 'enriched'   , @(x) enriched_observable_phi(x,3)      , 49; ...
  'D2', 'basic+const', @(x) [observable_phi(x,3); 1]          , 41; ...
  'D3', 'basic+const', @(x) [observable_phi(x,3); 1]          , 41; ...
  'D4', 'basic+const', @(x) [observable_phi(x,3); 1]          , 41  };

fprintf('\n====================== Table I ======================\n');
fprintf('beta_z = B*(v_z, f_t), units 1e-3 s/kg;  true upright gain -Ts/m = %.2f\n', ...
        -1e3*ts/params.m);
fprintf('%-4s %-12s %4s %10s %10s %10s\n', 'Data', 'Dictionary', 'Nz', ...
        'E[(R''e3)_3]', 'predicted', 'fitted');
tab1 = struct([]);
for i = 1:size(fits,1)
    d = D.(fits{i,1});
    Bfull = fit_edmd(d, fits{i,3}, fits{i,4}, lam);
    r.dataset = fits{i,1};  r.dict = fits{i,2};  r.Nz = fits{i,4};
    r.mu_b3z  = d.b3_mean(3);
    r.pred    = -(ts/params.m) * d.b3_mean(3);           % eq. (betaz)
    r.fitted  = Bfull(i_vz, j_ft);
    tab1 = [tab1, r]; %#ok<AGROW>
    fprintf('%-4s %-12s %4d %+10.2f %+10.3f %+10.3f\n', r.dataset, r.dict, ...
            r.Nz, r.mu_b3z, 1e3*r.pred, 1e3*r.fitted);
end

%% ---------- Table II: error floor on D4 ---------------------------------
d  = D.D4;
Ut = d.U - mean(d.U, 2);                                   % input fluctuation
floor_emp = (ts/params.m)^2 * mean( Ut(1,:).^2 .* sum((d.b3 - d.b3_mean).^2, 1) );

t2 = { '[x;1]'          , @(x) [x; 1]                            , 14; ...
       'basic + const'  , @(x) [observable_phi(x,3); 1]          , 41; ...
       'enriched + const', @(x) [enriched_observable_phi(x,3); 1], 50 };

fprintf('\n====================== Table II =====================\n');
fprintf('one-step velocity MSE on D4, units 1e-4 m^2/s^2;  floor = %.3f\n', 1e4*floor_emp);
fprintf('%-17s %4s %8s %11s %13s\n', 'Dictionary', 'Nz', 'Total', 'Input part', 'Input/floor');
tab2 = struct([]);
for i = 1:size(t2,1)
    [Bfull, mse_v] = fit_edmd(d, t2{i,2}, t2{i,3}, lam);
    % input part on the velocity rows: || (Ts G(x) - C B)_v u ||^2
    err = Bfull(4:6,:)*Ut + (ts/params.m) * (d.b3 .* Ut(1,:));
    r2.dict = t2{i,1};  r2.Nz = t2{i,3};
    r2.total = mse_v;  r2.input = mean(sum(err.^2, 1));  r2.ratio = r2.input/floor_emp;
    tab2 = [tab2, r2]; %#ok<AGROW>
    fprintf('%-17s %4d %8.2f %11.3f %13.3f\n', r2.dict, r2.Nz, 1e4*r2.total, ...
            1e4*r2.input, r2.ratio);
end

resultsDir = fullfile(thisDir, 'results');
save(fullfile(resultsDir, 'diagnostics_results.mat'), 'tab1', 'tab2', 'floor_emp');
fprintf('\nsaved results/diagnostics_results.mat\n');

%% ===================== local functions ==================================
function d = gen_data(cfg, params, ts, dt, steps, mg, u_min, u_max)
    % cfg.hold: >0 = new input every cfg.hold steps; -1 = alternate 3 / 7
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',cfg.seed));
    n = 13; nu = 4;
    Nsteps = round(cfg.T/ts);  Ktot = cfg.Ntraj * Nsteps;
    X = zeros(n,Ktot); Xp = zeros(n,Ktot); U = zeros(nu,Ktot);
    kk = 0;
    for tr = 1:cfg.Ntraj
        if strcmp(cfg.ic,'hover')
            x = [zeros(6,1); 1;0;0;0; zeros(3,1)];
        else
            p0 = cfg.p0*(2*rand(3,1)-1);
            v0 = cfg.v0*(2*rand(3,1)-1);
            ax = randn(3,1); ax = ax/max(norm(ax),1e-9);
            ang = cfg.max_ang*rand;
            q0 = normalize_quat([cos(ang/2); sin(ang/2)*ax]);
            w0 = cfg.w0*(2*rand(3,1)-1);
            x = [p0; v0; q0; w0];
        end
        x(7:10) = normalize_quat(x(7:10));
        if cfg.hold == -1
            if mod(tr,2)==0, hold_len = 3; else, hold_len = 7; end
        else
            hold_len = cfg.hold;
        end
        u = [mg;0;0;0];
        for k = 1:Nsteps
            if (k==1) || (mod(k-1,hold_len)==0)
                u = mvnrnd(zeros(4,1), cfg.Sigma).';
                u(1) = u(1) + mg;
                u = min(max(u,u_min),u_max);
            end
            x_k = x;
            for j = 1:steps
                x = rk4_step(x, u, params, dt);
                x(7:10) = normalize_quat(x(7:10));
            end
            kk = kk+1;
            X(:,kk) = x_k; Xp(:,kk) = x; U(:,kk) = u;
        end
    end
    d.X = X; d.Xp = Xp; d.U = U; d.K = kk;
end

function [Bfull, mse_v] = fit_edmd(d, f, Nz, lam)
    % least-squares EDMD fit (eq. edmd_lsq) with a small ridge
    K = d.K;
    Z = zeros(Nz,K); Zp = zeros(Nz,K);
    for k = 1:K
        Z(:,k)  = f(nq(d.X(:,k)));
        Zp(:,k) = f(nq(d.Xp(:,k)));
    end
    ZU = [Z; d.U];
    G  = Zp*ZU' / (ZU*ZU' + lam*eye(size(ZU,1)));
    Bfull = G(:, Nz+1:end);
    P = G*ZU;
    mse_v = mean(sum((Zp(4:6,:) - P(4:6,:)).^2, 1));   % one-step velocity MSE
end

function x = nq(x)
    x(7:10) = normalize_quat(x(7:10));
end
