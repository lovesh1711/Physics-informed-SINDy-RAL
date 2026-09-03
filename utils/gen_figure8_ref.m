function [Xref, Uref, tvec] = gen_figure8_ref(params, ts, T_total, A, w, z0)
% GEN_FIGURE8_REF  Differentially-flat figure-8 (lemniscate) reference.
%   Dynamically-consistent reference {x_r,u_r} for quad_dynamics.m
%   (NED, x=[p;v;q;omega], u=[ft;M]).  Flat outputs: figure-8 position in the
%   North-East plane at height z0, and yaw held at 0.
%
%   From p,v,a the thrust magnitude and body-z axis follow algebraically; the
%   body frame is completed from yaw; body rates omega are the exact inverse of
%   quat_kinematics (omega = 2 W(q)' q_dot), and moments M come from Euler's
%   equation.  All derivatives are evaluated on a fine grid and downsampled to
%   the control rate ts, so the reference is feasible to finite-difference
%   accuracy O(dt^2).
%
%   Inputs:  params(m,J,g), ts, T_total (s), A (amplitude m), w (rad/s),
%            z0 (height coord, NED).   Outputs: Xref(13xN), Uref(4xN), tvec.

    m = params.m; J = params.J; g = params.g;
    e3 = [0;0;1];  B = A/2;
    dt = 1e-3;                          % fine grid for accurate derivatives
    Nf = round(T_total/dt) + 1;
    tf = (0:Nf-1)*dt;
    psi = 0; xc = [cos(psi); sin(psi); 0];

    % --- position/velocity/attitude/thrust on the fine grid ---
    pf = zeros(3,Nf); vf = zeros(3,Nf); qf = zeros(4,Nf); ftf = zeros(1,Nf);
    for i = 1:Nf
        t = tf(i);
        p = [ A*sin(w*t);        B*sin(2*w*t);       z0 ];
        v = [ A*w*cos(w*t);      2*B*w*cos(2*w*t);    0 ];
        a = [-A*w^2*sin(w*t);   -4*B*w^2*sin(2*w*t);  0 ];
        c = g*e3 - a;  ftf(i) = m*norm(c);  b3 = c/norm(c);
        b2 = cross(b3,xc); b2 = b2/norm(b2);  b1 = cross(b2,b3);
        R_wb = [b1, b2, b3];               % body-to-world
        qk = normalize_quat( RotToQuat(R_wb') );   % QuatToRot(q)=R_bw=R_wb'
        if i > 1 && dot(qk, qf(:,i-1)) < 0, qk = -qk; end   % continuous q
        pf(:,i) = p;  vf(:,i) = v;  qf(:,i) = qk;
    end

    % --- omega = 2 W(q)' q_dot  (exact inverse of quat_kinematics) ---
    omf = zeros(3,Nf);
    for i = 1:Nf
        ip = min(i+1,Nf); im = max(i-1,1); dtk = (ip-im)*dt;
        qdot = (qf(:,ip) - qf(:,im))/dtk;
        qw = qf(1,i); qv = qf(2:4,i);
        W = [-qv'; qw*eye(3) + hat(qv)];
        omf(:,i) = 2*(W'*qdot);
    end

    % --- moments M = J omega_dot + omega x J omega ---
    Mf = zeros(3,Nf);
    for i = 1:Nf
        ip = min(i+1,Nf); im = max(i-1,1); dtk = (ip-im)*dt;
        omdot = (omf(:,ip) - omf(:,im))/dtk;
        Mf(:,i) = J*omdot + cross(omf(:,i), J*omf(:,i));
    end

    % --- downsample to control rate ---
    stride = round(ts/dt);
    idx  = 1:stride:Nf;
    tvec = tf(idx);
    Xref = [pf(:,idx); vf(:,idx); qf(:,idx); omf(:,idx)];
    Uref = [ftf(idx); Mf(:,idx)];
end
