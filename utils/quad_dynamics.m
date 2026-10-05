function x_dot = quad_dynamics(~, x, u, params)
% QUAD_DYNAMICS  Continuous-time quadrotor dynamics with quaternion orientation
%   State: x = [p(3); v(3); q(4); omega(3)]  (13-dim)
%   Input: u = [ft; M1; M2; M3]
%
%   Convention:
%     R_bw = QuatToRot(q)   world-to-body rotation
%     R_wb = R_bw'           body-to-world rotation
%     Thrust along body +z axis, mapped to world frame

    m = params.m;
    J = params.J;
    g = params.g;
    e3 = [0;0;1];

    v     = x(4:6);
    q     = x(7:10) / norm(x(7:10));   % unit norm only: no sign flip inside RK4 stages
    omega = x(11:13);

    ft = u(1);
    M  = u(2:4);

    R_bw = QuatToRot(q);      % world -> body
    R_wb = R_bw';             % body -> world

    p_dot     = v;
    Fw        = R_wb * (ft * e3);
    v_dot     = g*e3 - (1/m)*Fw;
    q_dot     = quat_kinematics(q, omega);
    omega_dot = J \ (M - cross(omega, J*omega));

    x_dot = [p_dot; v_dot; q_dot; omega_dot];
end
