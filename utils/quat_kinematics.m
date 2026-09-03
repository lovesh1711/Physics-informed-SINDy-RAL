function q_dot = quat_kinematics(q, omega)
% QUAT_KINEMATICS  Quaternion time derivative from angular velocity
%   q_dot = 0.5 * Omega(omega) * q
    wx = omega(1); wy = omega(2); wz = omega(3);

    Omega = [  0   -wx  -wy  -wz;
              wx    0    wz  -wy;
              wy  -wz    0    wx;
              wz   wy  -wx    0 ];

    q_dot = 0.5 * Omega * q;
end
