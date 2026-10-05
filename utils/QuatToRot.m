function bRw = QuatToRot(q)
% QUATTOROT  Converts a unit quaternion to rotation matrix (world-to-body)
%   q = [qw; qx; qy; qz],  bRw maps world frame to body frame.
%
%   The quaternion follows the Hamilton body-rate kinematics of
%   quat_kinematics.m (q_dot = 0.5*Omega(omega)*q), under which q rotates
%   body vectors into the world frame:  wRb = I + 2*qw*[qv]x + 2*[qv]x^2.
%   The world-to-body matrix is its transpose, which is eq. (R_from_q) of
%   the paper; the body z-axis in the world frame is bRw' * e3.
%   (Adapted from D. Mellinger's QuatToRot.)

    q = q ./ sqrt(sum(q.^2));

    qahat(1,2) = -q(4);
    qahat(1,3) =  q(3);
    qahat(2,3) = -q(2);
    qahat(2,1) =  q(4);
    qahat(3,1) = -q(3);
    qahat(3,2) =  q(2);

    wRb = eye(3) + 2*qahat*qahat + 2*q(1)*qahat;   % body-to-world
    bRw = wRb.';                                    % world-to-body
end
