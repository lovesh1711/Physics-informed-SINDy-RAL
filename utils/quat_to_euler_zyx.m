function eul = quat_to_euler_zyx(Q)
% QUAT_TO_EULER_ZYX  ZYX (yaw-pitch-roll) Euler angles of a quaternion series.
%   Q   : 4 x N quaternions [qw;qx;qy;qz]; R(q) = QuatToRot(q) is world-to-body
%   eul : 3 x N [roll; pitch; yaw] in rad, defined by the body-to-world matrix
%         R(q)' = Rz(yaw) * Ry(pitch) * Rx(roll)
%   Roll and yaw are unwrapped along the series, so full turns appear as
%   continuous curves instead of jumps of 2*pi. The angles are singular at
%   pitch = +-pi/2 (gimbal lock).
    N   = size(Q, 2);
    eul = zeros(3, N);
    for k = 1:N
        Rwb = QuatToRot(normalize_quat(Q(:, k))).';      % body-to-world
        eul(:, k) = [atan2(Rwb(3,2), Rwb(3,3));            % roll
                     -asin(max(-1, min(1, Rwb(3,1))));     % pitch
                     atan2(Rwb(2,1), Rwb(1,1))];           % yaw
    end
    eul([1 3], :) = unwrap(eul([1 3], :), [], 2);
end
