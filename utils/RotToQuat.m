function q = RotToQuat(R)
% ROTTOQUAT  Converts a rotation matrix to a unit quaternion [qw;qx;qy;qz]
%   Written by Daniel Mellinger
%   Reference: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

    tr = R(1,1) + R(2,2) + R(3,3);

    if (tr > 0)
        S = sqrt(tr+1.0) * 2;
        qw = 0.25 * S;
        qx = (R(3,2) - R(2,3)) / S;
        qy = (R(1,3) - R(3,1)) / S;
        qz = (R(2,1) - R(1,2)) / S;
    elseif ((R(1,1) > R(2,2)) && (R(1,1) > R(3,3)))
        S = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3)) * 2;
        qw = (R(3,2) - R(2,3)) / S;
        qx = 0.25 * S;
        qy = (R(1,2) + R(2,1)) / S;
        qz = (R(1,3) + R(3,1)) / S;
    elseif (R(2,2) > R(3,3))
        S = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3)) * 2;
        qw = (R(1,3) - R(3,1)) / S;
        qx = (R(1,2) + R(2,1)) / S;
        qy = 0.25 * S;
        qz = (R(2,3) + R(3,2)) / S;
    else
        S = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2)) * 2;
        qw = (R(2,1) - R(1,2)) / S;
        qx = (R(1,3) + R(3,1)) / S;
        qy = (R(2,3) + R(3,2)) / S;
        qz = 0.25 * S;
    end

    q = [qw; qx; qy; qz];
    q = q * sign(qw);
end
