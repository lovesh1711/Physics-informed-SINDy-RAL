function psi = observable_phi(x, p)
% OBSERVABLE_PHI  Koopman physics-informed observable (quaternion-based)
%   Builds lifted state: [p; v; q; omega; R*what^1(:); R*what^2(:); ...; R*what^p(:)]
%   Input:  x = [p(3); v(3); q(4); omega(3)]  (13-dim)
%           p = number of hat-power basis functions
%   Output: psi in R^{13 + 9*p}

    p_lin = x(1:3);
    v     = x(4:6);
    q     = normalize_quat(x(7:10));
    omega = x(11:13);

    R     = QuatToRot(q);
    w_hat = hat(omega);

    h_stack = [];
    A = w_hat;
    for i = 1:p
        Hi = R * A;
        h_stack = [h_stack; Hi(:)];
        A = A * w_hat;
    end

    psi = [p_lin; v; q; omega; h_stack];
end
