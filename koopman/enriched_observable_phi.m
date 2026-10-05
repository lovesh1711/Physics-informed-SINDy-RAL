function psi = enriched_observable_phi(x, p)
% ENRICHED_OBSERVABLE_PHI  Koopman observables with explicit R(:).
%
%   Difference vs. paper/koopman/observable_phi.m:
%     - Adds R(:)          (9 entries)  — matches the existing convention:
%                                         R = QuatToRot(q) (the matrix used in
%                                         the original `R*hat(w)^i` terms).
%     - The thrust direction b3 = R'*e3 is NOT added separately: it is the
%       third row of R, whose entries are already in R(:). (An earlier
%       version appended b3, giving 3 duplicate observables, Nz = 52.)
%
%   psi = [ p(3); v(3); q(4); omega(3);             % 13
%           R(:)(9);                                % 9   (NEW)
%           vec(R*hat(w)^1); ...; vec(R*hat(w)^p) ] % 9*p
%
%   For p=3, Nz = 13 + 9 + 27 = 49.
%
%   C = [eye(13), zeros(13, Nz-13)] still extracts the original state.

    p_lin = x(1:3);
    v     = x(4:6);
    q     = normalize_quat(x(7:10));
    omega = x(11:13);

    R  = QuatToRot(q);       % same matrix used in paper/koopman/observable_phi.m
    w_hat = hat(omega);
    h_stack = [];
    A_pow = w_hat;
    for i = 1:p
        Hi = R * A_pow;
        h_stack = [h_stack; Hi(:)]; %#ok<AGROW>
        A_pow = A_pow * w_hat;
    end

    psi = [p_lin; v; q; omega; R(:); h_stack];
end
