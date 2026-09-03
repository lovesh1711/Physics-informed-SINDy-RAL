function [Ac, Bc] = sindy_jacobian_fd(xr, ur, Xi)
% SINDY_JACOBIAN_FD  Finite-difference Jacobians of continuous-time SINDy model
%   Computes Ac = df/dx, Bc = df/du at (xr, ur) using central differences

    n  = length(xr);
    nu = length(ur);

    Ac = zeros(n, n);
    Bc = zeros(n, nu);

    epsx = 1e-6;
    for i = 1:n
        dx = zeros(n, 1); dx(i) = epsx;
        xp = xr + dx;  xm = xr - dx;
        xp(7:10) = normalize_quat(xp(7:10));
        xm(7:10) = normalize_quat(xm(7:10));
        fp = sindy_rhs(xp, ur, Xi);
        fm = sindy_rhs(xm, ur, Xi);
        Ac(:, i) = (fp - fm) / (2*epsx);
    end

    epsu = 1e-6;
    for j = 1:nu
        du = zeros(nu, 1); du(j) = epsu;
        fp = sindy_rhs(xr, ur + du, Xi);
        fm = sindy_rhs(xr, ur - du, Xi);
        Bc(:, j) = (fp - fm) / (2*epsu);
    end
end
