function [Ak, Bk] = local_koopman_linearization(z_ref, u_ref, A, B)
% LOCAL_KOOPMAN_LINEARIZATION  Finite-difference Jacobians of Koopman predictor
%   Linearizes z_{k+1} = A*z_k + B*u_k around (z_ref, u_ref)
%   Returns local Jacobians Ak (dz+/dz), Bk (dz+/du)

    eps_fd = 1e-4;
    Nz = size(A, 1);
    nu = size(B, 2);

    z_base = A*z_ref + B*u_ref;

    Ak = zeros(Nz, Nz);
    for i = 1:Nz
        ei = zeros(Nz, 1); ei(i) = eps_fd;
        z_plus = A*(z_ref + ei) + B*u_ref;
        Ak(:, i) = (z_plus - z_base) / eps_fd;
    end

    Bk = zeros(Nz, nu);
    for j = 1:nu
        ej = zeros(nu, 1); ej(j) = eps_fd;
        z_plus = A*z_ref + B*(u_ref + ej);
        Bk(:, j) = (z_plus - z_base) / eps_fd;
    end
end
