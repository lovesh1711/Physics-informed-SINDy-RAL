function x_next = rk4_step_sindy(x, u, Xi, dt, steps)
% RK4_STEP_SINDY  Integrate SINDy model over multiple sub-steps using RK4
%   Integrates xdot = Xi * phi(x,u) for 'steps' substeps of size dt

    for j = 1:steps
        k1 = sindy_rhs(x, u, Xi);
        k2 = sindy_rhs(x + 0.5*dt*k1, u, Xi);
        k3 = sindy_rhs(x + 0.5*dt*k2, u, Xi);
        k4 = sindy_rhs(x + dt*k3, u, Xi);
        x  = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        x(7:10) = normalize_quat(x(7:10));
    end
    x_next = x;
end
