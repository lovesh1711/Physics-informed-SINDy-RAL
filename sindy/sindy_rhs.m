function xdot = sindy_rhs(x, u, Xi)
% SINDY_RHS  Evaluate learned SINDy continuous-time model: xdot = Xi * phi(x,u)
    [Theta1, ~] = sindy_library(x, u);
    phi = Theta1(:, 1);
    xdot = Xi * phi;
end
