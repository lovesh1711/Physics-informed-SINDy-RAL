function [Theta, feature_names] = sindy_library(Xdata, Udata)
% SINDY_LIBRARY  Physics-informed feature library for quadrotor SINDy
%   Features: {1, p, v, q, w, u, R(:), ft*b3w, q*w, w_i*w_j}  (48 total)
%   Input:  Xdata (13 x K) state snapshots, Udata (4 x K) input snapshots
%   Output: Theta (Nfeat x K) feature matrix, feature_names (cell array)

    K  = size(Xdata, 2);
    e3 = [0;0;1];

    Theta = [];
    feature_names = {};

    for k = 1:K
        x = Xdata(:, k);
        u = Udata(:, k);

        p  = x(1:3);
        v  = x(4:6);
        q  = x(7:10) / norm(x(7:10));   % unit norm only (sign-consistent with x)
        w  = x(11:13);
        ft = u(1);

        R_bw = QuatToRot(q);
        R_wb = R_bw';
        b3w  = R_wb * e3;

        qomega = [
            q(1)*w(1); q(1)*w(2); q(1)*w(3);
            q(2)*w(1); q(2)*w(2); q(2)*w(3);
            q(3)*w(1); q(3)*w(2); q(3)*w(3);
            q(4)*w(1); q(4)*w(2); q(4)*w(3)
        ];

        % Gyroscopic products w_i*w_j (the omega-omega coupling in w_dot)
        wprod = [
            w(1)*w(1); w(2)*w(2); w(3)*w(3);
            w(1)*w(2); w(1)*w(3); w(2)*w(3)
        ];

        phi = [
            1;
            p; v; q; w;
            u;
            R_wb(:);
            ft*b3w;
            qomega;
            wprod
        ];

        if isempty(Theta)
            feature_names = make_feature_names();
            Theta = zeros(length(feature_names), K);
        end

        Theta(:, k) = phi;
    end
end
