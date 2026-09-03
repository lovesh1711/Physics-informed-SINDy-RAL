function Theta = so3_log(R)
% SO3_LOG  Matrix logarithm on SO(3)
%   Returns Theta in so(3) such that expm(Theta) ≈ R

    % reproject to SO(3) via SVD
    [U,~,V] = svd(R);
    R = U*V';

    cos_phi = (trace(R) - 1)/2;
    cos_phi = max(min(cos_phi,1),-1);
    phi = acos(cos_phi);

    if abs(phi) < 1e-6
        Theta = 0.5*(R - R');
    else
        Theta = (phi/(2*sin(phi))) * (R - R');
    end
end
