function q = normalize_quat(q)
% NORMALIZE_QUAT  Normalize quaternion to unit norm with canonical sign (qw >= 0)
    q = q / max(norm(q), 1e-12);
    if q(1) < 0
        q = -q;
    end
end
