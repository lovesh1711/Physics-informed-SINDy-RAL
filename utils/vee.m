function a = vee(S)
% VEE  Inverse hat map: so(3) -> R^3
%   Extracts the 3-vector from a skew-symmetric matrix
    a = zeros(3,1);
    a(1) = -S(2,3);
    a(2) =  S(1,3);
    a(3) = -S(1,2);
end
