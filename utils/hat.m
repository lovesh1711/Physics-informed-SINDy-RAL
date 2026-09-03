function S = hat(w)
% HAT  Skew-symmetric (hat) map: R^3 -> so(3)
%   hat(w) returns the 3x3 skew-symmetric matrix such that hat(a)*b = cross(a,b)
    S = [  0   -w(3)  w(2);
          w(3)   0   -w(1);
         -w(2)  w(1)   0  ];
end
