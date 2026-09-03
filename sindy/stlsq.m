function Xi = stlsq(Theta, Xdot, lambda, n_iter)
% STLSQ  Sequentially Thresholded Least Squares for SINDy
%   Theta: Nfeat x K  (feature library)
%   Xdot:  n x K      (state derivatives)
%   lambda: sparsity threshold
%   n_iter: number of thresholding iterations
%   Returns Xi: n x Nfeat  (sparse coefficient matrix)

    Nfeat = size(Theta, 1);
    n     = size(Xdot, 1);
    ridge = 1e-8;

    % initial least-squares solve
    Xi = (Theta*Theta' + ridge*eye(Nfeat)) \ (Theta*Xdot');
    Xi = Xi';                                                % n x Nfeat

    for it = 1:n_iter
        small = abs(Xi) < lambda;
        Xi(small) = 0;

        for i = 1:n
            big = find(~small(i,:));
            if isempty(big), continue; end
            Th = Theta(big, :);
            Xi(i, big) = (Th*Th' + ridge*eye(length(big))) \ (Th*Xdot(i,:)');
        end
    end
end
