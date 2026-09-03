function print_sindy_model(Xi, feature_names, state_names, tol)
% PRINT_SINDY_MODEL  Pretty-print the identified SINDy sparse model
    n = size(Xi, 1);
    fprintf('\n================ SINDy identified model ================\n');

    for i = 1:n
        coeffs = Xi(i, :);
        nz = find(abs(coeffs) > tol);

        fprintf('\n%s_dot =\n', state_names{i});

        if isempty(nz)
            fprintf('  0\n');
            continue;
        end

        [~, ord] = sort(abs(coeffs(nz)), 'descend');
        nz = nz(ord);

        for k = 1:length(nz)
            j = nz(k);
            fprintf('  %+ .6e * %s\n', coeffs(j), feature_names{j});
        end
    end

    fprintf('\n========================================================\n');
end
