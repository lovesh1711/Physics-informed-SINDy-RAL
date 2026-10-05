function names = make_feature_names()
% MAKE_FEATURE_NAMES  Returns ordered cell array of SINDy library feature names
%   Must match the order of features in sindy_library.m exactly

    names = {};

    names{end+1} = '1';

    for i = 1:3, names{end+1} = sprintf('p%d', i); end
    for i = 1:3, names{end+1} = sprintf('v%d', i); end

    names{end+1} = 'qw'; names{end+1} = 'qx';
    names{end+1} = 'qy'; names{end+1} = 'qz';
    for i = 1:3, names{end+1} = sprintf('w%d', i); end

    names{end+1} = 'ft'; names{end+1} = 'M1';
    names{end+1} = 'M2'; names{end+1} = 'M3';

    names{end+1} = 'R11'; names{end+1} = 'R21'; names{end+1} = 'R31';
    names{end+1} = 'R12'; names{end+1} = 'R22'; names{end+1} = 'R32';
    names{end+1} = 'R13'; names{end+1} = 'R23'; names{end+1} = 'R33';

    names{end+1} = 'ft*b3w1'; names{end+1} = 'ft*b3w2'; names{end+1} = 'ft*b3w3';

    qnames = {'qw','qx','qy','qz'};
    wnames = {'w1','w2','w3'};
    for qi = 1:4
        for wi = 1:3
            names{end+1} = sprintf('%s*%s', qnames{qi}, wnames{wi});
        end
    end

    % Gyroscopic products w_i*w_j (must match wprod order in sindy_library.m)
    names{end+1} = 'w1*w1'; names{end+1} = 'w2*w2'; names{end+1} = 'w3*w3';
    names{end+1} = 'w1*w2'; names{end+1} = 'w1*w3'; names{end+1} = 'w2*w3';

    names = names(:);
end
