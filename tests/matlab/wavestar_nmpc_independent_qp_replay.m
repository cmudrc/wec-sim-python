function wavestar_nmpc_independent_qp_replay(inputMat, applicationDir, outputCsv, horizonCsv)
% Replay pinned NMPC.m on independently generated Python controller inputs.
addpath(applicationDir);
data = load(inputMat);
n = numel(data.time);
assert(isequal(size(data.states), [5, n]));
assert(isequal(size(data.forecasts), [40, n]));
assert(isequal(size(data.python_commands), [n, 1]));
clear NMPC;
alpha = (1/0.7 + 0.7)/2;
beta = (1/0.7 - 0.7)/2;
commands = zeros(n, 1);
saveHorizons = nargin >= 4 && ~isempty(horizonCsv);
if saveHorizons
    horizons = zeros(n, 40);
end
for i = 1:n
    horizon = NMPC(data.states(:,i), data.forecasts(:,i), ...
        data.A, data.Bu, data.Bw, data.R, 40, alpha, beta, ...
        1000, [0 1; 1 0], 12, -12);
    commands(i) = horizon(1);
    if saveHorizons
        horizons(i,:) = horizon(:)';
    end
end
writematrix([data.time, commands, data.python_commands], outputCsv);
if saveHorizons
    writematrix([data.time, horizons], horizonCsv);
end
fprintf('max MATLAB/Python independent-input QP gap = %.9g N m\n', ...
    max(abs(commands - data.python_commands)));
end
