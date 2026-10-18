function wavestar_nmpc_independent_qp_replay(inputMat, applicationDir, outputCsv)
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
for i = 1:n
    horizon = NMPC(data.states(:,i), data.forecasts(:,i), ...
        data.A, data.Bu, data.Bw, data.R, 40, alpha, beta, ...
        1000, [0 1; 1 0], 12, -12);
    commands(i) = horizon(1);
end
writematrix([data.time, commands, data.python_commands], outputCsv);
fprintf('max MATLAB/Python independent-input QP gap = %.9g N m\n', ...
    max(abs(commands - data.python_commands)));
end
