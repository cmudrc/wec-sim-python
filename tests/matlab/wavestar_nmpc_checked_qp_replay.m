function wavestar_nmpc_checked_qp_replay(inputMat, applicationDir, outputMat)
% Replay the pinned published QP on the checked compatible source's inputs.
addpath(applicationDir);
data = load(inputMat);
n = numel(data.time);
assert(isequal(size(data.states), [5, n]));
assert(isequal(size(data.forecasts), [40, n]));
assert(isequal(size(data.source_commands), [n, 1]));
assert(n == 4201);
clear NMPC;
alpha = (1/0.7 + 0.7)/2;
beta = (1/0.7 - 0.7)/2;
commands = zeros(n, 1);
plans = zeros(n, 40);
for i = 1:n
    horizon = NMPC(data.states(:,i), data.forecasts(:,i), ...
        data.A, data.Bu, data.Bw, data.R, 40, alpha, beta, ...
        1000, [0 1; 1 0], 12, -12);
    commands(i) = horizon(1);
    plans(i,:) = horizon(:)';
end
time = data.time;
source_commands = data.source_commands;
save(outputMat, 'time', 'commands', 'plans', 'source_commands', '-v7');
fprintf('max source/direct-MATLAB request gap = %.9g N m\n', ...
    max(abs(commands - source_commands)));
end
