function wavestar_ode8_modal_probe(modesFile, outputFile)
% Measure one published-size Simulink ode8 step on every active WaveStar
% radiation pole. The 2x2 real blocks represent complex conjugate modes;
% this isolates solver amplification from the coupled plant and controller.
load(modesFile, 'eigenvalues');
eigenvalues = eigenvalues(:);
assert(numel(eigenvalues) == 92 && all(isfinite(eigenvalues)), ...
    'Expected the pinned 92 active WaveStar radiation poles');

count = numel(eigenvalues);
matrixA = zeros(2 * count);
initialState = zeros(2 * count, 1);
for index = 1:count
    pair = 2 * index - 1;
    realPart = real(eigenvalues(index));
    imaginaryPart = imag(eigenvalues(index));
    matrixA(pair:pair+1, pair:pair+1) = ...
        [realPart, -imaginaryPart; imaginaryPart, realPart];
    initialState(pair) = 1;
end

model = 'wecsimWaveStarOde8ModalProbe';
assert(~bdIsLoaded(model), 'Close the prior modal probe before rerunning');
new_system(model);
cleanup = onCleanup(@() close_system(model, 0));
workspace = get_param(model, 'ModelWorkspace');
assignin(workspace, 'matrixA', matrixA);
assignin(workspace, 'matrixB', zeros(2 * count, 1));
assignin(workspace, 'matrixC', eye(2 * count));
assignin(workspace, 'matrixD', zeros(2 * count, 1));
assignin(workspace, 'initialState', initialState);
add_block('simulink/Sources/Constant', [model '/zero'], 'Value', '0');
add_block('simulink/Continuous/State-Space', [model '/state'], ...
    'A', 'matrixA', 'B', 'matrixB', ...
    'C', 'matrixC', 'D', 'matrixD', ...
    'X0', 'initialState');
add_block('simulink/Sinks/To Workspace', [model '/output'], ...
    'VariableName', 'modalOutput', 'SaveFormat', 'Array');
add_line(model, 'zero/1', 'state/1');
add_line(model, 'state/1', 'output/1');
set_param(model, 'SolverType', 'Fixed-step', 'Solver', 'ode8', ...
    'FixedStep', '0.05', 'StopTime', '0.05');
result = sim(model);
values = result.get('modalOutput');
final = values(end, :);
amplification = final(1:2:end).' + 1i * final(2:2:end).';
assert(numel(amplification) == count && all(isfinite(amplification)), ...
    'The modal probe did not finish its single step');
save(outputFile, 'eigenvalues', 'amplification', '-v7');
fprintf('WAVESTAR_ODE8_MODAL_PROBE modes=%d stable=%d max_amplitude=%.12g\n', ...
    count, sum(abs(amplification) < 1), max(abs(amplification)));
end
