function oswec_hinge_body_damping_reference(caseDir, sourceDir, outputDir, fixedStep)
% Run a derived Nonhydro_Body input with flap body damping, not PTO damping.
if nargin < 4
    fixedStep = 0.1;
end
assert(ismember(fixedStep, [0.1, 0.025]), 'Unexpected solver-only step');
addpath(genpath(sourceDir));
priorDir = pwd;
cd(caseDir);
run('initializeWecSim');
assert(simu.dt == 0.1 && simu.endTime == 150 && ...
    strcmp(waves.type, 'regular') && waves.height == 2.5 && ...
    waves.period == 8 && pto(1).damping == 0, ...
    'Derived input changed more than the declared time and damping');
assert(isequal(size(body(1).linearDamping), [6, 6]), ...
    'Flap damping must be a 6x6 matrix');
assert(all(isfinite(body(1).linearDamping), 'all'), ...
    'Flap damping must be finite');
[~, modelName, ~] = fileparts(simu.simMechanicsFile);
assert(abs(str2double(get_param(modelName, 'FixedStep')) - .1) < 1e-12, ...
    'Unexpected initialized fixed step');
set_param(modelName, 'FixedStep', num2str(fixedStep, 16));
variables = who;
for k = 1:numel(variables)
    assignin('base', variables{k}, eval(variables{k}));
end
sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
assert(abs(str2double(get_param(modelName, 'FixedStep')) - fixedStep) < 1e-12, ...
    'The source solver step changed during simulation');
run('stopWecSim');
assert(numel(output.bodies(1).time) == 1501 && ...
    numel(output.bodies(2).time) == 1501 && ...
    numel(output.ptos(1).time) == 1501, 'Incomplete source run');
if ~isfolder(outputDir)
    mkdir(outputDir);
end
for k = 1:2
    record = output.bodies(k);
    writematrix([record.time(:), record.position, record.velocity], ...
        fullfile(outputDir, sprintf('body%d.csv', k)));
end
record = output.ptos(1);
writematrix([record.time(:), record.position, record.velocity, ...
    record.forceTotal], fullfile(outputDir, 'pto.csv'));
writematrix(body(1).linearDamping, fullfile(outputDir, 'damping.csv'));
writematrix(fixedStep, fullfile(outputDir, 'fixed_step.csv'));
close_system(modelName, 0);
cd(priorDir);
end
