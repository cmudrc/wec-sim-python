function sphere_linear_damping_reference(caseDir, sourceDir, outputCsv)
% Run the pinned Sphere 1 m free decay with only body heave damping added.
addpath(genpath(sourceDir));
priorDir = pwd;
cd(caseDir);
run('initializeWecSim');
assert(simu.dt == .01 && simu.endTime == 40 && ...
    strcmp(waves.type, 'noWaveCIC'), 'Published case settings changed');
assert(isequal(body(1).linearDamping, diag([0 0 100000 0 0 0])), ...
    'Expected only the derived heave-damping setting');
[~, modelName, ~] = fileparts(simu.simMechanicsFile);
variables = who;
for k = 1:numel(variables)
    assignin('base', variables{k}, eval(variables{k}));
end
sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
run('stopWecSim');
record = output.bodies(1);
assert(numel(record.time) == 4001, 'Incomplete source run');
writematrix([record.time(:), record.position, record.velocity], outputCsv);
close_system(modelName, 0);
cd(priorDir);
end
