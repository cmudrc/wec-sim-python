function sphere_hydrostatic_stiffness_reference(caseDir, sourceDir, outputCsv)
% Run the pinned Sphere 1 m free decay with a derived heave restoring term.
addpath(genpath(sourceDir));
priorDir = pwd;
cd(caseDir);
run('initializeWecSim');
assert(simu.dt == .01 && simu.endTime == 40 && ...
    strcmp(waves.type, 'noWaveCIC'), 'Published case settings changed');
assert(isequal(body(1).hydroStiffness, ...
    diag([0 0 1155583.665 0 0 0])), ...
    'Expected only the derived hydrostatic stiffness setting');
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
