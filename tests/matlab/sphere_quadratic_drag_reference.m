function sphere_quadratic_drag_reference(caseDir, sourceDir, motionCsv, forceCsv, expectedDrag, expectedCd, expectedArea)
% Run the pinned Sphere 1 m free decay with either source body-drag input.
if nargin < 5
    expectedDrag = diag([0 0 250000 0 0 0]);
end
if nargin < 6
    expectedCd = zeros(1, 6);
    expectedArea = zeros(1, 6);
end
% Keep unrelated example post-processors off the path. stopWecSim invokes
% userDefinedFunctions.m whenever one is visible.
addpath(genpath(fullfile(sourceDir, 'source')));
priorDir = pwd;
cd(caseDir);
run('initializeWecSim');
assert(simu.dt == .01 && simu.endTime == 40 && ...
    strcmp(waves.type, 'noWaveCIC'), 'Published case settings changed');
assert(isequal(body(1).quadDrag.drag, expectedDrag), ...
    'Quadratic body drag differs from the expected derived setting');
assert(isequal(body(1).quadDrag.cd, expectedCd) && ...
    isequal(body(1).quadDrag.area, expectedArea), ...
    'Quadratic body coefficient or area differs from the expected setting');
assert(isequal(body(1).linearDamping, zeros(6)), ...
    'This reference has no extra linear body damping');
[~, modelName, ~] = fileparts(simu.simMechanicsFile);
variables = who;
for k = 1:numel(variables)
    assignin('base', variables{k}, eval(variables{k}));
end
sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
run('stopWecSim');
record = output.bodies(1);
assert(numel(record.time) == 4001, 'Incomplete source run');
writematrix([record.time(:), record.position, record.velocity], motionCsv);
writematrix([record.time(:), record.forceMorisonAndViscous], forceCsv);
close_system(modelName, 0);
cd(priorDir);
end
