function wavestar_added_mass_delay_probe(coreDir, applicationsDir, repoDir, outDir, massFactor)
% Reproduce the first 9.95 s of the pinned WaveStar NMPC application.
% The alternate mass factor is a source sensitivity diagnostic only.
if nargin < 5
    massFactor = 2;
end
assert(ismember(massFactor, [1, 2]), 'Use the default or factor-one diagnostic');
assert(~isfolder(outDir), 'Use a fresh output directory');
mkdir(outDir);
source = fullfile(applicationsDir, 'WECCCOMP');
scratch = fullfile(outDir, 'WECCCOMP');
[copied, message] = copyfile(source, scratch);
assert(copied, message);
[copied, message] = copyfile(fullfile(repoDir, 'tests', 'fixtures', ...
    'matlab', 'wavestar_fixed_input_step', 'wavestar.h5'), ...
    fullfile(scratch, 'hydroData', 'wavestar.h5'));
assert(copied, message);
caseDir = fullfile(scratch, 'WECCCOMP_Nonlinear_Model_Predictive');
inputFile = fullfile(caseDir, 'wecSimInputFile.m');
input = fileread(inputFile);
assert(contains(input, 'simu.dt             = 50/1000;') && ...
    contains(input, "simu.solver         = 'ode8';") && ...
    contains(input, 'simu.cicEndTime     = 2;'), ...
    'The pinned WaveStar source settings changed');
for name = ["Float", "Arm", "Frame", "BC", "Motor"]
    oldName = name + '.stl';
    newName = name + '.STL';
    assert(contains(input, oldName) && ...
        isfile(fullfile(scratch, 'geometry', char(newName))), ...
        'The pinned geometry name changed');
    input = strrep(input, char(oldName), char(newName));
end
input = strrep(input, 'simu.cicEndTime     = 2;', ...
    sprintf('simu.endTime = 9.95;\n    simu.cicEndTime     = 2;'));
if massFactor == 1
    needle = "body(1) = bodyClass('../hydroData/wavestar.h5');";
    assert(contains(input, needle), 'The pinned body definition changed');
    input = strrep(input, needle, needle + newline + ...
        'body(1).adjMassFactor = 1;');
end
write_text(inputFile, input);
write_text(fullfile(caseDir, 'userDefinedFunctions.m'), ...
    '% Post-run plotting disabled for this batch audit.');

assert(~isempty(which('quadprog')) && ~isempty(which('dlqr')), ...
    'The published controller toolboxes are unavailable');
addpath(genpath(fullfile(coreDir, 'source')));
cd(caseDir);
wecSim;
assert(strcmp(simu.solver, 'ode8') && simu.dt == .05 && ...
    simu.endTime == 9.95 && numel(output.bodies(1).time) == 200 && ...
    body(1).adjMassFactor == massFactor, ...
    'The source probe did not run at the selected published settings');

record = output.bodies(1);
hydro = body(1).hydroForce.hf1;
time = record.time;
pitch = record.position(:, 5);
speed = record.velocity(:, 5);
acceleration = record.acceleration;
appliedAddedMass = hydro.storage.output_forceAddedMass;
modifiedAddedMass = hydro.storage.hydroForce_fAddedMass;
originalAddedMass = hydro.fAddedMass;
adjustedMass = hydro.adjustedMass;
rigidMass = body(1).mass;
aux4Position = output.bodies(4).position;
aux4Velocity = output.bodies(4).velocity;
aux4Acceleration = output.bodies(4).acceleration;
aux5Position = output.bodies(5).position;
aux5Velocity = output.bodies(5).velocity;
aux5Acceleration = output.bodies(5).acceleration;
assert(isequal(size(acceleration), [200, 6]) && ...
        isequal(size(appliedAddedMass), [200, 6]) && ...
        isequal(size(modifiedAddedMass), [6, 6]) && ...
        isequal(size(aux4Position), [200, 6]) && ...
        isequal(size(aux5Position), [200, 6]) && ...
        all(isfinite([pitch, speed, acceleration, appliedAddedMass]), 'all'), ...
    'The source acceleration or applied-force history is incomplete');
save(fullfile(outDir, 'added_mass_probe.mat'), 'time', 'pitch', ...
    'speed', 'acceleration', 'appliedAddedMass', ...
    'modifiedAddedMass', 'originalAddedMass', ...
    'adjustedMass', 'rigidMass', 'massFactor', ...
    'aux4Position', 'aux4Velocity', 'aux4Acceleration', ...
    'aux5Position', 'aux5Velocity', 'aux5Acceleration', '-v7');
fprintf('WAVESTAR_ADDED_MASS_PROBE factor=%d alpha0=%.12g mass=%.12g\n', ...
    massFactor, acceleration(1, 5), adjustedMass);
end

function write_text(path, value)
file = fopen(path, 'w');
assert(file ~= -1, 'Could not write disposable WaveStar file');
cleanup = onCleanup(@() fclose(file));
fprintf(file, '%s', value);
end
