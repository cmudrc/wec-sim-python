function wavestar_auxiliary_mass_probe(coreDir, applicationsDir, repoDir, outDir)
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
        'The pinned WaveStar geometry name changed');
    input = strrep(input, char(oldName), char(newName));
end
input = strrep(input, 'simu.cicEndTime     = 2;', ...
    sprintf('simu.endTime = 9.95;\n    simu.cicEndTime     = 2;'));
old = 'body(4).mass            = 0.0001;';
new = 'body(4).mass            = 1e-9;';
assert(contains(input, old));
input = strrep(input, old, new);
old = 'body(4).inertia         = [0.0001 0.0001 0.0001];';
new = 'body(4).inertia         = [1e-9 1e-9 1e-9];';
assert(contains(input, old));
input = strrep(input, old, new);
old = 'body(5).mass            = 0.0001;';
new = 'body(5).mass            = 1e-9;';
assert(contains(input, old));
input = strrep(input, old, new);
old = 'body(5).inertia         = [0.0001 0.0001 0.0001];';
new = 'body(5).inertia         = [1e-9 1e-9 1e-9];';
assert(contains(input, old));
input = strrep(input, old, new);
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
    body(4).mass == 1e-9 && body(5).mass == 1e-9, ...
    'The source probe did not run at the selected settings');
record = output.bodies(1);
time = record.time;
pitch = record.position(:, 5);
speed = record.velocity(:, 5);
acceleration = record.acceleration;
save(fullfile(outDir, 'aux_mass_probe.mat'), 'time', 'pitch', ...
    'speed', 'acceleration', '-v7');
fprintf('AUX_MASS_PROBE_PASSED samples=%d alpha0=%.12g\n', ...
    numel(time), acceleration(1, 5));
end

function write_text(path, value)
file = fopen(path, 'w');
assert(file ~= -1, 'Could not write disposable WaveStar file');
cleanup = onCleanup(@() fclose(file));
fprintf(file, '%s', value);
end
