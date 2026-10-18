function wavestar_auxiliary_mass_probe(coreDir, applicationsDir, repoDir, outDir, duration)
if nargin < 5
    duration = 9.95;
end
assert(ismember(duration, [9.95, 20]), 'Expected a bounded prefix');
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
    sprintf('simu.endTime = %.17g;\n    simu.cicEndTime     = 2;', duration));
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
    simu.endTime == duration && ...
    numel(output.bodies(1).time) == round(duration/.05) + 1 && ...
    body(4).mass == 1e-9 && body(5).mass == 1e-9, ...
    'The source probe did not run at the selected settings');
record = output.bodies(1);
time = record.time;
pitch = record.position(:, 5);
speed = record.velocity(:, 5);
acceleration = record.acceleration;
commandTime = cmd_ptoM.time;
command = cmd_ptoM.signals.values;
components = [waves.omega(:), waves.amplitude(:), ...
    waves.dOmega(:), waves.phase(:)];
wave = [output.wave.time(:), output.wave.elevation(:)];
assert(isequal(size(commandTime), size(time)) && ...
    isequal(size(command), size(time)) && ...
    isequal(size(components), [500, 4]) && ...
    isequal(size(wave), [numel(time), 2]) && ...
    all(isfinite([time, pitch, speed, commandTime, command]), 'all'), ...
    'The controller or wave record is incomplete');
save(fullfile(outDir, 'aux_mass_probe.mat'), 'time', 'pitch', ...
    'speed', 'acceleration', 'commandTime', 'command', ...
    'components', 'wave', '-v7');
fprintf('AUX_MASS_PROBE_PASSED samples=%d alpha0=%.12g\n', ...
    numel(time), acceleration(1, 5));
end

function write_text(path, value)
file = fopen(path, 'w');
assert(file ~= -1, 'Could not write disposable WaveStar file');
cleanup = onCleanup(@() fclose(file));
fprintf(file, '%s', value);
end
