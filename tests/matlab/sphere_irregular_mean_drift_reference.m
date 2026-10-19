function sphere_irregular_mean_drift_reference(caseSource, coreSource, hydroFile, outDir, driftOption)
% Derive a short PM-wave input from the pinned Sphere Mean_Drift application.
if nargin < 5
    driftOption = 1;
end
assert(ismember(driftOption, [1, 2, 3]), ...
    'Expected a supported HDF5 mean-drift coefficient');
addpath(genpath(coreSource));
if ~isfolder(outDir)
    mkdir(outDir);
end
caseDir = fullfile(outDir, 'Mean_Drift_irregular');
assert(~isfolder(caseDir), 'The derived case directory already exists');
[copied, message] = copyfile(caseSource, caseDir);
assert(copied, message);
if ~isfolder(fullfile(caseDir, 'hydroData'))
    mkdir(fullfile(caseDir, 'hydroData'));
end
copyfile(hydroFile, fullfile(caseDir, 'hydroData', 'sphere.h5'));
% The application's optional plotting script assumes regular-wave scalar
% amplitude; it is unrelated to the simulated force or body response.
if isfile(fullfile(caseDir, 'userDefinedFunctions.m'))
    delete(fullfile(caseDir, 'userDefinedFunctions.m'));
end

inputFile = fullfile(caseDir, 'wecSimInputFile.m');
input = fileread(inputFile);
edits = {
    'simu.rampTime           = 20;', 'simu.rampTime           = 10;';
    'simu.endTime            = 100;', 'simu.endTime            = 30;';
    'simu.explorer           = ''on'';', 'simu.explorer           = ''off'';';
    'waves = waveClass(''regularCIC'');', [ ...
        'waves = waveClass(''irregular'');' newline ...
        'waves.spectrumType = ''PM'';' newline ...
        'waves.phaseSeed = 1;'];
    'waves.period = 2;', 'waves.period = 8;';
    'waves.height = 0.1;', 'waves.height = 0.5;';
};
for index = 1:size(edits, 1)
    assert(contains(input, edits{index, 1}), ...
        'Pinned Sphere Mean_Drift input changed');
    input = strrep(input, edits{index, 1}, edits{index, 2});
end
if driftOption ~= 1
    oldDrift = 'body(1).meanDrift = 1;';
    assert(contains(input, oldDrift), 'Pinned Sphere drift option changed');
    input = strrep(input, oldDrift, ...
        sprintf('body(1).meanDrift = %d;', driftOption));
end
geometry = fullfile(fileparts(caseSource), '_Common_Input_Files', ...
    'Sphere', 'geometry', 'sphere.stl');
oldGeometry = '''../_Common_Input_Files/Sphere/geometry/sphere.stl''';
assert(contains(input, oldGeometry), 'Pinned Sphere geometry path changed');
input = strrep(input, oldGeometry, ['''' geometry '''']);
fid = fopen(inputFile, 'w');
assert(fid ~= -1, 'Could not configure Sphere irregular drift input');
fprintf(fid, '%s', input);
fclose(fid);

priorDir = pwd;
cd(caseDir);
wecSim;
assert(simu.dt == 0.01 && simu.endTime == 30 && ...
    simu.rampTime == 10 && simu.cicEndTime == 10 && ...
    strcmp(waves.type, 'irregular') && ...
    strcmp(waves.spectrumType, 'PM') && waves.phaseSeed == 1 && ...
    waves.height == 0.5 && waves.period == 8 && ...
    body(1).meanDrift == driftOption, 'Derived Sphere settings changed');
record = output.bodies(1);
assert(numel(record.time) == 3001, 'Incomplete Sphere source run');
assert(all(isfinite(record.forceExcitation), 'all'), ...
    'Nonfinite source excitation');
writematrix([record.time(:), record.position, record.velocity, ...
    record.forceExcitation], fullfile(outDir, 'body.csv'));
writematrix([output.wave.time(:), output.wave.elevation(:)], ...
    fullfile(outDir, 'wave.csv'));
components = [waves.omega(:), waves.amplitude(:), ...
    waves.dOmega(:), waves.phase(:)];
assert(size(components, 2) == 4 && ...
    all(isfinite(components), 'all'), 'Invalid source sea components');
writematrix(components, fullfile(outDir, 'components.csv'));
writematrix([body(1).mass, body(1).inertia], ...
    fullfile(outDir, 'mass.csv'));
close_system('sphere', 0);
cd(priorDir);
end
