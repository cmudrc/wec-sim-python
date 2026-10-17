function rm3_two_column_motion_reference(sourceRoot, outputDir)
% Advance RM3 with a seeded two-column MAT spectrum and save native motion.
% sourceRoot is a checkout of pinned WEC-Sim; outputDir receives raw CSVs.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
sourceCase = fullfile(sourceRoot, 'examples', 'RM3');
fixture = fullfile(repoRoot, 'tests', 'fixtures', 'matlab', ...
    'imported_spectrum_two_column', 'spectrum.mat');
assert(isfile(fixture), 'The two-column spectrum fixture is missing');
addpath(genpath(fullfile(sourceRoot, 'source')));
if ~isfolder(outputDir)
    mkdir(outputDir);
end
caseDir = tempname;
mkdir(caseDir);
mkdir(fullfile(caseDir, 'geometry'));
mkdir(fullfile(caseDir, 'hydroData'));
for name = {'RM3.slx', 'userDefinedFunctions.m'}
    copyfile(fullfile(sourceCase, name{1}), fullfile(caseDir, name{1}));
end
for name = {'float.stl', 'plate.stl'}
    copyfile(fullfile(sourceCase, 'geometry', name{1}), ...
        fullfile(caseDir, 'geometry', name{1}));
end
copyfile(fullfile(sourceCase, 'hydroData', 'rm3.h5'), ...
    fullfile(caseDir, 'hydroData', 'rm3.h5'));

contents = fileread(fullfile(sourceCase, 'wecSimInputFile.m'));
contents = replaceOnce(contents, 'simu.explorer = ''on'';', ...
    'simu.explorer = ''off'';');
contents = replaceOnce(contents, 'simu.rampTime = 100;', ...
    'simu.rampTime = 1;');
contents = replaceOnce(contents, 'simu.endTime = 400;', ...
    'simu.endTime = 30;');
contents = replaceOnce(contents, 'waves = waveClass(''regular'');', ...
    'waves = waveClass(''spectrumImport'');');
contents = replaceOnce(contents, 'waves.height = 2.5;', ...
    sprintf('waves.spectrumFile = ''%s'';\nwaves.phaseSeed = 1;', fixture));
contents = replaceOnce(contents, 'waves.period = 8;', '');
fid = fopen(fullfile(caseDir, 'wecSimInputFile.m'), 'w');
assert(fid > 0, 'Cannot write temporary RM3 input');
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, contents);
clear cleanup;

cd(caseDir);
wecSim;
assert(simu.dt == 0.1 && simu.endTime == 30 && ...
    simu.rampTime == 1 && strcmp(simu.solver, 'ode4') && ...
    strcmp(waves.type, 'spectrumImport') && waves.phaseSeed == 1 && ...
    numel(output.bodies(1).time) == 301 && numel(waves.omega) == 5, ...
    'Unexpected source settings or output shape');
writematrix([output.wave.time(:), output.wave.elevation(:)], ...
    fullfile(outputDir, 'wave.csv'));
for iBody = 1:2
    record = output.bodies(iBody);
    writematrix([record.time(:), record.position, record.velocity], ...
        fullfile(outputDir, sprintf('body%d.csv', iBody)));
end
record = output.ptos(1);
writematrix([record.time(:), record.position, record.velocity, ...
    record.forceInternalMechanics, record.powerInternalMechanics], ...
    fullfile(outputDir, 'pto.csv'));
writematrix([waves.omega(:), waves.phase(:)], ...
    fullfile(outputDir, 'components.csv'));
close_system('RM3', 0);
end

function result = replaceOnce(contents, before, after)
starts = regexp(contents, ['(?m)^' regexptranslate('escape', before)], ...
    'start');
assert(numel(starts) == 1, ...
    'The pinned RM3 input changed: %s', before);
first = starts(1);
result = [contents(1:first-1), after, ...
    contents(first+numel(before):end)];
end
