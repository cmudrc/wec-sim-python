function rm3_moordyn_full_solver_step
% Full published sea with only ode45 MaxStep refined independently.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
addpath(fullfile(repoRoot, 'tests', 'matlab'));
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end
sourceDir = fullfile(repoRoot, 'applications', ...
    'Paraview_Visualization', 'RM3_MoorDyn_Viz');
caseDir = fullfile(fileparts(sourceDir), ...
    'RM3_MoorDyn_Viz_full_solver_step');
outDir = fullfile(repoRoot, 'matlab-rm3-full-solver-step');
mkdir(caseDir);
mkdir(outDir);
for filename = {'wecSimInputFile.m', 'RM3MoorDyn.slx', ...
        'userDefinedFunctions.m'}
    copyfile(fullfile(sourceDir, filename{1}), ...
        fullfile(caseDir, filename{1}));
end
mkdir(fullfile(caseDir, 'Mooring'));
copyfile(fullfile(sourceDir, 'Mooring', 'lines.txt'), ...
    fullfile(caseDir, 'Mooring', 'lines.txt'));
inputFile = fullfile(caseDir, 'wecSimInputFile.m');
contents = fileread(inputFile);
changes = {
    'simu.dtOut = 0.1;', 'simu.dtOut = 0.01;';
    'simu.paraview.option = 1;', 'simu.paraview.option = 0;';
    'waves.period = 8;', sprintf('waves.period = 8;\nwaves.phaseSeed = 1;')
};
for i = 1:size(changes, 1)
    assert(contains(contents, changes{i, 1}), ...
        'The pinned RM3 visualization input changed');
    contents = strrep(contents, changes{i, 1}, changes{i, 2});
end
assert(contains(contents, 'simu.dt = 0.01;') && ...
    contains(contents, 'simu.endTime = 80;'), ...
    'The pinned step or duration changed');
fid = fopen(inputFile, 'w');
assert(fid > 0);
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, contents);
clear cleanup;
cd(caseDir);
% Simulink's SrcWorkspace='parent' requires the script call context used by
% pinned wecSim.m. Keep simu.dt at 0.01 s for the MoorDyn trigger and
% source wave/radiation grids; override only the model's ode45 MaxStep.
rm3_moordyn_full_solver_simulate;
assert(size(waves.phase, 1) == 1000, 'Unexpected phase grid');
writematrix(waves.phase(:), fullfile(outDir, 'phase.csv'));
writematrix([output.wave.time(:), output.wave.elevation(:)], ...
    fullfile(outDir, 'wave.csv'));
for iBody = 1:2
    record = output.bodies(iBody);
    values = [record.time(:), record.position, record.velocity, ...
        record.forceTotal, record.forceExcitation];
    assert(isequal(size(values), [8001, 25]) && ...
        all(isfinite(values), 'all'), 'Incomplete body trace');
    writematrix(values, fullfile(outDir, ...
        sprintf('body%d.csv', iBody)));
end
record = output.mooring(1);
values = [record.time(:), record.position, ...
    record.velocity, record.forceMooring];
assert(isequal(size(values), [8001, 19]) && ...
    all(isfinite(values), 'all'), 'Incomplete MoorDyn trace');
writematrix(values, fullfile(outDir, 'mooring.csv'));
record = output.ptos(1);
values = [record.time(:), record.position, record.velocity, ...
    record.forceInternalMechanics, record.powerInternalMechanics];
assert(isequal(size(values), [8001, 25]) && ...
    all(isfinite(values), 'all'), 'Incomplete PTO trace');
writematrix(values, fullfile(outDir, 'pto.csv'));
close_system('RM3MoorDyn', 0);
end
