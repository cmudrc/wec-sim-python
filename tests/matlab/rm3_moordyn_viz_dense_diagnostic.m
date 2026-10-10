function rm3_moordyn_viz_dense_diagnostic
% Save every MoorDyn coupling step from a short derived source run.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
caseDir = fullfile(repoRoot, 'applications', 'Paraview_Visualization', ...
    'RM3_MoorDyn_Viz');
denseDir = fullfile(fileparts(caseDir), 'RM3_MoorDyn_Viz_dense');
copyfile(caseDir, denseDir);
inputFile = fullfile(denseDir, 'wecSimInputFile.m');
contents = fileread(inputFile);
changes = {
    'simu.endTime = 80;', 'simu.endTime = 10;';
    'simu.dtOut = 0.1;', 'simu.dtOut = 0.01;';
    'simu.paraview.option = 1;', 'simu.paraview.option = 0;';
    'waves.period = 8;', sprintf('waves.period = 8;\nwaves.phaseSeed = 1;')
};
for i = 1:size(changes, 1)
    assert(contains(contents, changes{i, 1}), ...
        'The pinned RM3 MoorDyn visualization input changed');
    contents = strrep(contents, changes{i, 1}, changes{i, 2});
end
fid = fopen(inputFile, 'w');
assert(fid > 0);
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, contents);
clear cleanup;
cd(denseDir);
wecSim;
assert(simu.dt == 0.01 && simu.dtOut == 0.01 && ...
    simu.endTime == 10 && strcmp(simu.solver, 'ode45') && ...
    waves.phaseSeed == 1 && mooring(1).moorDyn == 1, ...
    'The dense coupling diagnostic settings changed');
record = output.mooring(1);
values = [record.time(:), record.position, record.velocity, ...
    record.forceMooring];
assert(size(values, 1) == 1001 && size(values, 2) == 19 && ...
    all(isfinite(values), 'all'), ...
    'The dense source mooring record is incomplete');
outDir = fullfile(repoRoot, 'matlab-rm3-moordyn-viz');
writematrix(values, fullfile(outDir, 'dense_mooring.csv'));
writematrix(waves.phase(:), fullfile(outDir, 'dense_phase.csv'));
writematrix([output.wave.time(:), output.wave.elevation(:)], ...
    fullfile(outDir, 'dense_wave.csv'));
for iBody = 1:2
    bodyRecord = output.bodies(iBody);
    bodyValues = [bodyRecord.time(:), bodyRecord.position, ...
        bodyRecord.velocity, bodyRecord.forceTotal, ...
        bodyRecord.forceExcitation];
    assert(size(bodyValues, 1) == 1001 && ...
        size(bodyValues, 2) == 25 && all(isfinite(bodyValues), 'all'), ...
        'The dense source body record is incomplete');
    writematrix(bodyValues, fullfile(outDir, ...
        sprintf('dense_body%d.csv', iBody)));
end
end
