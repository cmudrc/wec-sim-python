function rm3_moordyn_viz_full_fine
% Repeat the seeded published 80 s case at two finer coupling steps.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end
caseDir = fullfile(repoRoot, 'applications', 'Paraview_Visualization', ...
    'RM3_MoorDyn_Viz');
outDir = fullfile(repoRoot, 'matlab-rm3-moordyn-viz-full-fine');
mkdir(outDir);
steps = [0.000625, 0.0003125];
labels = {'fullfine_000625', 'fullfine_0003125'};

for iStep = 1:numel(steps)
    step = steps(iStep);
    label = labels{iStep};
    targetDir = fullfile(fileparts(caseDir), ['RM3_MoorDyn_Viz_' label]);
    mkdir(targetDir);
    for name = {'wecSimInputFile.m', 'RM3MoorDyn.slx', ...
            'userDefinedFunctions.m'}
        copyfile(fullfile(caseDir, name{1}), fullfile(targetDir, name{1}));
    end
    mkdir(fullfile(targetDir, 'Mooring'));
    copyfile(fullfile(caseDir, 'Mooring', 'lines.txt'), ...
        fullfile(targetDir, 'Mooring', 'lines.txt'));
    inputFile = fullfile(targetDir, 'wecSimInputFile.m');
    contents = fileread(inputFile);
    changes = {
        'simu.dt = 0.01;', sprintf('simu.dt = %.7f;', step);
        'simu.dtOut = 0.1;', 'simu.dtOut = 0.01;';
        'simu.paraview.option = 1;', 'simu.paraview.option = 0;';
        'waves.period = 8;', ...
            sprintf('waves.period = 8;\nwaves.phaseSeed = 1;')
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

    cd(targetDir);
    wecSim;
    assert(abs(simu.dt - step) < 1e-12 && ...
        simu.dtOut == 0.01 && simu.endTime == 80 && ...
        simu.rampTime == 0 && simu.paraview.option == 0 && ...
        strcmp(simu.solver, 'ode45') && waves.phaseSeed == 1 && ...
        mooring(1).moorDyn == 1 && size(waves.phase, 1) == 1000, ...
        'The full-duration refined source settings changed');
    writematrix(waves.phase(:), fullfile(outDir, [label '_phase.csv']));
    writematrix([output.wave.time(:), output.wave.elevation(:)], ...
        fullfile(outDir, [label '_wave.csv']));
    for iBody = 1:2
        record = output.bodies(iBody);
        bodyValues = [record.time(:), record.position, ...
            record.velocity, record.forceTotal, record.forceExcitation];
        assert(size(bodyValues, 1) == 8001 && ...
            size(bodyValues, 2) == 25 && ...
            all(isfinite(bodyValues), 'all'), ...
            'The full-duration refined body record is incomplete');
        writematrix(bodyValues, fullfile(outDir, ...
            sprintf('%s_body%d.csv', label, iBody)));
    end
    record = output.ptos(1);
    ptoValues = [record.time(:), record.position, record.velocity, ...
        record.forceInternalMechanics, record.powerInternalMechanics];
    assert(size(ptoValues, 1) == 8001 && ...
        size(ptoValues, 2) == 25 && all(isfinite(ptoValues), 'all'), ...
        'The full-duration refined PTO record is incomplete');
    writematrix(ptoValues, fullfile(outDir, [label '_pto.csv']));
    record = output.mooring(1);
    mooringValues = [record.time(:), record.position, ...
        record.velocity, record.forceMooring];
    assert(size(mooringValues, 1) == 8001 && ...
        size(mooringValues, 2) == 19 && ...
        all(isfinite(mooringValues), 'all'), ...
        'The full-duration refined mooring record is incomplete');
    writematrix(mooringValues, ...
        fullfile(outDir, [label '_mooring.csv']));
    close_system('RM3MoorDyn', 0);
end
end
