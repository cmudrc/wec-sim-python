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
    % MoorDyn writes at the coupling step, while the body comparison uses
    % 0.01 s samples. Keep the three independently reported line tensions
    % on that common grid rather than inferring them from net mooring load.
    lines = output.moorDyn(1).Lines;
    stride = round(0.01 / step);
    assert(abs(stride * step - 0.01) < 1e-12 && ...
        numel(lines.Time) == round(simu.endTime / step) + 1, ...
        'The refined MoorDyn output grid changed');
    samples = 1:stride:numel(lines.Time);
    fairlead = [reshape(lines.Time(samples), [], 1), ...
        reshape(lines.FairTen1(samples), [], 1), ...
        reshape(lines.FairTen2(samples), [], 1), ...
        reshape(lines.FairTen3(samples), [], 1)];
    assert(size(fairlead, 1) == 8001 && ...
        max(abs(fairlead(:, 1) - mooringValues(:, 1))) < 1e-8 && ...
        all(isfinite(fairlead), 'all'), ...
        'The refined fairlead record is incomplete');
    writematrix(fairlead, fullfile(outDir, [label '_fairlead.csv']));
    % Save complete line geometry and segment loads at the 0.1 s published
    % visualization cadence. These histories come from this refined MATLAB
    % trajectory, not from the published-step source-pose VTP baseline.
    assert(isequal(mooring(1).moorDynNodes, [21 21 21]), ...
        'The three-line MoorDyn topology changed');
    frameStride = round(0.1 / step);
    assert(abs(frameStride * step - 0.1) < 1e-12);
    for iLine = 1:3
        lineRecord = output.moorDyn(1).(sprintf('Line%d', iLine));
        frames = 1:frameStride:numel(lineRecord.Time);
        nodeCount = mooring(1).moorDynNodes(iLine);
        history = zeros(numel(frames), 4 * nodeCount);
        history(:, 1) = reshape(lineRecord.Time(frames), [], 1);
        column = 2;
        for iNode = 0:nodeCount - 1
            for axis = {'px', 'py', 'pz'}
                field = sprintf('Node%d%s', iNode, axis{1});
                values = lineRecord.(field);
                history(:, column) = reshape(values(frames), [], 1);
                column = column + 1;
            end
        end
        for iSegment = 1:nodeCount - 1
            values = lineRecord.(sprintf('Seg%dTe', iSegment));
            history(:, column) = reshape(values(frames), [], 1);
            column = column + 1;
        end
        assert(column == size(history, 2) + 1 && ...
            size(history, 1) == 801 && ...
            max(abs(history(:, 1) - (0:0.1:80)')) < 1e-8 && ...
            all(isfinite(history), 'all'), ...
            'The refined MoorDyn line history is incomplete');
        writematrix(history, fullfile(outDir, ...
            sprintf('%s_line%d.csv', label, iLine)));
    end
    close_system('RM3MoorDyn', 0);
end
end
