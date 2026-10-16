function rm3_moordyn_initial_step_80s_audit
% Compare ode45 InitialStep=auto and 0.01 over the full seeded RM3 case.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end
sourceDir = fullfile(repoRoot, 'applications', ...
    'Paraview_Visualization', 'RM3_MoorDyn_Viz');
outDir = fullfile(repoRoot, 'matlab-rm3-initial-step-80s-audit');
mkdir(outDir);
for arm = 1:2
    if arm == 1
        label = 'auto';
    else
        label = 'initial_001';
    end
    caseDir = fullfile(fileparts(sourceDir), ...
        ['RM3_MoorDyn_Viz_initial_step_80s_' label]);
    assert(~isfolder(caseDir), 'The disposable audit case already exists');
    mkdir(caseDir);
    for name = {'wecSimInputFile.m', 'RM3MoorDyn.slx', ...
            'userDefinedFunctions.m'}
        copyfile(fullfile(sourceDir, name{1}), fullfile(caseDir, name{1}));
    end
    mkdir(fullfile(caseDir, 'Mooring'));
    copyfile(fullfile(sourceDir, 'Mooring', 'lines.txt'), ...
        fullfile(caseDir, 'Mooring', 'lines.txt'));
    inputFile = fullfile(caseDir, 'wecSimInputFile.m');
    contents = fileread(inputFile);
    changes = {
        'simu.paraview.option = 1;', 'simu.paraview.option = 0;';
        'simu.mode = ''accelerator'';', 'simu.mode = ''normal'';';
        'waves.period = 8;', sprintf('waves.period = 8;\nwaves.phaseSeed = 1;')
    };
    for i = 1:size(changes, 1)
        assert(contains(contents, changes{i, 1}), 'Pinned input changed');
        contents = strrep(contents, changes{i, 1}, changes{i, 2});
    end
    fid = fopen(inputFile, 'w');
    assert(fid > 0);
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, contents);
    clear cleanup;
    cd(caseDir);
    run('initializeWecSim');
    assert(simu.dt == .01 && simu.dtOut == .1 && simu.endTime == 80 && ...
        strcmp(simu.solver, 'ode45') && strcmp(simu.mode, 'normal') && ...
        waves.phaseSeed == 1 && mooring(1).moorDyn == 1, ...
        'Full source settings changed');
    [~, modelName, ~] = fileparts(simu.simMechanicsFile);
    assert(strcmp(get_param(modelName, 'InitialStep'), 'auto') && ...
        abs(str2double(get_param(modelName, 'MaxStep')) - .01) < 1e-12, ...
        'Pinned solver settings changed');
    if arm == 2
        set_param(modelName, 'InitialStep', '0.01');
    end
    % Simulink resolves model variant controls in its global workspace.
    names = who;
    for j = 1:numel(names)
        assignin('base', names{j}, eval(names{j}));
    end
    sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
    run('stopWecSim');
    assert(numel(output.bodies(1).time) == 801 && ...
        numel(output.moorDyn(1).Lines.Time) == 8001, ...
        'Incomplete full source run');
    writematrix(waves.phase(:), fullfile(outDir, [label '_phase.csv']));
    writematrix([output.wave.time(:), output.wave.elevation(:)], ...
        fullfile(outDir, [label '_wave.csv']));
    for iBody = 1:2
        record = output.bodies(iBody);
        writematrix([record.time(:), record.position, record.velocity], ...
            fullfile(outDir, sprintf('%s_body%d.csv', label, iBody)));
    end
    record = output.ptos(1);
    writematrix([record.time(:), record.position, record.velocity, ...
        record.forceInternalMechanics, record.powerInternalMechanics], ...
        fullfile(outDir, [label '_pto.csv']));
    record = output.mooring(1);
    writematrix([record.time(:), record.position, record.velocity, ...
        record.forceMooring], fullfile(outDir, [label '_mooring.csv']));
    lines = output.moorDyn(1).Lines;
    writematrix([lines.Time(:), lines.FairTen1(:), ...
        lines.FairTen2(:), lines.FairTen3(:)], ...
        fullfile(outDir, [label '_fairlead.csv']));
    close_system(modelName, 0);
end
end
