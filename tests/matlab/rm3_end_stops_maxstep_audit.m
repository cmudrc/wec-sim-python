function rm3_end_stops_maxstep_audit
% Isolate ode45 MaxStep from the published End_Stops force/output grid.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end

sourceDir = fullfile(repoRoot, 'applications', 'End_Stops');
outDir = fullfile(repoRoot, 'matlab-rm3-end-stops-maxstep');
mkdir(outDir);
for arm = 1:2
    if arm == 1
        label = 'published';
    else
        label = 'maxstep_0025';
    end
    caseDir = fullfile(repoRoot, 'applications', ...
        ['End_Stops_maxstep_' label]);
    assert(~isfolder(caseDir), 'Disposable End_Stops audit case exists');
    mkdir(caseDir);
    for name = {'wecSimInputFile.m', 'RM3.slx'}
        copyfile(fullfile(sourceDir, name{1}), ...
            fullfile(caseDir, name{1}));
    end
    inputFile = fullfile(caseDir, 'wecSimInputFile.m');
    contents = fileread(inputFile);
    assert(contains(contents, 'simu.endTime = 400;') && ...
        contains(contents, 'simu.dt = 0.1;'), ...
        'Pinned End_Stops input changed');
    contents = strrep(contents, 'simu.endTime = 400;', ...
        'simu.endTime = 120;');
    fid = fopen(inputFile, 'w');
    assert(fid > 0, 'Could not write bounded End_Stops input');
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, contents);
    clear cleanup;

    cd(caseDir);
    run('initializeWecSim');
    assert(simu.dt == .1 && simu.endTime == 120 && ...
        simu.rampTime == 100 && strcmp(simu.solver, 'ode45'), ...
        'Published End_Stops settings changed');
    [~, modelName, ~] = fileparts(simu.simMechanicsFile);
    assert(abs(str2double(get_param(modelName, 'MaxStep')) - .1) < 1e-12, ...
        'Published ode45 maximum step changed');
    if arm == 2
        set_param(modelName, 'MaxStep', '0.025');
    end
    maxStep = str2double(get_param(modelName, 'MaxStep'));
    writematrix([simu.dt, simu.endTime, simu.rampTime, maxStep], ...
        fullfile(outDir, [label '_settings.csv']));

    % Simulink resolves model variant controls from its global workspace.
    names = who;
    for j = 1:numel(names)
        assignin('base', names{j}, eval(names{j}));
    end
    sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
    run('stopWecSim');
    assert(numel(output.bodies(1).time) == 1201 && ...
        numel(output.ptos(1).time) == 1201, ...
        'Incomplete 120 s End_Stops audit');
    writematrix([output.wave.time(:), output.wave.elevation(:)], ...
        fullfile(outDir, [label '_wave.csv']));
    for iBody = 1:2
        record = output.bodies(iBody);
        writematrix([record.time(:), record.position(:, [1 3 5]), ...
            record.velocity(:, [1 3 5])], ...
            fullfile(outDir, sprintf('%s_body%d.csv', label, iBody)));
    end
    record = output.ptos(1);
    writematrix([record.time(:), record.position(:, 3), ...
        record.velocity(:, 3), record.forceInternalMechanics(:, 3)], ...
        fullfile(outDir, [label '_pto.csv']));
    close_system(modelName, 0);
    clear output;
end
end
