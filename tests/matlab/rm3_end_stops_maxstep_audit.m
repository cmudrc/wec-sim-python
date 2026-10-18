function rm3_end_stops_maxstep_audit(auditEndTime)
% Isolate ode45 MaxStep from the published End_Stops force/output grid.
if nargin == 0
    auditEndTime = 120;
end
assert(auditEndTime == 120 || auditEndTime == 400, ...
    'Unsupported audit duration');
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end

sourceDir = fullfile(repoRoot, 'applications', 'End_Stops');
if auditEndTime == 120
    outDir = fullfile(repoRoot, 'matlab-rm3-end-stops-maxstep');
else
    outDir = fullfile(repoRoot, 'matlab-rm3-end-stops-maxstep-400s');
end
mkdir(outDir);
for arm = 1:2
    if arm == 1
        label = 'published';
    else
        label = 'maxstep_0025';
    end
    if auditEndTime == 120
        caseName = ['End_Stops_maxstep_' label];
    else
        caseName = ['End_Stops_maxstep_400s_' label];
    end
    caseDir = fullfile(repoRoot, 'applications', caseName);
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
    if auditEndTime ~= 400
        contents = strrep(contents, 'simu.endTime = 400;', ...
            sprintf('simu.endTime = %d;', auditEndTime));
    end
    fid = fopen(inputFile, 'w');
    assert(fid > 0, 'Could not write bounded End_Stops input');
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, contents);
    clear cleanup;

    cd(caseDir);
    run('initializeWecSim');
    assert(simu.dt == .1 && simu.endTime == auditEndTime && ...
        simu.rampTime == 100 && strcmp(simu.solver, 'ode45'), ...
        'Published End_Stops settings changed');
    [~, modelName, ~] = fileparts(simu.simMechanicsFile);
    assert(abs(str2double(get_param(modelName, 'MaxStep')) - .1) < 1e-12, ...
        'Published ode45 maximum step changed');
    if arm == 2
        set_param(modelName, 'MaxStep', '0.025');
    elseif auditEndTime == 400
        % Retain the accepted solver times so the saved-grid added-mass
        % recurrence can be checked against actual off-grid steps.
        set_param(modelName, 'SaveTime', 'on', 'TimeSaveName', 'tout', ...
            'Decimation', '1', 'Refine', '1');
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
    if arm == 1 && auditEndTime == 400
        assert(exist('tout', 'var') == 1, 'Solver times were not logged');
        writematrix(tout(:), fullfile(outDir, 'published_internal_times.csv'));
    end
    run('stopWecSim');
    sampleCount = round(auditEndTime / simu.dt) + 1;
    assert(numel(output.bodies(1).time) == sampleCount && ...
        numel(output.ptos(1).time) == sampleCount, ...
        'Incomplete End_Stops audit');
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
    if arm == 1 && auditEndTime == 400
        for iBody = 1:2
            writematrix(output.bodies(iBody).acceleration, ...
                fullfile(outDir, sprintf('published_body%d_acceleration.csv', iBody)));
            applied = body(iBody).hydroForce.hf1.storage;
            writematrix(applied.output_forceAddedMass, ...
                fullfile(outDir, sprintf('published_body%d_added_force.csv', iBody)));
            writematrix(applied.hydroForce_fAddedMass, ...
                fullfile(outDir, sprintf('published_body%d_added_matrix.csv', iBody)));
        end
    end
    close_system(modelName, 0);
    clear output;
end
end
