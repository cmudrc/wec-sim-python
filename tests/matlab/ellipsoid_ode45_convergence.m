function ellipsoid_ode45_convergence(waveType, selectedLabel, selectedOutputStep)
% Screen the pinned ellipsoid ode45 application at successively finer clocks.
% The finest physical arm samples output every 0.005 s after the all-clock
% finest arm exhausted hosted runners. Work only in copied input folders.
waveType = string(waveType);
assert(any(waveType == ["Regular", "RegularCIC"]), ...
    'Choose a pinned ellipsoid ode45 wave case');
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(repoRoot, 'applications', 'Nonlinear_Hydro', 'hydroData');
cd(hydroDir);
if ~isfile('ellipsoid.h5')
    bemio;
end
caseDir = fullfile(repoRoot, 'applications', 'Nonlinear_Hydro', ...
    'ode45', waveType);
outDir = fullfile(repoRoot, 'matlab-ellipsoid-ode45-convergence', waveType);
mkdir(outDir);
steps = [0.005, 0.0025, 0.00125, 0.000625];
labels = {'step005', 'step0025', 'step00125', 'step000625'};
if nargin >= 2
    % Each workflow invocation uses a fresh MATLAB process so histories
    % accumulate on disk without retaining earlier simulation workspaces.
    selectedLabel = string(selectedLabel);
    selection = strcmp(labels, selectedLabel);
    assert(nnz(selection) == 1, 'Choose a supported source step');
    steps = steps(selection);
    labels = labels(selection);
end

for iStep = 1:numel(steps)
    step = steps(iStep);
    label = labels{iStep};
    outputStep = step;
    if nargin >= 3
        assert(step == 0.000625 && selectedOutputStep == 0.005, ...
            'Only the screened finest output clock may differ');
        outputStep = selectedOutputStep;
    end
    targetDir = fullfile(fileparts(caseDir), waveType + "_" + label);
    mkdir(targetDir);
    copyfile(fullfile(caseDir, 'wecSimInputFile.m'), ...
        fullfile(targetDir, 'wecSimInputFile.m'));
    copyfile(fullfile(caseDir, 'ellipsoid.slx'), ...
        fullfile(targetDir, 'ellipsoid.slx'));
    inputFile = fullfile(targetDir, 'wecSimInputFile.m');
    contents = fileread(inputFile);
    assert(contains(contents, 'simu.dt = 0.05;') && ...
        ~contains(contents, 'simu.dtOut') && ...
        ~contains(contents, 'simu.nonlinearDt') && ...
        ~contains(contents, 'simu.cicDt'), ...
        'The pinned ellipsoid ode45 sample settings changed');
    replacement = sprintf('simu.dt = %.6f;', step);
    if outputStep ~= step
        replacement = sprintf('%s\nsimu.dtOut = %.6f;', ...
            replacement, outputStep);
    end
    contents = strrep(contents, 'simu.dt = 0.05;', replacement);
    fid = fopen(inputFile, 'w');
    assert(fid > 0);
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, contents);
    clear cleanup;

    cd(targetDir);
    wecSim;
    maxStep = str2double(get_param('ellipsoid', 'MaxStep'));
    assert(strcmp(simu.solver, 'ode45') && ...
        abs(simu.dt - step) < 1e-12 && ...
        abs(simu.dtOut - outputStep) < 1e-12 && ...
        abs(simu.nonlinearDt - step) < 1e-12 && ...
        abs(simu.cicDt - step) < 1e-12 && ...
        abs(maxStep - step) < 1e-12 && ...
        simu.endTime == 150 && simu.rampTime == 50 && ...
        strcmp(waves.type, lowerFirst(waveType)) && ...
        waves.height == 4 && waves.period == 6 && ...
        body(1).nonlinearHydro == 2 && pto(1).damping == 1200000, ...
        'The refined ellipsoid ode45 settings changed');
    n = round(simu.endTime / outputStep) + 1;
    ratio = round(outputStep / step);
    assert(abs(ratio * step - outputStep) < 1e-12);
    bodyRecord = output.bodies(1);
    ptoRecord = output.ptos(1);
    waveRecord = output.wave;
    if numel(waveRecord.time) == n
        waveTime = waveRecord.time;
        waveElevation = waveRecord.elevation;
    else
        assert(numel(waveRecord.time) == round(simu.endTime / step) + 1);
        waveTime = waveRecord.time(1:ratio:end);
        waveElevation = waveRecord.elevation(1:ratio:end);
    end
    values = [bodyRecord.time(:), bodyRecord.position(:, 3), ...
        bodyRecord.velocity(:, 3), ptoRecord.position(:, 3), ...
        ptoRecord.velocity(:, 3), ptoRecord.forceInternalMechanics(:, 3), ...
        -ptoRecord.forceInternalMechanics(:, 3) .* ptoRecord.velocity(:, 3), ...
        waveElevation(:), bodyRecord.forceRestoring(:, 3), ...
        bodyRecord.forceExcitation(:, 3), bodyRecord.acceleration(:, 3)];
    assert(size(values, 1) == n && size(values, 2) == 11 && ...
        all(isfinite(values), 'all') && ...
        max(abs(bodyRecord.time(:) - ptoRecord.time(:))) < 1e-10 && ...
        numel(waveTime) == n && ...
        max(abs(bodyRecord.time(:) - waveTime(:))) < 1e-10, ...
        'The refined ellipsoid output is incomplete');
    writematrix(values, fullfile(outDir, label + "_history.csv"));
    writematrix([simu.dt, simu.dtOut, simu.nonlinearDt, ...
        simu.cicDt, maxStep], fullfile(outDir, label + "_settings.csv"));
    close_system('ellipsoid', 0);
end
end

function value = lowerFirst(value)
value = char(value);
value(1) = lower(value(1));
end
