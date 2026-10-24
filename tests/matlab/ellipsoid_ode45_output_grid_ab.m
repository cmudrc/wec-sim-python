function ellipsoid_ode45_output_grid_ab(waveType)
% Compare two copies of the pinned ode45 input with identical physics clocks.
% Output sampling changes from 0.00125 to 0.005 s; the simulation class
% also derives Paraview output timing from dtOut.
waveType = string(waveType);
assert(any(waveType == ["Regular", "RegularCIC"]));
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(repoRoot, 'applications', 'Nonlinear_Hydro', 'hydroData');
cd(hydroDir);
if ~isfile('ellipsoid.h5')
    bemio;
end
caseDir = fullfile(repoRoot, 'applications', 'Nonlinear_Hydro', ...
    'ode45', waveType);
outDir = fullfile(repoRoot, 'matlab-ellipsoid-ode45-output-grid-ab', waveType);
mkdir(outDir);
dt = 0.00125;
outputSteps = [dt, 0.005];
labels = {'all125', 'out005'};

for iArm = 1:numel(labels)
    dtOut = outputSteps(iArm);
    label = labels{iArm};
    targetDir = fullfile(fileparts(caseDir), waveType + "_ab_" + label);
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
    replacement = 'simu.dt = 0.00125;';
    if iArm == 2
        replacement = sprintf('%s\nsimu.dtOut = 0.005;', replacement);
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
        abs(simu.dt - dt) < 1e-12 && ...
        abs(simu.dtOut - dtOut) < 1e-12 && ...
        abs(simu.nonlinearDt - dt) < 1e-12 && ...
        abs(simu.cicDt - dt) < 1e-12 && ...
        abs(maxStep - dt) < 1e-12 && ...
        simu.endTime == 150 && simu.rampTime == 50 && ...
        strcmp(waves.type, lowerFirst(waveType)) && ...
        waves.height == 4 && waves.period == 6 && ...
        body(1).nonlinearHydro == 2 && pto(1).damping == 1200000, ...
        'The output-only A/B changed a physical case setting');
    n = round(simu.endTime / dtOut) + 1;
    ratio = round(dtOut / dt);
    bodyRecord = output.bodies(1);
    ptoRecord = output.ptos(1);
    waveRecord = output.wave;
    if numel(waveRecord.time) == n
        waveTime = waveRecord.time;
        waveElevation = waveRecord.elevation;
    else
        assert(numel(waveRecord.time) == round(simu.endTime / dt) + 1);
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
        'The output-only A/B history is incomplete');
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
