% Isolate nonlinear-force sampling from the ellipsoid ode45 solver step.
% Invoke as a script: waveCase='Regular'; mode='fine'; ellipsoid_ode45_nonlinear_clock
assert(any(strcmp(waveCase, {'Regular', 'RegularCIC'})));
assert(any(strcmp(mode, {'fine', 'coarse_force', 'coarse_force_cic'})));
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
auditRoot = getenv('WEC_SIM_ELLIPSOID_CLOCK_ROOT');
if isempty(auditRoot)
    auditRoot = repoRoot;
end
addpath(genpath(fullfile(auditRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(auditRoot, 'applications', 'Nonlinear_Hydro', 'hydroData');
originalDir = fullfile(auditRoot, 'applications', 'Nonlinear_Hydro', 'ode45', waveCase);
original = fileread(fullfile(originalDir, 'wecSimInputFile.m'));
assert(contains(original, 'simu.dt = 0.05;') && ...
    ~contains(original, 'simu.dtOut') && ...
    ~contains(original, 'simu.nonlinearDt') && ...
    ~contains(original, 'simu.cicDt') && ...
    contains(original, 'simu.endTime=150;') && ...
    contains(original, 'body(1).nonlinearHydro = 2;'), ...
    'Pinned ellipsoid input changed');
cd(hydroDir);
if ~isfile('ellipsoid.h5')
    bemio;
end

% First restore the published nonlinear-force clock while holding solver,
% wave, output, and radiation clocks fixed. For RegularCIC, a third arm
% restores the published convolution clock as a separate change.
caseDir = fullfile(fileparts(originalDir), [waveCase '_clock_' mode]);
assert(~isfolder(caseDir), 'Disposable clock-audit case already exists');
mkdir(caseDir);
copyfile(fullfile(originalDir, 'ellipsoid.slx'), fullfile(caseDir, 'ellipsoid.slx'));
nonlinearDt = 0.005;
if ~strcmp(mode, 'fine')
    nonlinearDt = 0.05;
end
cicDt = 0.005;
if strcmp(mode, 'coarse_force_cic')
    cicDt = 0.05;
end
contents = strrep(original, 'simu.dt = 0.05;', sprintf([ ...
    'simu.dt = 0.005;\nsimu.dtOut = 0.05;\n' ...
    'simu.cicDt = %.3f;\nsimu.nonlinearDt = %.3f;'], cicDt, nonlinearDt));
fid = fopen(fullfile(caseDir, 'wecSimInputFile.m'), 'w');
assert(fid > 0);
fprintf(fid, '%s', contents);
fclose(fid);
cd(caseDir);
run('initializeWecSim');
assert(strcmp(simu.solver, 'ode45') && ...
    abs(simu.dt - 0.005) < 1e-12 && ...
    abs(simu.dtOut - 0.05) < 1e-12 && ...
    abs(simu.cicDt - cicDt) < 1e-12 && ...
    abs(simu.nonlinearDt - nonlinearDt) < 1e-12 && ...
    simu.endTime == 150 && simu.rampTime == 50, ...
    'An unrelated sample clock changed');
model = 'ellipsoid';
step = str2double(get_param(model, 'MaxStep'));
assert(abs(step - 0.005) < 1e-12);
fprintf('SOURCE_SETTINGS case=%s mode=%s dt=%.9g dtOut=%.9g nonlinearDt=%.9g cicDt=%.9g MaxStep=%.9g\n', ...
    waveCase, mode, simu.dt, simu.dtOut, simu.nonlinearDt, simu.cicDt, step);
sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
run('stopWecSim');
record = output.bodies(1);
bodyValues = [record.time(:), record.position, record.velocity, ...
    record.forceRestoring, record.forceExcitation];
pto = output.ptos(1);
ptoValues = [pto.time(:), pto.position, pto.velocity, ...
    pto.forceInternalMechanics];
waveValues = [output.wave.time(:), output.wave.elevation(:)];
assert(isequal(size(bodyValues), [3001, 25]) && ...
    isequal(size(ptoValues), [3001, 19]) && ...
    isequal(size(waveValues), [3001, 2]) && ...
    all(isfinite(bodyValues), 'all') && ...
    all(isfinite(ptoValues), 'all') && ...
    all(isfinite(waveValues), 'all'), 'Incomplete source output');
outDir = fullfile(auditRoot, 'matlab-ellipsoid-ode45-nonlinear-clock', waveCase);
if ~isfolder(outDir)
    mkdir(outDir);
end
writematrix(bodyValues, fullfile(outDir, [mode '_body.csv']));
writematrix(ptoValues, fullfile(outDir, [mode '_pto.csv']));
writematrix(waveValues, fullfile(outDir, [mode '_wave.csv']));
writematrix([simu.dt, simu.dtOut, simu.nonlinearDt, simu.cicDt, step], ...
    fullfile(outDir, [mode '_settings.csv']));
close_system(model, 0);
