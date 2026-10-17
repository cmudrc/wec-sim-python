% Isolate ode4 FixedStep from WEC-Sim's published 0.1 s force/output clocks.
% Invoke as a script in the base workspace:
%   caseName='float'; mode='published'; rm3_pto_extension_fixedstep_only
assert(any(strcmp(caseName, {'float', 'spar'})));
assert(any(strcmp(mode, {'published', 'fixedstep0125'})));
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
auditRoot = getenv('WEC_SIM_RM3_PTO_FIXEDSTEP_ROOT');
if isempty(auditRoot)
    auditRoot = repoRoot;
end
addpath(genpath(fullfile(auditRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(auditRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end
originalDir = fullfile(auditRoot, 'applications', 'RM3_PTO_Extension', caseName);
contents = fileread(fullfile(originalDir, 'wecSimInputFile.m'));
assert(contains(contents, 'simu.dt = 0.1;') && ...
    contains(contents, 'simu.endTime = 30;') && ...
    contains(contents, "waves = waveClass('noWaveCIC');") && ...
    contains(contents, 'pto(1).extension.PositionTargetValue = 5;') && ...
    ~contains(contents, 'simu.dtOut =') && ...
    ~contains(contents, 'simu.cicDt ='), ...
    'Pinned RM3 PTO-extension input changed');
caseDir = fullfile(fileparts(originalDir), [caseName '_fixedstep_' mode]);
assert(~isfolder(caseDir), 'Disposable PTO-extension audit case already exists');
mkdir(caseDir);
copyfile(fullfile(originalDir, 'RM3.slx'), fullfile(caseDir, 'RM3.slx'));
copyfile(fullfile(originalDir, 'wecSimInputFile.m'), ...
    fullfile(caseDir, 'wecSimInputFile.m'));
cd(caseDir);
run('initializeWecSim');
assert(strcmp(simu.solver, 'ode4') && simu.dt == 0.1 && ...
    simu.dtOut == 0.1 && simu.cicDt == 0.1 && ...
    simu.nonlinearDt == 0.1 && simu.endTime == 30 && ...
    strcmp(waves.type, 'noWaveCIC') && ...
    pto(1).extension.PositionTargetValue == 5 && ...
    pto(1).stiffness == 0 && pto(1).damping == 0, ...
    'Published settings changed');
model = 'RM3';
assert(abs(str2double(get_param(model, 'FixedStep')) - 0.1) < 1e-12);
if strcmp(mode, 'fixedstep0125')
    set_param(model, 'FixedStep', '0.0125');
end
step = str2double(get_param(model, 'FixedStep'));
fprintf('SOURCE_SETTINGS case=%s mode=%s dt=%.9g dtOut=%.9g nonlinearDt=%.9g cicDt=%.9g FixedStep=%.9g\n', ...
    caseName, mode, simu.dt, simu.dtOut, simu.nonlinearDt, simu.cicDt, step);
sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
assert(abs(str2double(get_param(model, 'FixedStep')) - step) < 1e-12, ...
    'Simulink changed FixedStep during the run');
run('stopWecSim');
float = output.bodies(1);
spar = output.bodies(2);
pto = output.ptos(1);
values = [float.time(:), float.position(:, 3), float.velocity(:, 3), ...
    spar.position(:, 3), spar.velocity(:, 3), ...
    pto.position(:, 3), pto.velocity(:, 3), ...
    pto.forceInternalMechanics(:, 3), pto.powerInternalMechanics(:, 3)];
assert(isequal(size(values), [301, 9]) && all(isfinite(values), 'all'), ...
    'Incomplete PTO-extension source trajectory');
outDir = fullfile(auditRoot, 'matlab-rm3-pto-fixedstep-only', caseName);
if ~isfolder(outDir)
    mkdir(outDir);
end
writematrix(values, fullfile(outDir, [mode '.csv']));
writematrix([simu.dt, simu.dtOut, simu.nonlinearDt, simu.cicDt, step], ...
    fullfile(outDir, [mode '_settings.csv']));
close_system(model, 0);
