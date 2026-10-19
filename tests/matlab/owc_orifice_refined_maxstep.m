% Reproduce the pinned OWC OrificeModel with only ode23t MaxStep reduced.
% Invoke in MATLAB's base workspace after preparing a disposable case copy:
% coreRoot = '/path/to/wec-sim-matlab-reference';
% caseDir = '/tmp/owc-audit/OWC/OrificeModel';
% outputDir = '/tmp/owc-audit';
% auditMaxStep = 0.000625; % or 0.0003125 or 0.00015625
% owc_orifice_refined_maxstep
assert(exist('coreRoot', 'var') == 1 && exist('caseDir', 'var') == 1 && ...
    exist('outputDir', 'var') == 1 && exist('auditMaxStep', 'var') == 1);
assert(any(abs(auditMaxStep - [0.000625, 0.0003125, 0.00015625]) < 1e-12), ...
    'Use one of the documented refined MaxStep values');
assert(isfolder(fullfile(coreRoot, 'source')) && ...
    isfile(fullfile(caseDir, 'OWC_GBM.slx')) && ...
    isfile(fullfile(caseDir, 'hydroData', 'test17a_clean.h5')));
if ~isfolder(outputDir)
    mkdir(outputDir);
end
addpath(genpath(fullfile(coreRoot, 'source')));
previousDir = pwd;
cd(caseDir);
run('initializeWecSim');
assert(simu.endTime == 130 && simu.dt == .005 && ...
    simu.dtOut == .005 && strcmp(simu.solver, 'ode23t') && ...
    waves.phaseSeed == 128, 'Pinned OWC source settings changed');
[~, model, ~] = fileparts(simu.simMechanicsFile);
assert(abs(str2double(get_param(model, 'MaxStep')) - .005) < 1e-12);
set_param(model, 'MaxStep', sprintf('%.8f', auditMaxStep));
assert(abs(str2double(get_param(model, 'MaxStep')) - auditMaxStep) < 1e-12);
sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
assert(abs(str2double(get_param(model, 'MaxStep')) - auditMaxStep) < 1e-12);
run('stopWecSim');

moving = output.bodies(1);
assert(numel(moving.time) == 26001 && ...
    numel(output.wave.time) == 26001 && ...
    numel(orificeData.time) > 400000, 'Incomplete OWC source run');
writematrix([moving.time(:), moving.position, moving.velocity, ...
    moving.acceleration], fullfile(outputDir, 'body.csv'));
writematrix([Flex_out.time(:), Flex_out.signals.values], ...
    fullfile(outputDir, 'flex.csv'));
writematrix([orificeData.time(:), orificeData.signals.values], ...
    fullfile(outputDir, 'orifice.csv'));
writematrix([output.wave.time(:), output.wave.elevation(:)], ...
    fullfile(outputDir, 'wave.csv'));
writematrix([waves.omega(:), waves.amplitude(:), waves.dOmega(:), ...
    waves.phase(:)], fullfile(outputDir, 'components.csv'));
close_system(model, 0);
cd(previousDir);
