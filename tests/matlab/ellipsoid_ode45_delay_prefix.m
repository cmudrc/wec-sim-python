function ellipsoid_ode45_delay_prefix(sourceDir, applicationsDir, outputDir)
% Log the active acceleration Transport Delay in the published ode45 case.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1
    sourceDir = fullfile(repoRoot, 'matlab-ref', 'source');
end
if nargin < 2
    applicationsDir = fullfile(repoRoot, 'applications');
end
if nargin < 3
    outputDir = fullfile(repoRoot, 'matlab-ellipsoid-ode45-delay-prefix');
end
addpath(genpath(sourceDir));
hydroDir = fullfile(applicationsDir, 'Nonlinear_Hydro', 'hydroData');
cd(hydroDir);
if ~isfile('ellipsoid.h5')
    bemio;
end
caseDir = fullfile(applicationsDir, 'Nonlinear_Hydro', 'ode45', 'Regular');
cd(caseDir);
run('initializeWecSim');
assert(strcmp(simu.solver, 'ode45') && simu.endTime == 150 && ...
    abs(simu.dt - .05) < 1e-12 && abs(simu.dtOut - .05) < 1e-12 && ...
    abs(simu.nonlinearDt - .05) < 1e-12, ...
    'The pinned published input has changed');

model = 'ellipsoid';
assert(abs(str2double(get_param(model, 'MaxStep')) - .05) < 1e-12);
set_param(model, 'StopTime', '2', 'SignalLogging', 'on', ...
    'SignalLoggingName', 'logsout', 'ReturnWorkspaceOutputs', 'on');
block = [model '/Rigid Body/Hydrodynamic Body/Transport Delay'];
assert(abs(str2double(get_param(block, 'DelayTime')) - 1e-7) < 1e-16);
ports = get_param(block, 'PortHandles');
for item = 1:2
    if item == 1
        line = get_param(ports.Inport(1), 'Line');
        name = 'delay_input';
    else
        line = get_param(ports.Outport(1), 'Line');
        name = 'delay_output';
    end
    assert(line ~= -1, 'The source delay is disconnected');
    sourcePort = get_param(line, 'SrcPortHandle');
    set_param(sourcePort, 'DataLogging', 'on', ...
        'DataLoggingNameMode', 'Custom', 'DataLoggingName', name);
end

% Model variants resolve in the base workspace when run from a function.
names = who;
for i = 1:numel(names)
    assignin('base', names{i}, eval(names{i}));
end
simOut = sim(simu.simMechanicsFile);
logs = simOut.logsout;
assert(logs.numElements == 2, 'Unexpected signal log');
if ~isfolder(outputDir)
    mkdir(outputDir);
end
for name = {'delay_input', 'delay_output'}
    signal = logs.getElement(name{1});
    times = signal.Values.Time(:);
    values = reshape(signal.Values.Data, numel(times), []);
    assert(isequal(size(values), [41, 6]) && ...
        max(abs(times - (0:40)' * .05)) < 1e-13, ...
        'The source delay prefix is incomplete');
    writematrix([times, values], fullfile(outputDir, [name{1} '.csv']));
end
close_system(model, 0);
end
