function wavestar_controller_clock_probe(end_time)
% Hold the NMPC clock at 50 ms while the ode8 plant runs at 25 ms.
% Work only in a disposable copy of the pinned application. This is a
% source-timing diagnostic, not a published-case motion-parity gate.
if nargin == 0
    end_time = 15.2;
end
assert(isnumeric(end_time) && isscalar(end_time) && isfinite(end_time) ...
    && end_time > 0 && abs(end_time / 0.05 - round(end_time / 0.05)) < 1e-9, ...
    'The requested end time must be a positive multiple of 50 ms');
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
source = fullfile(repo, 'applications', 'WECCCOMP');
scratch = fullfile(repo, 'wavestar-controller-clock-work');
out = fullfile(repo, 'wavestar-controller-clock-output');
assert(isfolder(source) && ~isfolder(scratch) && ~isfolder(out), ...
    'The pinned application or a clean probe directory is missing');
[copied, message] = copyfile(source, scratch);
assert(copied, message);
mkdir(out);
[copied, message] = copyfile(fullfile(repo, 'tests', 'fixtures', 'matlab', ...
    'wavestar_fixed_input_step', 'wavestar.h5'), ...
    fullfile(scratch, 'hydroData', 'wavestar.h5'));
assert(copied, message);
caseDir = fullfile(scratch, 'WECCCOMP_Nonlinear_Model_Predictive');

inputFile = fullfile(caseDir, 'wecSimInputFile.m');
input = fileread(inputFile);
assert(contains(input, 'simu.dt             = 50/1000;') && ...
    contains(input, "simu.solver         = 'ode8';") && ...
    contains(input, 'simu.cicEndTime     = 2;'), ...
    'The pinned WaveStar timing inputs changed');
input = strrep(input, 'simu.dt             = 50/1000;', ...
    sprintf('simu.dt             = 0.025;\nsimu.dtOut = 0.05;'));
input = strrep(input, 'simu.cicEndTime     = 2;', ...
    sprintf('simu.endTime = %.12g;\n    simu.cicEndTime     = 2;', end_time));
for name = ["Float", "Arm", "Frame", "BC", "Motor"]
    assert(contains(input, name + '.stl'), ...
        'The pinned WaveStar geometry name changed');
    input = strrep(input, name + '.stl', name + '.STL');
end
write_text(inputFile, input);

controllerFile = fullfile(caseDir, 'controller_init.m');
controller = fileread(controllerFile);
assert(contains(controller, 'dt              = simu.dt;'), ...
    'The pinned WaveStar controller step changed');
write_text(controllerFile, strrep(controller, ...
    'dt              = simu.dt;', 'dt              = 0.05;'));
plotFile = fullfile(caseDir, 'userDefinedFunctions.m');
assert(contains(fileread(plotFile), 'ii = find(Output_power.time==25);'), ...
    'The pinned WaveStar plot script changed');
write_text(plotFile, '% No plots in this bounded timing probe.');

addpath(genpath(fullfile(repo, 'matlab-ref', 'source')));
cd(caseDir);
load_system('WaveStar');
% The top Controller contains a From block fed by the global time signal;
% making that whole subsystem atomic would sever the Goto/From connection.
% Give only the command-producing strategy an atomic 50 ms boundary.
block = 'WaveStar/Controller/controller_rot/Control strategy';
assert(strcmp(get_param(block, 'TreatAsAtomicUnit'), 'off') && ...
    strcmp(get_param(block, 'SystemSampleTime'), '-1'), ...
    'The pinned control strategy no longer inherits its sample time');
set_param(block, 'TreatAsAtomicUnit', 'on', 'SystemSampleTime', '0.05');
% The observer, predictor, and power charts would otherwise keep updating
% at the plant rate despite the held strategy output.
root = sfroot;
charts = root.find('-isa', 'Stateflow.EMChart');
names = ["RTI_NMPC", "observer", "vectorExcMpred", "resistive", "Power"];
selected = strings(0);
for iChart = 1:numel(charts)
    chart = charts(iChart);
    if ~startsWith(string(chart.Path), "WaveStar/")
        continue;
    end
    script = string(chart.Script);
    for name = names
        if startsWith(strtrim(script), "function ") && ...
                contains(extractBefore(script, newline), name + "(")
            chart.ChartUpdate = 'DISCRETE';
            chart.SampleTime = '0.05';
            selected(end + 1) = name; %#ok<AGROW>
            break;
        end
    end
end
assert(numel(selected) == numel(names) && isempty(setxor(selected, names)), ...
    'The pinned controller function set changed');
filters = find_system('WaveStar', 'LookUnderMasks', 'all', ...
    'BlockType', 'DiscreteTransferFcn');
assert(numel(filters) == 3, 'The pinned controller filter set changed');
for iFilter = 1:numel(filters)
    name = string(get_param(filters{iFilter}, 'Name'));
    if ~contains(name, "Cmd torque to measured torque dynamics")
        set_param(filters{iFilter}, 'SampleTime', '0.05');
    end
end
wecSim;
assert(strcmp(simu.solver, 'ode8') && abs(simu.dt - 0.025) < 1e-12 ...
    && abs(simu.dtOut - 0.05) < 1e-12 && simu.endTime == end_time, ...
    'The source solver or output settings changed');

commandTime = cmd_ptoM.time(:);
command = squeeze(cmd_ptoM.signals.values);
command = command(:);
assert(numel(commandTime) == numel(command) && ...
    all(isfinite(commandTime)) && all(isfinite(command)), ...
    'The source command trace is invalid');
output_count = round(end_time / 0.05) + 1;
plant_count = round(end_time / 0.025) + 1;
if numel(commandTime) == output_count
    assert(max(abs(commandTime - (0:0.05:end_time)')) < 1e-9, ...
        'The 50 ms controller output grid changed');
elseif numel(commandTime) == plant_count
    assert(max(abs(commandTime - (0:0.025:end_time)')) < 1e-9, ...
        'The 25 ms controller output grid changed');
    % Logging may occur on every plant step; the held command must not.
    middle = commandTime(2:2:end) >= 10;
    updates = abs(command(2:2:end) - command(1:2:end-1));
    assert(max(updates(middle)) < 1e-9, ...
        'The controller command still updates between 50 ms ticks');
else
    error('Unexpected controller command count: %d', numel(commandTime));
end

writematrix([commandTime, command], fullfile(out, 'command.csv'));
writematrix([output.bodies(1).time(:), output.bodies(1).position, ...
    output.bodies(1).velocity], fullfile(out, 'float.csv'));
% The command reaches the plant through forceActuation. Internal mechanics
% is zero in this application and would hide an actuator timing difference.
writematrix([output.ptos(1).time(:), output.ptos(1).forceActuation, ...
    output.ptos(1).forceTotal], ...
    fullfile(out, 'pto_force.csv'));
% Save only the axial conjugate force/speed needed for an independent PTO
% work comparison. The full source controller still sees the live plant.
writematrix([output.ptos(1).time(:), ...
    output.ptos(1).forceActuation(:, 3), ...
    output.ptos(1).velocity(:, 3)], fullfile(out, 'pto_axial.csv'));
writematrix([waves.omega(:), waves.amplitude(:), waves.dOmega(:), ...
    waves.phase(:)], fullfile(out, 'components.csv'));
writematrix([output.wave.time(:), output.wave.elevation(:)], ...
    fullfile(out, 'wave.csv'));
fprintf('CLOCK_PROBE_PASSED endTime=%g commandSamples=%d plantStep=0.025 controllerStep=0.05\n', ...
    end_time, numel(commandTime));
close_system('WaveStar', 0);
end

function write_text(path, value)
file = fopen(path, 'w');
assert(file ~= -1, 'Could not write disposable WaveStar input');
cleanup = onCleanup(@() fclose(file));
fprintf(file, '%s', value);
end
