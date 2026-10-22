% Re-run the pinned corrected OWC model proposed upstream in PR #119.
% Run as a script so Simulink mask and variant initialization can resolve the
% WEC-Sim variables in the same workspace as initializeWecSim.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
coreRoot = fullfile(repoRoot, 'matlab-ref');
publishedDir = fullfile(repoRoot, 'applications', 'OWC', 'OrificeModel');
caseDir = fullfile(repoRoot, 'owc-orifice-corrected-case');
outputDir = fullfile(repoRoot, 'owc-orifice-corrected-output');
assert(isfile(fullfile(publishedDir, 'OWC_GBM.slx')) && ...
    isfile(fullfile(publishedDir, 'hydroData', 'test17a_clean.h5')));
assert(isfolder(fullfile(coreRoot, 'source')));
assert(~isfolder(caseDir), 'Do not reuse a previously modified model');
[copied, message] = copyfile(publishedDir, caseDir);
assert(copied, message);
mkdir(outputDir);

% Avoid opening the Mechanics Explorer or post-run plots on the CI runner.
inputFile = fullfile(caseDir, 'wecSimInputFile.m');
input = fileread(inputFile);
old = "simu.explorer = 'on';";
assert(count(input, old) == 1, 'Published input changed');
input = replace(input, old, "simu.explorer = 'off';");
fid = fopen(inputFile, 'w');
assert(fid ~= -1);
fprintf(fid, '%s', input);
fclose(fid);
fid = fopen(fullfile(caseDir, 'userDefinedFunctions.m'), 'w');
assert(fid ~= -1);
fprintf(fid, '%% Plotting disabled for the source comparison.\n');
fclose(fid);

addpath(genpath(fullfile(coreRoot, 'source')));
modelFile = fullfile(caseDir, 'OWC_GBM.slx');
load_system(modelFile);
sub = 'OWC_GBM/Flex Body/FlexibleBody';
sumBlock = [sub '/Piston input to GBM'];
assert(strcmp(get_param(sumBlock, 'BlockType'), 'Sum'));
assert(strcmp(get_param(sumBlock, 'Inputs'), '++'));
productPorts = get_param([sub '/Product4'], 'PortHandles');
sumPorts = get_param(sumBlock, 'PortHandles');
expectedInputs = {[sub '/Excitation'], [sub '/Gain']};
for i = 1:2
    line = get_param(sumPorts.Inport(i), 'Line');
    source = get_param(line, 'SrcBlockHandle');
    assert(strcmp(getfullname(source), expectedInputs{i}));
end
line = get_param(productPorts.Inport(2), 'Line');
source = get_param(line, 'SrcBlockHandle');
assert(strcmp(getfullname(source), sumBlock));
close_system('OWC_GBM', 0);
previousDir = pwd;
cd(caseDir);
run('initializeWecSim');
assert(simu.endTime == 130 && simu.dt == .005 && ...
    simu.dtOut == .005 && strcmp(simu.solver, 'ode23t') && ...
    waves.phaseSeed == 128 && waves.height == 1 && waves.period == 4, ...
    'Published OWC settings changed');
[~, model, ~] = fileparts(simu.simMechanicsFile);
assert(abs(str2double(get_param(model, 'MaxStep')) - .005) < 1e-12);
set_param(model, 'MaxStep', '0.00015625');
try
    % matlab-actions runs this script in the base workspace; the model's
    % masked blocks resolve the initialized body and variant values there.
    sim(simu.simMechanicsFile);
catch error
    fprintf(2, '%s\n', getReport(error, 'extended', 'hyperlinks', 'off'));
    rethrow(error);
end
run('stopWecSim');

time = output.bodies(1).time(:);
assert(numel(time) == 26001 && numel(Flex_out.time) == numel(time));
assert(max(abs(Flex_out.time(:) - time)) < 1e-10);
assert(numel(orificeData.time) > 800000 && ...
    all(diff(orificeData.time(:)) > 0));
assert(size(orificeData.signals.values, 2) == 5);
orifice = interp1(orificeData.time(:), ...
    orificeData.signals.values, time, 'linear');
assert(all(isfinite(orifice), 'all'));
assert(~any(orificeData.signals.values(:, 2)), ...
    'Corrected source exceeded the configured Mach limit');

moving = output.bodies(1);
writematrix([time, moving.position, moving.velocity], ...
    fullfile(outputDir, 'body.csv'));
writematrix([time, Flex_out.signals.values], ...
    fullfile(outputDir, 'flex.csv'));
writematrix([time, orifice], fullfile(outputDir, 'orifice.csv'));
writematrix([output.wave.time(:), output.wave.elevation(:)], ...
    fullfile(outputDir, 'wave.csv'));
writematrix([waves.omega(:), waves.amplitude(:), waves.dOmega(:), ...
    waves.phase(:)], fullfile(outputDir, 'components.csv'));
writematrix([trapz(orificeData.time(:), ...
    orificeData.signals.values(:, 5)) * 1000, ...
    max(abs(orificeData.signals.values(:, 4) / A2) / 343.2)], ...
    fullfile(outputDir, 'summary.csv'));
close_system(model, 0);
cd(previousDir);
