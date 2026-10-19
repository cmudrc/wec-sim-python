% Refine Simulink ode4 while retaining every published WEC-Sim 0.1 s clock.
% Invoke as a script in the base workspace with coreRoot, outputDir, and
% fixedStep set. Simulink resolves model parameters from that workspace.
assert(exist('coreRoot', 'var') == 1 && exist('outputDir', 'var') == 1 && ...
    exist('fixedStep', 'var') == 1);
assert(ismember(fixedStep, [0.1, 0.05, 0.025]));
caseDir = fullfile(coreRoot, 'examples', 'OSWEC');
assert(isfolder(caseDir) && isfolder(fullfile(coreRoot, 'source')));
assert(~isfolder(outputDir), 'Use a new output directory');
mkdir(outputDir);
scratch = fullfile(outputDir, 'case');
[copied, message] = copyfile(caseDir, scratch);
assert(copied, message);

inputFile = fullfile(scratch, 'wecSimInputFile.m');
contents = fileread(inputFile);
assert(contains(contents, 'simu.dt = 0.1;') && ...
    contains(contents, "simu.explorer = 'on';") && ...
    contains(contents, 'waves.spread = [0.1,0.2,0.7];') && ...
    ~contains(contents, 'waves.phaseSeed'), ...
    'The pinned OSWEC input has changed');
contents = strrep(contents, "simu.explorer = 'on';", ...
    "simu.explorer = 'off';");
contents = strrep(contents, 'waves.spread = [0.1,0.2,0.7];', ...
    sprintf('waves.spread = [0.1,0.2,0.7];\nwaves.phaseSeed = 1;'));
fid = fopen(inputFile, 'w');
assert(fid ~= -1);
fprintf(fid, '%s', contents);
fclose(fid);
fid = fopen(fullfile(scratch, 'userDefinedFunctions.m'), 'w');
assert(fid ~= -1);
fprintf(fid, '%% No post-run plots in this solver-step audit.\n');
fclose(fid);

addpath(genpath(fullfile(coreRoot, 'source')));
previousDir = pwd;
cd(scratch);
run('initializeWecSim');
assert(simu.dt == 0.1 && simu.dtOut == 0.1 && ...
    simu.cicDt == 0.1 && simu.nonlinearDt == 0.1 && ...
    simu.endTime == 400 && simu.rampTime == 100 && ...
    strcmp(simu.solver, 'ode4') && waves.phaseSeed == 1, ...
    'Published WEC-Sim clocks or input changed');
[~, model, ~] = fileparts(simu.simMechanicsFile);
assert(abs(str2double(get_param(model, 'FixedStep')) - 0.1) < 1e-12);
set_param(model, 'FixedStep', num2str(fixedStep, 16));
assert(abs(str2double(get_param(model, 'FixedStep')) - fixedStep) < 1e-12);
sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
assert(abs(str2double(get_param(model, 'FixedStep')) - fixedStep) < 1e-12);
run('stopWecSim');

pto = output.ptos(1);
wave = output.wave;
flap = output.bodies(1);
assert(numel(pto.time) == 4001 && numel(wave.time) == 4001 && ...
    numel(flap.time) == 4001, 'Unexpected output clock');
writematrix([pto.time(:), pto.velocity(:, 5), ...
    pto.forceInternalMechanics(:, 5), ...
    pto.powerInternalMechanics(:, 5)], fullfile(outputDir, 'pto.csv'));
writematrix([wave.time(:), wave.elevation(:)], ...
    fullfile(outputDir, 'wave.csv'));
writematrix([flap.time(:), flap.position, flap.velocity], ...
    fullfile(outputDir, 'body.csv'));
writematrix([simu.dt, simu.dtOut, simu.cicDt, ...
    simu.nonlinearDt, fixedStep], fullfile(outputDir, 'settings.csv'));
close_system(model, 0);
cd(previousDir);
