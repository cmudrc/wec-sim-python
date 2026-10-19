function oswec_core_step_work_audit(coreRoot, outputDir, step)
% Reproduce the published OSWEC sea with only the fixed ode4 step changed.
% Call each step in a separate MATLAB process and fresh scratch directory.
assert(ismember(step, [0.05, 0.025, 0.0125]), ...
    'Unsupported OSWEC audit step');
caseDir = fullfile(coreRoot, 'examples', 'OSWEC');
assert(isfolder(caseDir) && isfolder(fullfile(coreRoot, 'source')));
assert(~isfolder(outputDir), 'Use a new output directory for each source run');
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
contents = strrep(contents, 'simu.dt = 0.1;', ...
    sprintf('simu.dt = %.12g;', step));
contents = strrep(contents, "simu.explorer = 'on';", ...
    "simu.explorer = 'off';");
contents = strrep(contents, 'waves.spread = [0.1,0.2,0.7];', ...
    sprintf('waves.spread = [0.1,0.2,0.7];\nwaves.phaseSeed = 1;'));
fid = fopen(inputFile, 'w');
assert(fid ~= -1);
fprintf(fid, '%s', contents);
fclose(fid);

% Plotting runs only after the simulation; it is irrelevant to trajectories.
fid = fopen(fullfile(scratch, 'userDefinedFunctions.m'), 'w');
assert(fid ~= -1);
fprintf(fid, '%% No post-run plots in this local step audit.\n');
fclose(fid);

addpath(genpath(fullfile(coreRoot, 'source')));
previousDir = pwd;
cd(scratch);
wecSim;
cd(previousDir);
assert(simu.dt == step && simu.endTime == 400 && ...
    strcmp(simu.solver, 'ode4') && waves.phaseSeed == 1);
assert(numel(output.ptos) == 1 && ...
    numel(output.ptos(1).time) == round(400 / step) + 1);
pto = output.ptos(1);
writematrix([pto.time(:), pto.velocity(:, 5), ...
    pto.forceInternalMechanics(:, 5), ...
    pto.powerInternalMechanics(:, 5)], fullfile(outputDir, 'pto.csv'));
writematrix([output.wave.time(:), output.wave.elevation(:)], ...
    fullfile(outputDir, 'wave.csv'));
flap = output.bodies(1);
writematrix([flap.time(:), flap.position, flap.velocity], ...
    fullfile(outputDir, 'body.csv'));
end
