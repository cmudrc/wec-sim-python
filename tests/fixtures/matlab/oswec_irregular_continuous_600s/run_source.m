function run_source(source_dir, case_dir, output_dir, expected_dt)
% Export the independently simulated flap, PTO, and incident-wave histories.
% case_dir is a copy of the pinned application with the matching archived
% wecSimInputFile.m, calcIndex.m, and userDefinedFunctions.m installed.
addpath(genpath(source_dir));
if ~isfolder(output_dir), mkdir(output_dir); end
cd(case_dir);
waveFlag = 'irregular';
clock = tic;
wecSim;
assert(strcmp(simu.solver, 'ode4') && simu.dt == expected_dt);
assert(simu.endTime == 600 && waves.phaseSeed == 1);
assert(body(1).yaw.option == 1 && body(1).yaw.threshold == 0);
assert(body(1).variableHydro.option == 0);
trace = [output.bodies(1).time(:), output.bodies(1).position(:, 6), ...
    output.bodies(1).velocity(:, 6), ...
    output.ptos(1).forceInternalMechanics(:, 5), ...
    output.ptos(1).powerInternalMechanics(:, 5)];
assert(size(trace, 1) == round(600 / expected_dt) + 1);
assert(all(isfinite(trace), 'all'));
writematrix(trace, fullfile(output_dir, 'trajectory.csv'));
writematrix(waves.waveAmpTime, fullfile(output_dir, 'wave.csv'));
metrics = struct('elapsed_s', toc(clock), ...
    'yaw_range_rad', [min(trace(:, 2)), max(trace(:, 2))], ...
    'pto_work_j', -trapz(trace(:, 1), trace(:, 5)));
fid = fopen(fullfile(output_dir, 'metrics.json'), 'w');
assert(fid ~= -1);
fprintf(fid, '%s\n', jsonencode(metrics));
fclose(fid);
close_system('OSWEC', 0);
end
