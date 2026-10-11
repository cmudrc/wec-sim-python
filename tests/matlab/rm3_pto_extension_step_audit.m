function rm3_pto_extension_step_audit()
% Hold the two published no-wave PTO-extension cases fixed while refining
% only the MATLAB output/solver step. This is a source sensitivity audit,
% not a Python trajectory parity claim.
reference_model_baseline('RM3_PTO_Extension');
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
baselineDir = fullfile(repoRoot, 'matlab-reference-model-output');
outDir = fullfile(repoRoot, 'matlab-rm3-pto-extension-step-audit');
mkdir(outDir);
metrics = struct();
metrics.channels = {'floatHeaveM', 'floatHeaveSpeedMPerS', ...
    'sparHeaveM', 'sparHeaveSpeedMPerS', 'ptoStrokeM', ...
    'ptoStrokeSpeedMPerS', 'ptoForceN', 'ptoPowerW'};
metrics.stepsS = [0.1, 0.05, 0.025];

for caseName = ["float", "spar"]
    label = char(caseName);
    caseDir = fullfile(repoRoot, 'applications', 'RM3_PTO_Extension', label);
    inputFile = fullfile(caseDir, 'wecSimInputFile.m');
    original = fileread(inputFile);
    assert(contains(original, 'simu.dt = 0.1;') && ...
        contains(original, 'simu.endTime = 30;') && ...
        contains(original, "waves = waveClass('noWaveCIC');") && ...
        contains(original, 'pto(1).extension.PositionTargetValue = 5;'), ...
        'The pinned PTO-extension input changed');
    coarse = baseline_trajectory(baselineDir, label);
    assert(isequal(size(coarse), [301, 9]), ...
        'The published PTO-extension baseline is incomplete');
    writematrix(coarse, fullfile(outDir, label + "_dt_0p1.csv"));
    previous = coarse;
    levelNames = ["coarseToMid", "midToFine"];
    steps = [0.05, 0.025];
    for level = 1:numel(steps)
        dt = steps(level);
        write_input(inputFile, strrep(original, 'simu.dt = 0.1;', ...
            sprintf('simu.dt = %.3f;', dt)));
        cd(caseDir);
        wecSim;
        assert(strcmp(simu.solver, 'ode4') && simu.dt == dt && ...
            simu.endTime == 30 && strcmp(waves.type, 'noWaveCIC') && ...
            pto(1).damping == 0 && pto(1).stiffness == 0 && ...
            pto(1).extension.PositionTargetValue == 5, ...
            'Step refinement changed a published input');
        current = trajectory(output);
        assert(isequal(size(current), [round(30 / dt) + 1, 9]) && ...
            max(abs(current(1:2:end, 1) - previous(:, 1))) < 1e-10, ...
            'The refined PTO-extension output is incomplete');
        if dt == 0.05
            suffix = '0p05';
        else
            suffix = '0p025';
        end
        writematrix(current, fullfile(outDir, ...
            label + "_dt_" + suffix + ".csv"));
        delta = current(1:2:end, 2:9) - previous(:, 2:9);
        metrics.(label).(char(levelNames(level))) = max(abs(delta), [], 1);
        previous = current;
        close_system('RM3', 0);
    end
    write_input(inputFile, original);
end

fid = fopen(fullfile(outDir, 'metrics.json'), 'w');
assert(fid ~= -1, 'Could not write PTO-extension step metrics');
fprintf(fid, '%s\n', jsonencode(metrics));
fclose(fid);
fprintf('RM3 PTO-extension source step audit: %s\n', jsonencode(metrics));
end

function values = baseline_trajectory(baselineDir, label)
prefix = "RM3_PTO_Extension_" + label;
float = readmatrix(fullfile(baselineDir, prefix + "_body1.csv"));
spar = readmatrix(fullfile(baselineDir, prefix + "_body2.csv"));
pto = readmatrix(fullfile(baselineDir, prefix + "_pto1.csv"));
assert(isequal(size(float), [301, 25]) && ...
    isequal(size(spar), [301, 25]) && ...
    isequal(size(pto), [301, 25]), ...
    'The published PTO-extension records are missing');
values = [float(:, 1), float(:, 4), float(:, 10), ...
    spar(:, 4), spar(:, 10), pto(:, 4), pto(:, 10), ...
    pto(:, 16), pto(:, 22)];
assert(max(abs(float(:, 1) - spar(:, 1))) < 1e-10 && ...
    max(abs(float(:, 1) - pto(:, 1))) < 1e-10 && ...
    all(isfinite(values), 'all'), ...
    'The published PTO-extension records are inconsistent');
end

function values = trajectory(output)
float = output.bodies(1);
spar = output.bodies(2);
pto = output.ptos(1);
values = [float.time(:), float.position(:, 3), float.velocity(:, 3), ...
    spar.position(:, 3), spar.velocity(:, 3), ...
    pto.position(:, 3), pto.velocity(:, 3), ...
    pto.forceInternalMechanics(:, 3), ...
    pto.powerInternalMechanics(:, 3)];
assert(all(isfinite(values), 'all') && ...
    max(abs(float.time(:) - spar.time(:))) < 1e-10 && ...
    max(abs(float.time(:) - pto.time(:))) < 1e-10, ...
    'The PTO-extension trajectory is invalid');
end

function write_input(inputFile, contents)
fid = fopen(inputFile, 'w');
assert(fid ~= -1, 'Could not write the disposable source input');
fprintf(fid, '%s', contents);
fclose(fid);
end
