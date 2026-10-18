function morison_moving_force_probe(sourceRoot, baselineCsv, outDir)
% Re-evaluate archived moving states with the pinned WEC-Sim Morison law.
% The archived option-1 trace anchors the inputs; this focused probe avoids
% rerunning the unrelated 400 s Morison reference-model simulation.
addpath(fullfile(sourceRoot, 'source', 'functions', 'BEMIO'));
addpath(fullfile(sourceRoot, 'source', 'functions', 'simulink', 'model'));
if ~isfolder(outDir)
    mkdir(outDir);
end
rows = readmatrix(baselineCsv);
assert(isequal(size(rows), [8 25]) && all(isfinite(rows), 'all'), ...
    'Archived moving Morison state grid is invalid');
t = rows(:, 1);
pose = rows(:, 2:7);
velocity = rows(:, 8:13);
acceleration = rows(:, 14:19);
point = [.8 -.5 1.2];
cd = [1.1 .7 1.3];
ca = [.8 1.2 .6];
area = [2 3 1.5];
volume = 1.7;
omega = 2*pi/8;
k = calcWaveNumber(omega, 30, 9.81, 0);

% Replaying option 1 first confirms that the archived states still match
% this exact source function before exporting the other source variants.
baseForce = zeros(numel(t), 6);
normalTangentialForce = zeros(numel(t), 6);
for i = 1:numel(t)
    baseForce(i,:) = reshape(regWaveMorison([0 0 1], point, velocity(i,:), ...
        acceleration(i,:), pose(i,:), area, cd, volume, ca, t(i), ...
        1025, 25, 30, k, omega, 1.25, 10, 9.81, 0, 0, 30, 0, 1), 1, 6);
    normalTangentialForce(i,:) = reshape(regWaveMorison([.3 -.7 1.1], ...
        point, velocity(i,:), acceleration(i,:), pose(i,:), area, cd, ...
        volume, ca, t(i), 1025, 25, 30, k, omega, 1.25, 10, 9.81, ...
        0, 0, 30, 0, 2), 1, 6);
end
assert(max(abs(baseForce - rows(:, 20:25)), [], 'all') < 1e-9, ...
    'Pinned option-1 source force changed');
copyfile(baselineCsv, fullfile(outDir, 'MORISON_MOVING_REGULAR_states.csv'));
writematrix([t pose velocity acceleration normalTangentialForce], ...
    fullfile(outDir, 'MORISON_MOVING_NORMAL_TANGENTIAL_states.csv'));

profiles = {'uniform', 'power', 'linear'};
for iProfile = 1:numel(profiles)
    currentForce = zeros(numel(t), 6);
    for i = 1:numel(t)
        currentForce(i,:) = reshape(regWaveMorison([0 0 1], point, ...
            velocity(i,:), acceleration(i,:), pose(i,:), area, cd, ...
            volume, ca, t(i), 1025, 25, 30, k, omega, 1.25, 10, 9.81, ...
            .8, 45, 30, iProfile - 1, 1), 1, 6);
    end
    writematrix([t pose velocity acceleration currentForce], ...
        fullfile(outDir, sprintf( ...
        'MORISON_MOVING_REGULAR_CURRENT_%s_states.csv', profiles{iProfile})));
end
end
