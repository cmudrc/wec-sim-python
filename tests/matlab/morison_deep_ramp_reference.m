function forces = morison_deep_ramp_reference()
% Reproduce the portable deep-water force values with pinned regWaveMorison.m.
% Add WEC-Sim's source/functions/simulink/model directory to the MATLAB path.
k = 0.062880176966172; % finite-depth dispersion at T=8 s, h=100 m
forces = zeros(2, 6);
rampTime = [2, 0];
for row = 1:2
    forces(row, :) = regWaveMorison( ...
        [0 0 1], [0 0 -2], zeros(1, 6), zeros(1, 6), ...
        [0 0 -1 0 0 0], [0 0 3], [0 0 2], 1, [0 0 0], ...
        0.5, 1000, 0, 100, k, 2*pi/8, 1, rampTime(row), 9.81, ...
        0, 0, 0, 0, 1);
end
end
