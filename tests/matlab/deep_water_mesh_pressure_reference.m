function deep_water_mesh_pressure_reference(coreRoot, outputFile)
% Evaluate the pinned regular-wave deep-water FK block on two simple facets.
% This is a prescribed-pressure source check, not a body-motion comparison.
addpath(genpath(fullfile(coreRoot, 'source')));

center = [0, 0, -2; 0, 0, -0.25];
normal = repmat([0, 0, 1], 2, 1);
area = [0.5; 0.5];
times = [0; 4; 0];
heave = [0; 0; 1];
omega = 2*pi/8;
values = zeros(numel(times), 13);

for index = 1:numel(times)
    elevation = 0.5*cos(omega*times(index));
    pose = [0; 0; heave(index); 0; 0; 0];
    [correction, nonlinear, linear] = nonFKForce( ...
        pose, elevation, 0, 1, center, normal, area, ...
        1000, 9.81, [0; 0; 0], 0, omega, 0, 200, 1, 0.1, 10, ...
        times(index), 0);
    values(index, :) = [times(index), heave(index), elevation, ...
        nonlinear(:)', linear(:)', correction(:)'];
end

assert(all(isfinite(values), 'all') && size(values, 2) == 13, ...
    'Deep-water pressure reference is incomplete');
writematrix(values, outputFile);
end
