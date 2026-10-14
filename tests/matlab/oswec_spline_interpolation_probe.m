function oswec_spline_interpolation_probe(hydroFile, outputFile)
% Small R2025b check of the spline excitation convention used by OSWEC.
% Frequencies fall between bundled BEM knots; output is time and six forces.
frequency = h5read(hydroFile, '/simulation_parameters/w');
realCoeff = h5read(hydroFile, '/body1/hydro_coeffs/excitation/re');
imagCoeff = h5read(hydroFile, '/body1/hydro_coeffs/excitation/im');
omega = (frequency([21 23 25]) + frequency([22 24 26])) / 2;
headingIndex = [1 4]; % 0 and 30 degrees in the bundled OSWEC HDF5.
spread = [.6 .4];
phase = [.1 .4; .2 .5; .3 .6];
weight = sqrt(.1 * .05 .* spread);
time = (0:.1:.5)';
force = zeros(numel(time), 6);
for dof = 1:6
    for direction = 1:2
        re = interp1(frequency, realCoeff(:, headingIndex(direction), dof), ...
                     omega, 'spline') * 1000 * 9.81;
        im = interp1(frequency, imagCoeff(:, headingIndex(direction), dof), ...
                     omega, 'spline') * 1000 * 9.81;
        for index = 1:3
            angle = time * omega(index) + phase(index, direction);
            force(:, dof) = force(:, dof) + weight(direction) * (...
                re(index) * cos(angle) - im(index) * sin(angle));
        end
    end
end
ramp = ones(size(time));
early = time < .2;
ramp(early) = (1 - cos(pi * time(early) / .2)) / 2;
writematrix([time, force .* ramp], outputFile);
end
