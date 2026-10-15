function zero_ramp_regular_gauges(outputFile)
% Reproduce the compact fixture using pinned WEC-Sim R2025b waveClass.
waves = waveClass('regular');
waves.period = 8;
waves.height = 2.5;
waves.direction = 30;
waves.marker.location = [4, 0; 0, 5; -3, 2];
time = (0:0.1:2)';
waves.setup([0.4, 2], 25, 0, 0.1, 20, time, 9.81, 1000);
assert(isequal(size(waves.waveAmpTime), [21, 2]) && ...
    isequal(size(waves.waveAmpTimeViz), [21, 4]));
writematrix([waves.waveAmpTime, waves.waveAmpTimeViz(:, 2:4)], outputFile);
end
