function paraview_wave_baseline
% Export two small regular-wave grids with the pinned WEC-Sim VTP writer.
waves = waveClass('regular');
waves.height = 2.5;
waves.period = 8;
waves.direction = 90;
times = [2 3];
waves.setup([0.3 1.5], 30, 0, 1, length(times), times, 9.81, 1000);
output = fullfile(pwd, 'matlab-paraview-wave');
mkdir(output);
mkdir(fullfile(output, 'waves'));
writeParaviewWave(waves, times, 3, 2, 6, 'OSWEC', 'pinned', ...
    0, output, times, 9.81);
writematrix([waves.amplitude waves.wavenumber waves.omega ...
    waves.direction waves.waterDepth], fullfile(output, 'parameters.csv'));

% A two-frequency imported record exercises WEC-Sim's one-sided FFT surface
% reconstruction, including phase and the deep-water wavenumber from g.
imported = waveClass('elevationImport');
imported.direction = 60;
imported.waterDepth = 30;
recordTime = (0:15)' * 0.25;
recordElevation = 0.7*cos(pi*recordTime + 0.35) ...
    + 0.2*cos(2*pi*recordTime - 0.8);
importedOutput = fullfile(output, 'imported');
mkdir(importedOutput);
mkdir(fullfile(importedOutput, 'waves'));
recordFile = fullfile(importedOutput, 'record.csv');
writematrix([recordTime recordElevation], recordFile);
imported.elevationFile = recordFile;
imported.setup([0.3 1.5], 30, 0, 0.25, 15, recordTime, 9.7, 1000);
writeParaviewWave(imported, 0.375, 3, 2, 2, 'Imported', 'pinned', ...
    0, importedOutput, 0.375, 9.7);
writematrix(imported.waveAmpTime, fullfile(importedOutput, 'record.csv'));
end
