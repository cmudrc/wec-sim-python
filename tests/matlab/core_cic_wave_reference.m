% Compact preprocessing reference for the no-wave and regular CIC paths.
% Run after adding the pinned WEC-Sim source tree to MATLAB's path.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outputDir = getenv('WEC_SIM_MATLAB_CIC_WAVE_DIR');
if isempty(outputDir)
    outputDir = fullfile(repoRoot, 'matlab-reference-output', 'core_cic_wave');
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

bemFrequency = [0.4, 2];
time = (0:0.1:2)';
for waveType = {'noWave', 'noWaveCIC', 'regularCIC'}
    label = waveType{1};
    waves = waveClass(label);
    if strcmp(label, 'noWave')
        waves.period = 8;
    elseif strcmp(label, 'regularCIC')
        waves.period = 8;
        waves.height = 2.5;
        waves.direction = 30;
        waves.marker.location = [4 0; 0 5; -3 2];
    end
    % noWaveCIC leaves period unset: the source selects the minimum BEM
    % frequency for its no-wave reference properties.
    waves.setup(bemFrequency, 25, 1, 0.1, 20, time, 9.81, 1000);
    writematrix([waves.omega, waves.period, waves.wavenumber, ...
                 waves.amplitude], ...
                fullfile(outputDir, [label '_properties.csv']));
    writematrix(waves.waveAmpTime, ...
                fullfile(outputDir, [label '_elevation.csv']));
    if strcmp(label, 'regularCIC')
        writematrix(waves.waveAmpTimeViz, ...
                    fullfile(outputDir, [label '_markers.csv']));
        writematrix(waves.power, ...
                    fullfile(outputDir, [label '_power.csv']));
    end
end
