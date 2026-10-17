% Regenerate the compact RM3 PM sea at two heights from pinned R2025b source.
assert(strcmp(version('-release'), '2025b'), ...
    'The compact RM3 scaling reference was generated with R2025b');
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
sourceDir = getenv('WEC_SIM_SOURCE_DIR');
if ~isempty(sourceDir)
    addpath(genpath(sourceDir));
end
outputDir = getenv('WEC_SIM_MATLAB_RM3_SCALING_DIR');
if isempty(outputDir)
    outputDir = fullfile(repoRoot, 'matlab-reference-output', ...
                         'rm3_small_wave_scaling');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

hydroFile = fullfile(repoRoot, 'examples', 'data', 'rm3.h5');
frequency = h5read(hydroFile, '/simulation_parameters/w');
bemRange = [min(frequency(:)), max(frequency(:))];
times = (0:0.1:20)';
for index = 1:2
    waves = waveClass('irregular');
    waves.period = 8;
    waves.height = 0.01 * index;
    waves.spectrumType = 'PM';
    waves.bem.option = 'EqualEnergy';
    waves.bem.count = 32;
    waves.phaseSeed = 4;
    waves.setup(bemRange, 'infinite', 5, 0.1, 200, times, 9.81, 1000);
    writematrix(waves.waveAmpTime, ...
                fullfile(outputDir, sprintf('wave_%03d.csv', index)));
    writematrix([waves.omega, waves.dOmega, waves.spectrum], ...
                fullfile(outputDir, sprintf('bins_%03d.csv', index)));
    if index == 1
        writematrix(waves.phase, fullfile(outputDir, 'phase_seed4.csv'));
    end
end
