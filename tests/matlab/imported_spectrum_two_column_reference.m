% Generate a seeded sea from frequency/density only, with no saved phases.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
fixtureDir = fullfile(repoRoot, 'tests', 'fixtures', 'matlab', ...
    'imported_spectrum_two_column');
outputDir = fullfile(repoRoot, 'matlab-reference-output');
if ~isfolder(outputDir)
    mkdir(outputDir);
end
input = importdata(fullfile(fixtureDir, 'spectrum.mat'));
assert(isequal(size(input), [7, 2]), ...
    'The two-column imported-spectrum input changed');
bemOmega = h5read(fullfile(repoRoot, 'examples', 'data', 'sphere.h5'), ...
    '/simulation_parameters/w');
times = (0:0.1:20)';
waves = waveClass('spectrumImport');
waves.spectrumFile = fullfile(fixtureDir, 'spectrum.mat');
waves.phaseSeed = 1;
waves.setup(bemOmega, 'infinite', 1, 0.1, 200, times, 9.81, 1000);
assert(numel(waves.omega) == 5, ...
    'BEM filtering did not retain the expected five frequencies');
writematrix([waves.omega, waves.dOmega, waves.spectrum, waves.phase], ...
    fullfile(outputDir, 'imported_two_column_bins_phase.csv'));
writematrix(waves.waveAmpTime, ...
    fullfile(outputDir, 'imported_two_column_elevation.csv'));
