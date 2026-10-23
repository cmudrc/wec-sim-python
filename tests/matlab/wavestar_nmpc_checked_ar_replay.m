function wavestar_nmpc_checked_ar_replay(inputPath, applicationsRoot, outputPath)
% Replay the pinned predictor on one checked excitation sequence.
input = load(inputPath);
assert(isequal(size(input.time), [4501, 1]) && ...
    isequal(size(input.moment), [4501, 1]) && ...
    all(isfinite(input.moment)) && ...
    max(abs(input.time - (0:4500)' * 0.05)) < 1e-10, ...
    'The checked WaveStar predictor input changed');
sourceDir = fullfile(applicationsRoot, 'WECCCOMP', ...
    'WECCCOMP_Nonlinear_Model_Predictive');
assert(isfile(fullfile(sourceDir, 'predictor.m')), ...
    'The pinned WaveStar predictor is absent');
addpath(sourceDir);
clear predictor

forecasts = zeros(4501, 40);
for index = 201:4501
    past = input.moment(index-179:index);
    predicted = predictor(past, 18, 40, input.time(index));
    forecasts(index,:) = predicted(:).';
end
assert(all(isfinite(forecasts), 'all'), ...
    'The checked WaveStar predictor produced a nonfinite forecast');
save(outputPath, 'forecasts', '-v7');
end
