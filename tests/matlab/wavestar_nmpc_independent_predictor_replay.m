function wavestar_nmpc_independent_predictor_replay(applicationsRoot)
% Replay MATLAB's published predictor on Python-generated estimates.
% Pass the checkout of pinned WEC-Sim_Applications commit d53d4d4.
fixture = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'fixtures', 'matlab', 'wavestar_nmpc_independent_predictor');
input = readmatrix(fullfile(fixture, 'python_moment.csv'));
assert(isequal(size(input), [4501, 2]));
assert(max(abs(input(:,1) - (0:4500)' * 0.05)) < 1e-10);

% The controller first fits at 10 s and refits every one second thereafter.
indices = 201:20:4501;
rows = zeros(numel(indices), 41);
for j = 1:numel(indices)
    index = indices(j);
    past = input(index-179:index, 2);
    model = ar(past, 18);
    predicted = forecast(model, past, 40);
    rows(j,:) = [input(index,1), predicted(:).'];
end
writematrix(rows, fullfile(fixture, 'matlab_forecast.csv'));

sourceDir = fullfile(applicationsRoot, 'WECCCOMP', ...
    'WECCCOMP_Nonlinear_Model_Predictive');
assert(isfile(fullfile(sourceDir, 'predictor.m')));
addpath(sourceDir);
clear predictor

% The NMPC starts at 15 s, but the predictor first fits at 10 s. Call
% every preceding sample so its persistent model has the published history.
activeRows = zeros(4201, 41);
for index = 201:4501
    past = input(index-179:index, 2);
    predicted = predictor(past, 18, 40, input(index,1));
    if index >= 301
        activeRows(index-300,:) = [input(index,1), predicted(:).'];
    end
end
save(fullfile(fixture, 'matlab_active_forecast.mat'), ...
    'activeRows', '-v7');
end
