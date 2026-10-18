function wavestar_nmpc_independent_predictor_replay
% Replay MATLAB's published AR/forecast call on Python-generated estimates.
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
end
