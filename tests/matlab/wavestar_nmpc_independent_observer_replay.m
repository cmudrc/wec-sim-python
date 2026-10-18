function wavestar_nmpc_independent_observer_replay(applicationsRoot)
% Replay the pinned Kalman observer on measurements from Python motion.
% Pass the checkout of pinned WEC-Sim_Applications commit d53d4d4.
fixture = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'fixtures', 'matlab', 'wavestar_nmpc_independent_predictor');
s = load(fullfile(fixture, 'python_observer_input.mat'));
assert(numel(s.time) == 4501 && isequal(size(s.measure), [4501, 2]));

sourceDir = fullfile(applicationsRoot, 'WECCCOMP', ...
    'WECCCOMP_Nonlinear_Model_Predictive');
assert(isfile(fullfile(sourceDir, 'observer_kf.m')));
addpath(sourceDir);
clear observer_kf

rows = zeros(4501, 5);
for k = 1:4501
    rows(k,:) = observer_kf(s.prior(k), s.measure(k,:)', ...
        s.Ad, s.Bd, s.Qkf, s.Rkf)';
end
save(fullfile(fixture, 'matlab_observer.mat'), 'rows', '-v7');
end
