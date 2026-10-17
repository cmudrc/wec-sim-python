function [startupGap, laterGap] = check_archived_forecast(controllerFile, applicationFolder)
% Check the local numeric AR helpers against an unchanged-source export.
% controllerFile is WECCCOMP_NMPC_SOURCE_controller.mat from the pinned run.
addpath(fileparts(mfilename('fullpath')));
addpath(applicationFolder);
saved = load(controllerFile, 'estimated_states', 'AR_excM_pred');
time = saved.estimated_states.time(:);
state = saved.estimated_states.signals.values;
source = squeeze(saved.AR_excM_pred.signals.values);
assert(numel(time) == 4501 && size(state, 2) == 5 && ...
       isequal(size(source), [40, 4501]), 'Unexpected controller export');
clear predictor
past = zeros(180, 1);
replay = zeros(size(source));
for k = 1:numel(time)
    past = [past(2:end); state(k, end)];
    if time(k) >= 10
        replay(:, k) = predictor(past, 18, 40, time(k));
    end
end
startupGap = max(abs(replay(:, 201:220) - source(:, 201:220)), [], 'all');
laterGap = max(abs(replay(:, 221:end) - source(:, 221:end)), [], 'all');
assert(max(abs(source(:, 201:220)), [], 'all') > 1e6 && ...
       startupGap < 0.1 && laterGap < 1e-9, ...
       'Local AR helper does not reproduce the archived source forecast');
end
