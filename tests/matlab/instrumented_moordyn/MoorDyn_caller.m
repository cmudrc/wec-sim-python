function FLines = MoorDyn_caller(X, XD, Time, CoupTime)
% Preserve the pinned MoorDyn call while recording each coupling exchange.
% The audit copies this file into a disposable application directory, so it
% shadows the upstream caller only for that run.
FLines_value = zeros(1, length(X));
FLines_p = libpointer('doublePtr', FLines_value);
calllib('libmoordyn', 'MoorDynStep', X, XD, FLines_p, Time, CoupTime);
FLines = FLines_p.value;

values = [Time, CoupTime, X(:).', XD(:).', FLines(:).'];
assert(numel(values) == 20 && all(isfinite(values)), ...
    'Unexpected RM3 MoorDyn callback data');
logFile = fullfile(fileparts(mfilename('fullpath')), 'moordyn_calls.csv');
fid = fopen(logFile, 'a');
assert(fid > 0, 'Cannot write MoorDyn callback trace');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%.17g,', values(1:end-1));
fprintf(fid, '%.17g\n', values(end));
end
