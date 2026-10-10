function most_constant_wind_source_baseline(wind_speed, end_time)
% Exercise the pinned MOST platform and above-rated pitch control.
if nargin < 1
    wind_speed = 12;
end
if nargin < 2
    end_time = 10;
end
assert(isscalar(wind_speed) && isfinite(wind_speed) && wind_speed > 0);
assert(isscalar(end_time) && isfinite(end_time) && end_time > 0);
root = pwd;
input_file = fullfile(root, 'applications', 'MOST', 'wecSimInputFile.m');
contents = fileread(input_file);
old_flag = 'wind.constantWindFlag = 0;';
old_speed = 'windSpeed0=8;';
assert(contains(contents, old_flag) && contains(contents, old_speed));
contents = strrep(contents, old_flag, sprintf([ ...
    'wind.constantWindFlag = 1;\n' ...
    'wind.modules = [%.15g %.15g];\n' ...
    'wind.timeBreakpoints = [0 %.15g];\n' ...
    'wind.dt = 0.01;'], wind_speed, wind_speed, end_time));
contents = strrep(contents, old_speed, ...
    sprintf('windSpeed0=%.15g;', wind_speed));
fid = fopen(input_file, 'w');
assert(fid ~= -1);
fprintf(fid, '%s', contents);
fclose(fid);

most_short_source_baseline(end_time);
movefile(fullfile(root, 'matlab-most-short.mat'), ...
    fullfile(root, 'matlab-most-constant-wind.mat'));
end
