function oswec_nonlinear_viz_pressure_baseline
% Export the actual published flap's per-facet nonlinear visualization data.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'OSWEC', 'hydroData');
cd(hydroDir);
if ~isfile('oswec.h5')
    bemio;
end
caseDir = fullfile(repoRoot, 'applications', 'Paraview_Visualization', ...
    'OSWEC_NonLinear_Viz');
inputFile = fullfile(caseDir, 'wecSimInputFile.m');
contents = fileread(inputFile);
oldSetting = 'simu.paraview.option = 1;';
assert(contains(contents, oldSetting), 'Pinned ParaView setting changed');
% Pressure logging is independent of VTP writing. Avoid 1,201 large files.
contents = strrep(contents, oldSetting, 'simu.paraview.option = 0;');
fid = fopen(inputFile, 'w');
assert(fid ~= -1, 'Could not disable VTP writing in disposable checkout');
fprintf(fid, '%s', contents);
fclose(fid);
cd(caseDir);
wecSim;
assert(simu.pressure == 1 && simu.paraview.option == 0 && ...
    simu.dt == 0.1 && simu.endTime == 120 && simu.rampTime == 40 && ...
    strcmp(waves.type, 'regular') && waves.height == 2.5 && ...
    waves.period == 8 && body(1).nonlinearHydro == 2 && ...
    body(1).mass == 127000, 'Published OSWEC pressure case changed');

time = output.bodies(1).cellPressures_time;
pose = output.bodies(1).position;
bodyTime = output.bodies(1).time;
wave = waves.waveAmpTime;
centers = body(1).geometry.center;
cg = body(1).centerGravity;
hydrostatic = output.bodies(1).cellPressures_hydrostatic;
nonlinearWave = output.bodies(1).cellPressures_waveNonLinear;
linearWave = output.bodies(1).cellPressures_waveLinear;
waterDepth = waves.waterDepth;
waveNumber = waves.wavenumber;
deepWater = waves.deepWater;
rho = simu.rho;
gravity = simu.gravity;
numFace = size(centers, 1);
assert(numFace == 1042 && size(hydrostatic, 2) == numFace && ...
    isequal(size(hydrostatic), size(nonlinearWave), size(linearWave)) && ...
    numel(time) == size(hydrostatic, 1) && ...
    all(isfinite(hydrostatic), 'all') && ...
    all(isfinite(nonlinearWave), 'all') && ...
    all(isfinite(linearWave), 'all'), ...
    'Published flap pressure arrays have unexpected dimensions or values');
outDir = fullfile(repoRoot, 'matlab-oswec-nonlinear-viz-pressure');
mkdir(outDir);
save(fullfile(outDir, 'source.mat'), 'time', 'pose', 'bodyTime', ...
    'wave', 'centers', 'cg', 'hydrostatic', 'nonlinearWave', ...
    'linearWave', 'waterDepth', 'waveNumber', 'deepWater', ...
    'rho', 'gravity', '-v7');
end
