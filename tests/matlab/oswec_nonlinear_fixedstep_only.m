function oswec_nonlinear_fixedstep_only
% Refine only Simulink ode4's FixedStep; leave WEC-Sim's grids at 0.1 s.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
addpath(fullfile(repoRoot, 'tests', 'matlab'));
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'OSWEC', 'hydroData');
cd(hydroDir);
if ~isfile('oswec.h5')
    bemio;
end
caseDir = fullfile(repoRoot, 'applications', 'Paraview_Visualization', ...
    'OSWEC_NonLinear_Viz');
inputFile = fullfile(caseDir, 'wecSimInputFile.m');
original = fileread(inputFile);
assert(contains(original, 'simu.dt = 0.1;') && ...
    contains(original, 'simu.pressure = 1;') && ...
    contains(original, 'simu.paraview.option = 1;') && ...
    ~contains(original, 'simu.nonlinearDt') && ...
    ~contains(original, 'simu.cicDt') && ...
    ~contains(original, 'simu.dtOut'), ...
    'Pinned nonlinear OSWEC input changed');
restoreInput = onCleanup(@() writeInput(inputFile, original));
contents = strrep(original, 'simu.pressure = 1;', ...
    'simu.pressure = 0;');
contents = strrep(contents, 'simu.paraview.option = 1;', ...
    'simu.paraview.option = 0;');
writeInput(inputFile, contents);
cd(caseDir);

% Simulink's SrcWorkspace='parent' resolves the calling script's workspace.
% Run the same script-level sequence as the established RM3 solver audit.
oswec_nonlinear_fixedstep_only_simulate;

time = output.bodies(1).time;
pose = output.bodies(1).position;
velocity = output.bodies(1).velocity;
acceleration = output.bodies(1).acceleration;
forceExcitation = output.bodies(1).forceExcitation;
forceRadiationDamping = output.bodies(1).forceRadiationDamping;
forceAddedMass = output.bodies(1).forceAddedMass;
forceRestoring = output.bodies(1).forceRestoring;
assert(numel(time) == 1201 && ...
    all(isfinite([pose, velocity, acceleration, forceExcitation, ...
        forceRadiationDamping, forceAddedMass, forceRestoring]), 'all'), ...
    'Solver-only audit has missing or nonfinite source samples');
outDir = fullfile(repoRoot, 'matlab-oswec-nonlinear-fixedstep-only');
mkdir(outDir);
dt = simu.dt;
dtOut = simu.dtOut;
cicDt = simu.cicDt;
nonlinearDt = simu.nonlinearDt;
save(fullfile(outDir, 'solver-only.mat'), 'time', 'pose', 'velocity', ...
    'acceleration', 'forceExcitation', 'forceRadiationDamping', ...
    'forceAddedMass', 'forceRestoring', 'dt', 'dtOut', 'cicDt', ...
    'nonlinearDt', 'fixedStep', '-v7');
end

function writeInput(path, contents)
fid = fopen(path, 'w');
assert(fid ~= -1, 'Could not write disposable OSWEC input');
fprintf(fid, '%s', contents);
fclose(fid);
end
