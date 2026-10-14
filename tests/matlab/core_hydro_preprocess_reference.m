% Compare the pinned MATLAB bodyClass with Python on bundled and current HDF5s.
% Run from the fork root with WEC_SIM_MATLAB_CORE_DIR set to the pinned core.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
coreRoot = getenv('WEC_SIM_MATLAB_CORE_DIR');
if isempty(coreRoot)
    error('Set WEC_SIM_MATLAB_CORE_DIR to the pinned WEC-Sim checkout');
end
if ~isfolder(coreRoot)
    coreRoot = fullfile(repoRoot, coreRoot);
end
addpath(fullfile(coreRoot, 'source', 'objects'));
addpath(fullfile(coreRoot, 'source', 'functions', 'BEMIO'));
outputDir = fullfile(repoRoot, 'matlab-reference-output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

dt = 0.1;
time = (0:dt:10)';
simu = struct('cicTime', time, 'dt', dt, 'rho', 1000, ...
    'gravity', 9.81, 'stateSpace', 0, 'b2b', 0);
for geometry = ["rm3", "oswec", "rm3_current"]
    if geometry == "rm3_current"
        % The current core RM3 HDF5 is distinct from the bundled RM3 file.
        file = fullfile(coreRoot, 'examples', 'RM3', 'hydroData', 'rm3.h5');
    else
        file = fullfile(repoRoot, 'examples', 'data', geometry + ".h5");
    end
    for number = 1:2
        body = bodyClass(file);
        body.setNumber(number);
        body.setDOF(2, 0);
        body.mass = 'equilibrium';
        body.inertia = [1 1 1];
        body.loadHydroData(readBEMIOH5(file, number, 0), 1);

        wave = waveClass('regularCIC');
        wave.period = 8;
        wave.height = 2.5;
        wave.direction = 0;
        hydro = body.hydroData(1);
        wave.setup(hydro.simulation_parameters.w, ...
            hydro.simulation_parameters.waterDepth, 0, dt, ...
            length(time)-1, time, 9.81, 1000);
        body.hydroForcePre(wave, simu, 1);
        force = body.hydroForce.hf1;
        linearRestoring = force.linearHydroRestCoef;
        addedMass = force.fAddedMass;
        radiationKernel = force.irkb;
        excitationReal = force.fExt.re;
        excitationImag = force.fExt.im;
        save(fullfile(outputDir, sprintf('core_hydro_%s_body%d.mat', ...
            geometry, number)), 'linearRestoring', 'addedMass', ...
            'radiationKernel', 'excitationReal', 'excitationImag', '-v7');
    end
end

% OSWEC's multi-heading PM cases also interpolate excitation over both
% frequency and BEM direction. Check that branch on a compact 64-bin sea.
file = fullfile(repoRoot, 'examples', 'data', 'oswec.h5');
for number = 1:2
    body = bodyClass(file);
    body.setNumber(number);
    body.setDOF(2, 0);
    body.mass = 'equilibrium';
    body.inertia = [1 1 1];
    body.loadHydroData(readBEMIOH5(file, number, 0), 1);
    wave = waveClass('irregular');
    wave.period = 8;
    wave.height = 2.5;
    wave.direction = [0 30 90];
    wave.spread = [0.1 0.2 0.7];
    wave.spectrumType = 'PM';
    wave.bem.option = 'Traditional';
    wave.bem.count = 64;
    wave.phaseSeed = 1;
    hydro = body.hydroData(1);
    wave.setup(hydro.simulation_parameters.w, ...
        hydro.simulation_parameters.waterDepth, 0, dt, ...
        length(time)-1, time, 9.81, 1000);
    body.hydroForcePre(wave, simu, 1);
    omega = wave.omega;
    force = body.hydroForce.hf1;
    excitationReal = force.fExt.re;
    excitationImag = force.fExt.im;
    save(fullfile(outputDir, sprintf('core_hydro_oswec_irregular_body%d.mat', ...
        number)), 'omega', 'excitationReal', 'excitationImag', '-v7');
end
