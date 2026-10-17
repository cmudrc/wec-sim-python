% Derived, unmoored RM3 negative control for ode45 MaxStep sensitivity.
% Run from the MATLAB base workspace after setting rm3CoreRoot and
% rm3OutputRoot. The latter must name a directory that does not exist.
assert(exist('rm3CoreRoot', 'var') == 1 && isfolder(rm3CoreRoot));
assert(exist('rm3OutputRoot', 'var') == 1 && ~isfolder(rm3OutputRoot));
rm3CoreRoot = char(rm3CoreRoot);
rm3OutputRoot = char(rm3OutputRoot);
rm3CaseDir = fullfile(rm3CoreRoot, 'examples', 'RM3');
assert(isfile(fullfile(rm3CaseDir, 'RM3.slx')));
mkdir(rm3OutputRoot);
addpath(genpath(fullfile(rm3CoreRoot, 'source')));
rm3OldDir = pwd;
rm3RestoreDir = onCleanup(@() cd(rm3OldDir));

% Change only the source case's wave/solver settings and spar initial pose
% needed to match the published RM3 MoorDyn visualization input. Its RM3.slx
% model is unmoored; this is not the original visualization application.
rm3Rules = {
    '(?m)^simu\.mode[^\r\n]*', 'simu.mode = ''normal'';';
    '(?m)^simu\.explorer[^\r\n]*', 'simu.explorer = ''off'';';
    '(?m)^simu\.rampTime[^\r\n]*', 'simu.rampTime = 0;';
    '(?m)^simu\.endTime[^\r\n]*', 'simu.endTime = 10;';
    '(?m)^simu\.solver[^\r\n]*', 'simu.solver = ''ode45'';';
    '(?m)^simu\.dt =[^\r\n]*', sprintf('simu.dt = 0.01;\nsimu.dtOut = 0.01;\nsimu.cicDt = 0.05;');
    '(?m)^waves = waveClass\(''regular''\)[^\r\n]*', 'waves = waveClass(''irregular'');';
    '(?m)^waves\.height = 2\.5;[^\r\n]*', sprintf(['waves.height = 2;\nwaves.spectrumType = ''JS'';\n' ...
        'waves.bem.option = ''Traditional'';\nwaves.waterDepth = 70;\nwaves.phaseSeed = 1;']);
    '(?m)^body\(2\)\.mass = ''equilibrium'';[^\r\n]*', sprintf(['body(2).mass = ''equilibrium'';\n' ...
        'body(2).initial.displacement = [0 0 -0.21];'])
};

for rm3Label = {'coarse', 'eighth', 'sixteenth', 'fine'}
    rm3Name = rm3Label{1};
    rm3RunDir = fullfile(rm3OutputRoot, rm3Name);
    copyfile(rm3CaseDir, rm3RunDir);
    delete(fullfile(rm3RunDir, 'userDefinedFunctions.m')); % plotting only
    rm3Input = fullfile(rm3RunDir, 'wecSimInputFile.m');
    rm3Text = fileread(rm3Input);
    for rm3Rule = 1:size(rm3Rules, 1)
        assert(numel(regexp(rm3Text, rm3Rules{rm3Rule, 1}, 'match')) == 1, ...
            'Pinned RM3 input changed');
        rm3Text = regexprep(rm3Text, rm3Rules{rm3Rule, 1}, rm3Rules{rm3Rule, 2});
    end
    rm3Handle = fopen(rm3Input, 'w');
    assert(rm3Handle > 0);
    fwrite(rm3Handle, rm3Text);
    fclose(rm3Handle);

    cd(rm3RunDir);
    run('initializeWecSim');
    assert(simu.dt == 0.01 && simu.endTime == 10 && ...
        strcmp(simu.solver, 'ode45') && waves.phaseSeed == 1 && ...
        numel(waves.phase) == 1000);
    [~, rm3Model] = fileparts(simu.simMechanicsFile);
    if strcmp(rm3Name, 'fine')
        rm3MaxStep = 0.0003125;
    elseif strcmp(rm3Name, 'sixteenth')
        rm3MaxStep = 0.000625;
    elseif strcmp(rm3Name, 'eighth')
        rm3MaxStep = 0.00125;
    else
        rm3MaxStep = 0.01;
    end
    set_param(rm3Model, 'MaxStep', num2str(rm3MaxStep, '%.15g'));
    assert(abs(str2double(get_param(rm3Model, 'MaxStep')) - rm3MaxStep) < 1e-12);
    sim(simu.simMechanicsFile);
    run('stopWecSim');
    writematrix(waves.phase(:), fullfile(rm3OutputRoot, [rm3Name '_phase.csv']));
    writematrix([output.wave.time(:), output.wave.elevation(:)], ...
        fullfile(rm3OutputRoot, [rm3Name '_wave.csv']));
    for rm3Body = 1:2
        rm3Record = output.bodies(rm3Body);
        rm3Values = [rm3Record.time(:), rm3Record.position, rm3Record.velocity];
        assert(isequal(size(rm3Values), [1001, 13]) && all(isfinite(rm3Values), 'all'));
        writematrix(rm3Values, fullfile(rm3OutputRoot, ...
            sprintf('%s_body%d.csv', rm3Name, rm3Body)));
    end
    close_system(rm3Model, 0);
end
