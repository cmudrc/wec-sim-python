function rm3_moordyn_solver_coupling_audit
% Separate ode45 MaxStep from MoorDyn's simu.dt-driven coupling pulse.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot, 'matlab-ref', 'source')));
hydroDir = fullfile(repoRoot, 'applications', '_Common_Input_Files', ...
    'RM3', 'hydroData');
cd(hydroDir);
if ~isfile('rm3.h5')
    bemio;
end
sourceDir = fullfile(repoRoot, 'applications', ...
    'Paraview_Visualization', 'RM3_MoorDyn_Viz');
outDir = fullfile(repoRoot, 'matlab-rm3-solver-coupling-audit');
mkdir(outDir);
runOne(sourceDir, outDir, 'published_step', 0.01, 0.01);
runOne(sourceDir, outDir, 'solver_only', 0.01, 0.0003125);
runOne(sourceDir, outDir, 'solver_and_moordyn', 0.0003125, 0.0003125);
end

function runOne(sourceDir, outDir, label, couplingStep, maximumStep)
        caseDir = fullfile(fileparts(sourceDir), ...
            ['RM3_MoorDyn_Viz_solver_audit_' label]);
        mkdir(caseDir);
        for filename = {'wecSimInputFile.m', 'RM3MoorDyn.slx', ...
                'userDefinedFunctions.m'}
            copyfile(fullfile(sourceDir, filename{1}), ...
                fullfile(caseDir, filename{1}));
        end
        mkdir(fullfile(caseDir, 'Mooring'));
        copyfile(fullfile(sourceDir, 'Mooring', 'lines.txt'), ...
            fullfile(caseDir, 'Mooring', 'lines.txt'));
        inputFile = fullfile(caseDir, 'wecSimInputFile.m');
        contents = fileread(inputFile);
        changes = {
            'simu.endTime = 80;', 'simu.endTime = 10;';
            'simu.dt = 0.01;', sprintf('simu.dt = %.7f;', couplingStep);
            'simu.dtOut = 0.1;', 'simu.dtOut = 0.01;';
            'simu.paraview.option = 1;', 'simu.paraview.option = 0;';
            'waves.period = 8;', sprintf('waves.period = 8;\nwaves.phaseSeed = 1;')
        };
        for i = 1:size(changes, 1)
            assert(contains(contents, changes{i, 1}), ...
                'The pinned RM3 visualization input changed');
            contents = strrep(contents, changes{i, 1}, changes{i, 2});
        end
        fid = fopen(inputFile, 'w');
        assert(fid > 0);
        cleanup = onCleanup(@() fclose(fid));
        fwrite(fid, contents);
        clear cleanup;
        cd(caseDir);

        % These are the three calls made by the pinned wecSim.m. Override
        % MaxStep after setup, leaving simu.dt, its wave grid, and the
        % MoorDynCaller pulse/step at the requested couplingStep.
        run('initializeWecSim');
        assert(abs(simu.dt - couplingStep) < 1e-12 && ...
            simu.dtOut == 0.01 && simu.endTime == 10 && ...
            strcmp(simu.solver, 'ode45') && waves.phaseSeed == 1 && ...
            mooring(1).moorDyn == 1, 'Audit settings changed');
        [~, modelName, ~] = fileparts(simu.simMechanicsFile);
        set_param(modelName, 'MaxStep', ...
            num2str(maximumStep, '%.15g'));
        assert(abs(str2double(get_param(modelName, ...
            'MaxStep')) - maximumStep) < 1e-12, ...
            'Could not set the independent ode45 maximum step');
        sim(simu.simMechanicsFile, [], simset('SrcWorkspace', 'parent'));
        run('stopWecSim');
        assert(size(waves.phase, 1) == 1000, 'Unexpected phase grid');
        writematrix(waves.phase(:), ...
            fullfile(outDir, [label '_phase.csv']));
        writematrix([output.wave.time(:), output.wave.elevation(:)], ...
            fullfile(outDir, [label '_wave.csv']));
        for iBody = 1:2
            record = output.bodies(iBody);
            values = [record.time(:), record.position, record.velocity];
            assert(isequal(size(values), [1001, 13]) && ...
                all(isfinite(values), 'all'), ...
                'Incomplete body trace');
            writematrix(values, fullfile(outDir, ...
                sprintf('%s_body%d.csv', label, iBody)));
        end
        record = output.mooring(1);
        values = [record.time(:), record.position, ...
            record.velocity, record.forceMooring];
        assert(isequal(size(values), [1001, 19]) && ...
            all(isfinite(values), 'all'), ...
            'Incomplete MoorDyn trace');
        writematrix(values, fullfile(outDir, ...
            [label '_mooring.csv']));
        close_system('RM3MoorDyn', 0);
end
