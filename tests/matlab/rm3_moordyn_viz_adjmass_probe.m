function rm3_moordyn_viz_adjmass_probe
% Change only the source's added-mass shift, on one reproducible 10 s sea.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root, 'matlab-ref', 'source')));
caseDir = fullfile(root, 'applications', 'Paraview_Visualization', 'RM3_MoorDyn_Viz');
outDir = fullfile(root, 'matlab-rm3-moordyn-viz', 'adjmass_probe');
if ~isfolder(outDir), mkdir(outDir); end
for factor = [2, 1]
    target = fullfile(fileparts(caseDir), sprintf('RM3_MoorDyn_Viz_adjmass%d', factor));
    if ~isfolder(target), mkdir(target); end
    for name = {'wecSimInputFile.m', 'RM3MoorDyn.slx', 'userDefinedFunctions.m'}
        copyfile(fullfile(caseDir, name{1}), fullfile(target, name{1}));
    end
    mooringDir = fullfile(target, 'Mooring');
    if ~isfolder(mooringDir), mkdir(mooringDir); end
    copyfile(fullfile(caseDir, 'Mooring', 'lines.txt'), fullfile(mooringDir, 'lines.txt'));
    input = fullfile(target, 'wecSimInputFile.m');
    source = fileread(input);
    assert(contains(source, 'simu.endTime = 80;') && ...
        contains(source, 'simu.dtOut = 0.1;') && ...
        contains(source, 'simu.paraview.option = 1;') && ...
        ~contains(source, 'waves.phaseSeed') && ...
        count(source, 'body(1).inertia = ') == 1 && ...
        count(source, 'body(2).inertia = ') == 1);
    source = strrep(source, 'simu.endTime = 80;', 'simu.endTime = 10;');
    source = strrep(source, 'simu.dtOut = 0.1;', 'simu.dtOut = 0.01;');
    source = strrep(source, 'simu.paraview.option = 1;', 'simu.paraview.option = 0;');
    source = strrep(source, 'waves.period = 8;', ...
        sprintf('waves.period = 8;\nwaves.phaseSeed = 1;'));
    for iBody = 1:2
        marker = sprintf('body(%d).inertia = ', iBody);
        start = strfind(source, marker);
        assert(isscalar(start));
        tail = source(start:end);
        statement = tail(1:find(tail == ';', 1));
        source = strrep(source, statement, sprintf('%s\nbody(%d).adjMassFactor = %d;', ...
            statement, iBody, factor));
    end
    assert(count(source, sprintf('.adjMassFactor = %d;', factor)) == 2);
    fid = fopen(input, 'w');
    assert(fid > 0); cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, source); clear cleaner;
    cd(target);
    wecSim;
    assert(simu.dt == 0.01 && simu.dtOut == 0.01 && ...
        simu.endTime == 10 && waves.phaseSeed == 1 && ...
        body(1).adjMassFactor == factor && body(2).adjMassFactor == factor);
    for iBody = 1:2
        record = output.bodies(iBody);
        writematrix([record.time(:), record.position, record.velocity], ...
            fullfile(outDir, sprintf('factor%d_body%d.csv', factor, iBody)));
    end
    record = output.mooring(1);
    writematrix([record.time(:), record.forceMooring], ...
        fullfile(outDir, sprintf('factor%d_mooring.csv', factor)));
    writematrix(waves.phase(:), fullfile(outDir, sprintf('factor%d_phase.csv', factor)));
    close_system('RM3MoorDyn', 0);
end
end
