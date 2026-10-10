function paraview_mooring_baseline
% Exercise the pinned MoorDyn VTP writer with two lines and time interpolation.
source_times = [0; 1; 2];
frame_times = [0.25; 1.25];
nodes_per_line = [3 4];
moorDyn = struct();
for line = 1:2
    name = sprintf('Line%d', line);
    for node = 0:nodes_per_line(line)-1
        moorDyn.(name).(sprintf('Node%dpx', node)) = ...
            10*line + node + 0.4*source_times;
        moorDyn.(name).(sprintf('Node%dpy', node)) = ...
            2*line - 0.25*node + 0.3*source_times;
        moorDyn.(name).(sprintf('Node%dpz', node)) = ...
            -line - 1.5*node + 0.2*source_times;
    end
    for segment = 1:nodes_per_line(line)-1
        moorDyn.(name).(sprintf('Seg%dTe', segment)) = ...
            100*line + 10*segment + [0; 2.5; 7];
    end
end
directory = fullfile(pwd, 'matlab-paraview-mooring');
mkdir(directory);
mkdir(fullfile(directory, 'mooring1'));
writeParaviewMooring(moorDyn, 'RM3', source_times, 'pinned', 2, ...
    nodes_per_line, directory, frame_times, frame_times, 1);
end
