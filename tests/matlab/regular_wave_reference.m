% Run from the Python fork root after adding the pinned MATLAB WEC-Sim paths.
% Produces direct waveClass outputs for comparison with the production port.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outputDir = fullfile(repoRoot, 'matlab-reference-output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

times = (0:0.1:200)';
waves = waveClass('regular');
waves.period = 8;
waves.height = 2.5;
waves.direction = 0;
waves.marker.location = [5 5; 10 0; 0 -10];
waves.setup([5.19999512307279, 0.0199999977946844], ...
    'infinite', 100, 0.1, 2000, times, 9.81, 1000);
writematrix(waves.waveAmpTime, fullfile(outputDir, 'regular_origin.csv'));
writematrix(waves.waveAmpTimeViz, fullfile(outputDir, 'regular_markers.csv'));

finite = waveClass('regular');
finite.period = 8;
finite.height = 2.5;
finite.direction = 0;
finite.setup([0.4, 2.0], 25, 0, 0.1, 2000, times, 9.81, 1000);
writematrix([finite.wavenumber, finite.power], ...
    fullfile(outputDir, 'finite_depth.csv'));

% Current IEC PM and JONSWAP spectra with replayable Threefry phases.
for spectrum = ["PM", "JS"]
    for discretization = ["Traditional", "EqualEnergy", "NarrowEqualEnergy"]
        irregular = waveClass('irregular');
        irregular.period = 8;
        irregular.height = 2.5;
        if spectrum == "JS"
            irregular.height = 4;
        end
        irregular.spectrumType = char(spectrum);
        if discretization == "NarrowEqualEnergy"
            irregular.bem.option = 'EqualEnergy';
            irregular.bem.range = [0.5, 1.5];
        else
            irregular.bem.option = char(discretization);
        end
        irregular.bem.count = 64;
        irregular.phaseSeed = 1;
        irregular.marker.location = [5 5; 10 0; 0 -10];
        if spectrum == "PM"
            irregular.direction = [0 30 90];
            irregular.spread = [0.1 0.2 0.7];
        end
        irregular.setup([0.4, 2.0], 'infinite', 1, 0.1, 20, ...
            (0:0.1:2)', 9.81, 1000);
        label = lower(char(spectrum));
        if discretization ~= "Traditional"
            if discretization == "EqualEnergy"
                label = label + "_equal";
            else
                label = label + "_narrow";
            end
            writematrix([irregular.omega, irregular.dOmega, irregular.spectrum], ...
                fullfile(outputDir, sprintf('%s_bins.csv', label)));
        else
            writematrix([irregular.omega, irregular.spectrum], ...
                fullfile(outputDir, sprintf('%s_spectrum.csv', label)));
        end
        writematrix(irregular.power, ...
            fullfile(outputDir, sprintf('%s_power.csv', label)));
        writematrix(irregular.phase, ...
            fullfile(outputDir, sprintf('%s_phase.csv', label)));
        writematrix(irregular.waveAmpTime, ...
            fullfile(outputDir, sprintf('%s_elevation.csv', label)));
        writematrix(irregular.waveAmpTimeViz, ...
            fullfile(outputDir, sprintf('%s_markers.csv', label)));
        if spectrum == "JS"
            writematrix(irregular.gamma, ...
                fullfile(outputDir, sprintf('%s_gamma.csv', label)));
        end
    end
end

% Exercise the two constant branches of WEC-Sim's automatic JS peak factor.
% The existing height=4 comparison covers the interpolated middle branch.
for height = [9, 1]
    irregular = waveClass('irregular');
    irregular.period = 8;
    irregular.height = height;
    irregular.spectrumType = 'JS';
    irregular.bem.option = 'EqualEnergy';
    irregular.bem.range = [0.4, 2.0];
    irregular.bem.count = 64;
    irregular.phaseSeed = 1;
    irregular.setup([0.4, 2.0], 'infinite', 1, 0.1, 20, ...
        (0:0.1:2)', 9.81, 1000);
    label = sprintf('js_auto_height_%g', height);
    writematrix([irregular.omega, irregular.dOmega, irregular.spectrum], ...
        fullfile(outputDir, sprintf('%s_bins.csv', label)));
    writematrix(irregular.phase, ...
        fullfile(outputDir, sprintf('%s_phase.csv', label)));
    writematrix(irregular.waveAmpTime, ...
        fullfile(outputDir, sprintf('%s_elevation.csv', label)));
    writematrix(irregular.gamma, ...
        fullfile(outputDir, sprintf('%s_gamma.csv', label)));
end

fid = fopen(fullfile(outputDir, 'matlab_release.txt'), 'w');
fprintf(fid, '%s\n', version('-release'));
fclose(fid);
