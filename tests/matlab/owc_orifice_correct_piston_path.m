% Patch only a disposable copy of the pinned OWC OrificeModel Simulink model.
% In MATLAB's base workspace, set coreRoot and caseDir before running this.
assert(exist('coreRoot', 'var') == 1 && exist('caseDir', 'var') == 1);
assert(~contains(caseDir, 'wec-sim-applications-reference'), ...
    'Do not edit the pinned Applications checkout');
addpath(genpath(fullfile(coreRoot, 'source')));
modelFile = fullfile(caseDir, 'OWC_GBM.slx');
assert(isfile(modelFile));
load_system(modelFile);

sub = 'OWC_GBM/Flex Body/FlexibleBody';
product = [sub '/Product4'];
ports = get_param(product, 'PortHandles');
sourceLine = get_param(ports.Inport(2), 'Line');
sourceBlock = get_param(sourceLine, 'SrcBlockHandle');
assert(strcmp(getfullname(sourceBlock), [sub '/Excitation']));
assert(strcmp(get_param([sub '/Add'], 'Inputs'), '--++'));
assert(strcmp(get_param([sub '/Gain'], 'Gain'), '1'));
assert(strcmp(get_param([sub '/GBM'], 'BlockType'), 'StateSpace'));

% The source reports the piston contribution in Add but leaves it out of
% Product4, which drives the GBM state-space force input. That state-space
% already contains flexible stiffness and damping; add only the piston force.
delete_line(sub, 'Excitation/1', 'Product4/2');
newBlock = [sub '/Piston input to GBM'];
position = get_param(product, 'Position');
add_block('simulink/Math Operations/Sum', newBlock, ...
    'Inputs', '++', ...
    'Position', [position(1)-110, position(2)+55, ...
                 position(1)-75, position(2)+90]);
add_line(sub, 'Excitation/1', 'Piston input to GBM/1', 'autorouting', 'on');
add_line(sub, 'Gain/1', 'Piston input to GBM/2', 'autorouting', 'on');
add_line(sub, 'Piston input to GBM/1', 'Product4/2', 'autorouting', 'on');
ports = get_param(product, 'PortHandles');
newLine = get_param(ports.Inport(2), 'Line');
assert(strcmp(getfullname(get_param(newLine, 'SrcBlockHandle')), newBlock));
save_system('OWC_GBM', modelFile);
close_system('OWC_GBM', 0);
