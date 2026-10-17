% Keep sim in a script called directly by the probe function. The pinned
% WEC-Sim runner uses this arrangement so SrcWorkspace='parent' can resolve
% body, waves, and simu from the function that ran initializeWecSim.
sim(model, [], simset('SrcWorkspace', 'parent'));
