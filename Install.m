%% Install

% Find folder of install file
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


disp('MPT toolbox is not automatically installed')

disp('Yalmip and Mosek are not automatically installed')