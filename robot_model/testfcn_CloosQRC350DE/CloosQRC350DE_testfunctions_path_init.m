% Pfade hinzufügen
% 
% For this file to work, it has to be executed as a complete file, not
% line-wise (because of the function `mfilename`)
% 
% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-03
% (C) Institut für Regelungstechnik, Universität Hannover

%% Pfade hinzufügen

this_path = fileparts( mfilename('fullpath') );
% In diesem Pfad liegen die Ordner matlabfcn, testfcn und simulink,
% die die entsprechenden Dateien der generierten Umgebung des Roboters beinhalten
main_path_robot = clean_absolute_path(fullfile(this_path, '..'));

% Verzeichnis mit generierten Matlab-Funktionen hinzufügen
matlabfcn_path = fullfile(main_path_robot, 'matlabfcn_CloosQRC350DE');
addpath(matlabfcn_path);
