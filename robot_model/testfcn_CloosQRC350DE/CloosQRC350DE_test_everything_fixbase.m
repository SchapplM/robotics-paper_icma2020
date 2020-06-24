% Test script for all fixed-base Matlab/Simulink functions for 
% CloosQRC350DE

% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-10
% (C) Institut für Regelungstechnik, Universität Hannover

clc
clear

CloosQRC350DE_testfunctions_path_init

CloosQRC350DE_varpar_fixbase_kinematics_test
CloosQRC350DE_varpar_fixbase_invdyn_test
CloosQRC350DE_varpar_fixbase_num_test

% Regressorform
CloosQRC350DE_varpar_fixbase_paramlin_test

CloosQRC350DE_compile_test

fprintf('Tests der Fixed-Base Funktionen für CloosQRC350DE abgeschlossen.\n');
