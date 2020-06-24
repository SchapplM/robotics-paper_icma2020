% Test script for all fixed-base Matlab/Simulink functions for 
% CloosQRC350OL

% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-10
% (C) Institut für Regelungstechnik, Universität Hannover

clc
clear

CloosQRC350OL_testfunctions_path_init

CloosQRC350OL_varpar_fixbase_kinematics_test
CloosQRC350OL_varpar_fixbase_invdyn_test
CloosQRC350OL_varpar_fixbase_num_test

% Regressorform
CloosQRC350OL_varpar_fixbase_paramlin_test

CloosQRC350OL_compile_test

fprintf('Tests der Fixed-Base Funktionen für CloosQRC350OL abgeschlossen.\n');
