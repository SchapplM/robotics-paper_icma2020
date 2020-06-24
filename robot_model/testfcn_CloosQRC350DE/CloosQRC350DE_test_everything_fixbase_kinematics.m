% Test script for fixed-base Matlab/Simulink kinematic-functions for 
% CloosQRC350DE

% Tim-David Job, 2019-11, HiWi bei
% Moritz Schappler, schappler@irt.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

clc
clear

CloosQRC350DE_testfunctions_path_init

CloosQRC350DE_varpar_fixbase_kinematics_test
CloosQRC350DE_compile_test

fprintf('Tests der Fixed-Base Kinematik-Funktionen für CloosQRC350DE abgeschlossen.\n');
