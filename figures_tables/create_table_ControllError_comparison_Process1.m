% This script creates Table VII
%
% This script creates a latex table to compare all models applied to
% process 1 in the experimental validation

% Bjoern Volkmann, bjoern.volkmann@imes.uni-hannover.de, 2020-06
% (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover

addpath(genpath('../'));
clear
close all
clc

%% Load Data
% Data for model A
load('Process 1/Model_A.mat');
error_A = rms(abs(deg_e));
% Data for model B
load('Process 1/Model_B.mat');
error_B = rms(abs(deg_e));
% Data for process model
load('Process 1/Process_Model.mat');
error_Pr = rms(abs(deg_e));
% Data for process model
load('Process 1/Without_Model.mat');
error_0 = rms(abs(deg_e));

%% Write File
fileID1 = fopen('Table_RefErrorProcess1.tex', 'w');
%header
fprintf(fileID1, 'joint\t & no Model A\t  & Model A\t & Model B\t& Process Model\t \\\\ \\hline \n');
fprintf(fileID1, 'j\t & $q_{e,RMS,j}$ ($^{\\circ}$)\t  & $q_{e,RMS,j}$ ($^{\\circ}$)\t & $q_{e,RMS,j}$ ($^{\\circ}$)\t & $q_{e,RMS,j}$ ($^{\\circ}$)\t \\\\ \\hline \n');
%entries
for i = 1:6
    fprintf(fileID1, '%d\t& %1.4f\t& %1.4f\t& %1.4f\t& %1.4f\t \\\\ \\hline \n', i, error_0(i), error_A(i), error_B(i), error_Pr(i));
end
%close file
fclose('all');
rmpath(genpath('./'));
fprintf('Created the file Table_RefErrorProcess1.tex for Tab. VII.\n');