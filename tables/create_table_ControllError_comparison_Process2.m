% This skript creates an equivalent of Table VII for process 2
%
% This skript creates a latex table to compare all models applied to
% process 2 in the experimental validation

addpath(genpath('../'));
clear
close all
clc

%% Load Data

% Data for model A
load('Process 2/Model_A.mat');
error_A = rms(abs(deg_e));


% Data for model B
load('Process 2/Model_B.mat');
error_B = rms(abs(deg_e));


% Data for process model
load('Process 2/Process_Model.mat');
error_Pr = rms(abs(deg_e));


% Data for process model
load('Process 2/Without_Model.mat');
error_0 = rms(abs(deg_e));


%% Write File

fileID1 = fopen('Table_RefErrorProcess2.tex', 'w');

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