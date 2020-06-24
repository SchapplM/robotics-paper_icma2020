% This script creates Table IV
%
% The script Model_Reduction.m must be run for both precesses to create
% the results in ModelReduction_process1.mat and ModelReduction_process2.mat

% Bjoern Volkmann, bjoern.volkmann@imes.uni-hannover.de, 2020-06
% (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover

fileID1 = fopen('Table_ParametersProcess.tex', 'w');

addpath(genpath('../'));
%load results of model reduction
if ~exist('ModelReduction_process1.mat', 'file')
  error('Run Model_Reduction.m for process 1 first!');
end
load('ModelReduction_process1.mat');
theta_P1 = theta_red;

if ~exist('ModelReduction_process2.mat', 'file')
  error('Run Model_Reduction.m for process 2 first!');
end
load('ModelReduction_process2.mat');
theta_P2 = theta_red;

fprintf(fileID1, '      & \\multicolumn{2}{c||}{Process 1} & \\multicolumn{2}{c|}{Process 2}		\\\\ \\hline \n');
fprintf(fileID1, '$i$     & $\\hat{\\theta}_{0,i}$ & $\\hat{\\theta}_{\\text{r},i}$ & $\\hat{\\theta}_{0,i}$ & $\\hat{\\theta}_{\\text{r},i}$ \\\\ \\hline \n');

for vecIndex = 1:35
    
fprintf(fileID1, '%d\t&', vecIndex);
fprintf(fileID1, '\t%1.3f\t&', theta_P1(vecIndex,1));
if theta_P1(vecIndex,end) == 0
    fprintf(fileID1, '\t--\t&');
else
    fprintf(fileID1, '\t%1.3f\t&', theta_P1(vecIndex,end));
end

fprintf(fileID1, '\t%1.3f\t&', theta_P2(vecIndex,1));
if theta_P2(vecIndex,end) == 0
    fprintf(fileID1, '\t--\t');
else
    fprintf(fileID1, '\t%1.3f\t', theta_P2(vecIndex,end));
end
fprintf(fileID1, '\\\\ \n');


end
fclose('all');

rmpath(genpath('../'));
fprintf('Created the file Table_ParametersProcess.tex for Tab. IV.\n');