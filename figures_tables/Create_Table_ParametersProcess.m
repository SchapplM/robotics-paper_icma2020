% This skript creates Table IV
%
% The skript Model_Reduction.m must be run for both precesses to create
% the results in ModelReduction_process1.mat and ModelReduction_process2.mat

fileID1 = fopen('Table_ParametersProcess.tex', 'w');

addpath(genpath('../'));
%load results of model reduction
load('ModelReduction_process1.mat');
theta_P1 = theta_red;

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