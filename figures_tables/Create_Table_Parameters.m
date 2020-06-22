% This skript creates Table I
%
% The skript Optimal_Exitation.m must be run at least once to create
% the results in optimal_exitation.mat

fileID1 = fopen('Table_Parameters.tex', 'w');

addpath(genpath('../'));
load('optimal_exitation.mat');

fprintf(fileID1, '		& \\multicolumn{2}{c||}{Model A}		& \\multicolumn{2}{c|}{Model B}		\\\\ \\hline \n');
fprintf(fileID1, '$i$ & $\\hat{\\theta}_{i}$   & $\\sigma_{i}^*$ (\\%%) & $\\hat{\\theta}_{i}$   & $\\sigma_{i}^*$ (\\%%) \\\\ \\hline \n');

for vecIndex = 1:35
    
fprintf(fileID1, '%d\t& \t%1.3f\t& \t%1.3f\t& \t%1.3f\t& \t%1.3f\t \\\\ \n', vecIndex, theta_A(vecIndex), sigma_A(vecIndex), theta_B(vecIndex), sigma_B(vecIndex));

end
fclose('all');

rmpath(genpath('../'));