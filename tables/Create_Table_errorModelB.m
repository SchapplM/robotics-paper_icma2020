% This skript creates Table III
%
% The skript Optimal_Exitation.m must be run at least once to create
% the results in optimal_exitation.mat

fileID1 = fopen('Table_errorModelB.tex', 'w');

addpath(genpath('../'));
load('optimal_exitation.mat');


fprintf(fileID1, 'joint & \\multicolumn{2}{c||}{Tracectory B} & \\multicolumn{2}{c|}{Tracectory C} \\\\ \\hline	 \n');
fprintf(fileID1, '$j$     & $e_j$ (Nm) & $e_j^*$ (\\%%)& $e_j$ (Nm) & $e_j^* (\\%%)$ \\\\ \\hline	 \n');

for i = 1:6
    
    fprintf(fileID1, '%d\t & %1.3f\t & %1.3f\t & %1.3f\t & %1.3f\t \\\\	 \n', i, eB_B.axes(i), eB_B.axes_rel(i)*100, eC_B.axes(i), eC_B.axes_rel(i)*100);
    
end

fclose('all');

rmpath(genpath('../'));