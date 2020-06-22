% This skript creates Table VI
%
% The skript Optimal_Exitation.m must be run at least once to create
% the results in optimal_exitation.mat

fileID1 = fopen('Table_errorModelAonProcess.tex', 'w');

addpath(genpath('../'));
load('optimal_exitation.mat');

fprintf(fileID1, 'joint & \\multicolumn{2}{c||}{Process 1} & \\multicolumn{2}{c|}{Process 2} \\\\ \\hline	 \n');
fprintf(fileID1, '$j$     & $e_j$ (Nm) & $e_j^*$ (\\%%)& $e_j$ (Nm) & $e_j^* (\\%%)$ \\\\ \\hline	 \n');

for i = 1:6
    
    fprintf(fileID1, '%d\t & %1.3f\t & %1.3f\t & %1.3f\t & %1.3f\t \\\\	 \n', i, eP1_A.axes(i), eP1_A.axes_rel(i)*100, eP2_A.axes(i), eP2_A.axes_rel(i)*100);
    
end

fclose('all');

rmpath(genpath('../'));