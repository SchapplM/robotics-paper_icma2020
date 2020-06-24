% This script creates Table V
%
% The script Model_Reduction.m must be run for both precesses to create
% the results in ModelReduction_process1.mat and ModelReduction_process2.mat

% Bjoern Volkmann, bjoern.volkmann@imes.uni-hannover.de, 2020-06
% (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover

addpath(genpath('../'));
clear
fileID1 = fopen('Table_errorProcessOnProcess.tex', 'w');

%load results of model reduction
if ~exist('ModelReduction_process1.mat', 'file')
  error('Run Model_Reduction.m for process 1 first!');
end
load('ModelReduction_process1.mat');
eP1 = eV;

if ~exist('ModelReduction_process2.mat', 'file')
  error('Run Model_Reduction.m for process 2 first!');
end
load('ModelReduction_process2.mat');
eP2 = eV;

fprintf(fileID1, 'joint & \\multicolumn{2}{c||}{Process 1} & \\multicolumn{2}{c|}{Process 2} \\\\ \\hline	 \n');
fprintf(fileID1, '$j$     & $e_j$ (Nm) & $e_j^*$ (\\%%)& $e_j$ (Nm) & $e_j^* (\\%%)$ \\\\ \\hline	 \n');

for i = 1:6
    fprintf(fileID1, '%d\t & %1.3f\t & %1.3f\t & %1.3f\t & %1.3f\t \\\\	 \n', i, eP1.result(i), eP1.result_rel(i)*100, eP2.result(i), eP2.result_rel(i)*100);
end
fclose('all');
rmpath(genpath('../'));
fprintf('Created the file Table_errorProcessOnProcess.tex for Tab. V.\n');