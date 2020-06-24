% This script creates Table I
%
% The script Optimal_Exitation.m must be run at least once to create
% the results in optimal_exitation.mat

% Bjoern Volkmann, bjoern.volkmann@imes.uni-hannover.de, 2020-06
% (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover

fileID1 = fopen('Table_Parameters.tex', 'w');

addpath(genpath('../'));
if ~exist('optimal_exitation.mat', 'file')
  error('Run Optimal_Exitation.m first!');
end
load('optimal_exitation.mat');
fprintf(fileID1, '		& \\multicolumn{2}{c||}{Model A}		& \\multicolumn{2}{c|}{Model B}		\\\\ \\hline \n');
fprintf(fileID1, '$i$ & $\\hat{\\theta}_{i}$   & $\\sigma_{i}^*$ (\\%%) & $\\hat{\\theta}_{i}$   & $\\sigma_{i}^*$ (\\%%) \\\\ \\hline \n');

for vecIndex = 1:35   
fprintf(fileID1, '%d\t& \t%1.3f\t& \t%1.3f\t& \t%1.3f\t& \t%1.3f\t \\\\ \n', vecIndex, theta_A(vecIndex), sigma_A(vecIndex), theta_B(vecIndex), sigma_B(vecIndex));
end
fclose('all');

rmpath(genpath('../'));
fprintf('Created the file Table_Parameters.tex for Tab. I.\n');