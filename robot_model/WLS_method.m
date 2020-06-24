% Return inverse dynamics regressor matrix (minimal parameter form)
% Input:
% C
%   Design matrix (stacked regressor matrices)
% y
%   Measurements (motor torques)
% variance
% 
% Output:
% theta
%   Estimation of parameter vector

% Bjoern Volkmann, bjoern.volkmann@imes.uni-hannover.de, 2020-06
% (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover

function theta = WLS_method(C, y, variance)

[N, ~] = size(y);%number of measurement points
[~, numParameter] = size(C);

%rowwise multiplication to reduce complexity since W is a diagonal matrix
w = reshape(variance, [], 1);
w = 1./w;
WC = zeros(size(C));
for row = 1:N
    WC(row,:) = w(row)*C(row,:);
end
CT_C = transpose(C) * WC;% C^T*W*C
CT_y = transpose(C) * (w .* y);% C^T*W*y

rang = rank(CT_C, 1e-8);

%Test for full rank of C^T*W*C
if rang < numParameter
   fprintf('ERROR: The Watrix C^T*W*C is not of full rank \n');
    return
end

%estimate parameters
theta = CT_C \ CT_y;