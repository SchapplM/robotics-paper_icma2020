% Convert vector of kinematic parameters to modified DH parameters of
% CloosQRC350DE
%
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [6x1]
%   Rotation around z
% b_mdh [6x1]
%   Translation along z
% alpha_mdh [6x1]
%   Rotation around x
% a_mdh [6x1]
%   Translation along x
% theta_mdh [6x1]
%   Rotation around z
% d_mdh [6x1]
%   Translation along z
% qoffset_mdh [6x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = CloosQRC350DE_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(7,1)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_pkin2mdhparam: Kinematic parameters pkin have to be [7x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = -pi / 0.2e1;
t2 = [0; t1; 0; t1; pi / 0.2e1; t1;];
alpha_mdh = t2;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(2); pkin(3); pkin(4); 0; 0;];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(1); 0; 0; pkin(5); 0; pkin(6);];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; -pi / 0.2e1; 0; 0; 0; pi;];
qoffset_mdh = t1;
