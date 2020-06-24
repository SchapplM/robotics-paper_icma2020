% Calculate homogenous joint transformation matrices for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(6+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = CloosQRC350DE_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 20:58:47
% EndTime: 2020-06-23 20:58:47
% DurationCPUTime: 0.07s
% Computational Cost: add. (14->9), mult. (4->1), div. (0->0), fcn. (24->12), ass. (0->14)
t61 = cos(qJ(1));
t60 = cos(qJ(2));
t59 = cos(qJ(3));
t58 = cos(qJ(4));
t57 = cos(qJ(5));
t56 = sin(qJ(1));
t55 = sin(qJ(2));
t54 = sin(qJ(3));
t53 = sin(qJ(4));
t52 = sin(qJ(5));
t51 = -pkin(7) * qJ(5) + qJ(6);
t50 = cos(t51);
t49 = sin(t51);
t1 = [t61, t56, 0, 0; -t56, t61, 0, 0; 0, 0, 1, pkin(1); t55, t60, 0, pkin(2); 0, 0, 1, 0; t60, -t55, 0, 0; t59, -t54, 0, pkin(3); t54, t59, 0, 0; 0, 0, 1, 0; t58, -t53, 0, pkin(4); 0, 0, 1, pkin(5); -t53, -t58, 0, 0; t57, -t52, 0, 0; 0, 0, -1, 0; t52, t57, 0, 0; -t50, t49, 0, 0; 0, 0, 1, pkin(6); t49, t50, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
