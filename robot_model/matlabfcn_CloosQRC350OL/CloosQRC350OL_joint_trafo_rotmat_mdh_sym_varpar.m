% Calculate homogenous joint transformation matrices for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(6+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = CloosQRC350OL_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:55:31
% EndTime: 2020-06-23 21:55:31
% DurationCPUTime: 0.13s
% Computational Cost: add. (8->8), mult. (0->0), div. (0->0), fcn. (24->12), ass. (0->13)
t61 = cos(qJ(1));
t60 = cos(qJ(2));
t59 = cos(qJ(3));
t58 = cos(qJ(4));
t57 = cos(qJ(5));
t56 = cos(qJ(6));
t55 = sin(qJ(1));
t54 = sin(qJ(2));
t53 = sin(qJ(3));
t52 = sin(qJ(4));
t51 = sin(qJ(5));
t50 = sin(qJ(6));
t1 = [t61, -t55, 0, 0; t55, t61, 0, 0; 0, 0, 1, pkin(1); t54, t60, 0, pkin(2); 0, 0, 1, 0; t60, -t54, 0, 0; t59, -t53, 0, pkin(3); t53, t59, 0, 0; 0, 0, 1, 0; t58, -t52, 0, pkin(4); 0, 0, 1, pkin(5); -t52, -t58, 0, 0; t57, -t51, 0, 0; 0, 0, -1, 0; t51, t57, 0, 0; -t56, t50, 0, 0; 0, 0, 1, pkin(6); t50, t56, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
