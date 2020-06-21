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
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
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
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:19
% EndTime: 2020-06-19 21:37:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (14->11), mult. (4->1), div. (0->0), fcn. (24->12), ass. (0->85)
unknown=NaN(18,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = sin(qJ(2));
t4 = cos(qJ(2));
t5 = cos(qJ(3));
t6 = sin(qJ(3));
t7 = cos(qJ(4));
t8 = sin(qJ(4));
t9 = cos(qJ(5));
t10 = sin(qJ(5));
t12 = pkin(7) * qJ(5) - qJ(6);
t13 = cos(t12);
t14 = sin(t12);
unknown(1,1) = t1;
unknown(1,2) = t2;
unknown(1,3) = 0.0e0;
unknown(1,4) = 0.0e0;
unknown(2,1) = -t2;
unknown(2,2) = t1;
unknown(2,3) = 0.0e0;
unknown(2,4) = 0.0e0;
unknown(3,1) = 0.0e0;
unknown(3,2) = 0.0e0;
unknown(3,3) = 0.1e1;
unknown(3,4) = pkin(1);
unknown(4,1) = t3;
unknown(4,2) = t4;
unknown(4,3) = 0.0e0;
unknown(4,4) = pkin(2);
unknown(5,1) = 0.0e0;
unknown(5,2) = 0.0e0;
unknown(5,3) = 0.1e1;
unknown(5,4) = 0.0e0;
unknown(6,1) = t4;
unknown(6,2) = -t3;
unknown(6,3) = 0.0e0;
unknown(6,4) = 0.0e0;
unknown(7,1) = t5;
unknown(7,2) = -t6;
unknown(7,3) = 0.0e0;
unknown(7,4) = pkin(3);
unknown(8,1) = t6;
unknown(8,2) = t5;
unknown(8,3) = 0.0e0;
unknown(8,4) = 0.0e0;
unknown(9,1) = 0.0e0;
unknown(9,2) = 0.0e0;
unknown(9,3) = 0.1e1;
unknown(9,4) = 0.0e0;
unknown(10,1) = t7;
unknown(10,2) = -t8;
unknown(10,3) = 0.0e0;
unknown(10,4) = pkin(4);
unknown(11,1) = 0.0e0;
unknown(11,2) = 0.0e0;
unknown(11,3) = 0.1e1;
unknown(11,4) = pkin(5);
unknown(12,1) = -t8;
unknown(12,2) = -t7;
unknown(12,3) = 0.0e0;
unknown(12,4) = 0.0e0;
unknown(13,1) = t9;
unknown(13,2) = -t10;
unknown(13,3) = 0.0e0;
unknown(13,4) = 0.0e0;
unknown(14,1) = 0.0e0;
unknown(14,2) = 0.0e0;
unknown(14,3) = -0.1e1;
unknown(14,4) = 0.0e0;
unknown(15,1) = t10;
unknown(15,2) = t9;
unknown(15,3) = 0.0e0;
unknown(15,4) = 0.0e0;
unknown(16,1) = -t13;
unknown(16,2) = -t14;
unknown(16,3) = 0.0e0;
unknown(16,4) = 0.0e0;
unknown(17,1) = 0.0e0;
unknown(17,2) = 0.0e0;
unknown(17,3) = 0.1e1;
unknown(17,4) = pkin(6);
unknown(18,1) = -t14;
unknown(18,2) = t13;
unknown(18,3) = 0.0e0;
unknown(18,4) = 0.0e0;
T_stack = unknown;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
