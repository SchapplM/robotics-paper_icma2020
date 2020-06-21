% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
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
% Tc_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)
% T_c_stack [(6+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = CloosQRC350DE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:19
% EndTime: 2020-06-19 21:37:19
% DurationCPUTime: 0.26s
% Computational Cost: add. (195->61), mult. (214->58), div. (0->0), fcn. (295->12), ass. (0->128)
unknown=NaN(21,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = (pkin(1) + 0);
t4 = sin(qJ(2));
t6 = cos(qJ(2));
t14 = qJ(2) + qJ(3);
t15 = sin(t14);
t16 = t1 * t15;
t17 = cos(t14);
t18 = t1 * t17;
t20 = t4 * pkin(3) + pkin(2);
t21 = t1 * t20;
t23 = t2 * t15;
t24 = t2 * t17;
t25 = t2 * t20;
t27 = t6 * pkin(3);
t29 = cos(qJ(4));
t31 = sin(qJ(4));
t33 = t16 * t29 - t2 * t31;
t36 = -t16 * t31 - t2 * t29;
t37 = t16 * pkin(4);
t38 = t18 * pkin(5);
t39 = t37 + t38 + t21 + 0;
t42 = -t1 * t31 - t23 * t29;
t45 = -t1 * t29 + t23 * t31;
t46 = t23 * pkin(4);
t47 = t24 * pkin(5);
t48 = -t46 - t47 - t25 + 0;
t49 = t17 * t29;
t50 = t17 * t31;
t51 = t17 * pkin(4);
t52 = t15 * pkin(5);
t53 = t51 - t52 + t27 + pkin(1) + 0;
t54 = cos(qJ(5));
t56 = sin(qJ(5));
t58 = t18 * t56 + t33 * t54;
t61 = t18 * t54 - t33 * t56;
t64 = -t24 * t56 + t42 * t54;
t67 = -t24 * t54 - t42 * t56;
t70 = -t15 * t56 + t49 * t54;
t73 = -t15 * t54 - t49 * t56;
t75 = pkin(7) * qJ(5) - qJ(6);
t76 = cos(t75);
t78 = sin(t75);
unknown(1,1) = 1;
unknown(1,2) = 0;
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(2,1) = 0;
unknown(2,2) = 1;
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 1;
unknown(3,4) = 0;
unknown(4,1) = t1;
unknown(4,2) = t2;
unknown(4,3) = 0;
unknown(4,4) = 0;
unknown(5,1) = -t2;
unknown(5,2) = t1;
unknown(5,3) = 0;
unknown(5,4) = 0;
unknown(6,1) = 0;
unknown(6,2) = 0;
unknown(6,3) = 1;
unknown(6,4) = t3;
unknown(7,1) = (t1 * t4);
unknown(7,2) = (t1 * t6);
unknown(7,3) = t2;
unknown(7,4) = (t1 * pkin(2) + 0);
unknown(8,1) = -(t2 * t4);
unknown(8,2) = -(t2 * t6);
unknown(8,3) = t1;
unknown(8,4) = (-t2 * pkin(2) + 0);
unknown(9,1) = t6;
unknown(9,2) = -t4;
unknown(9,3) = 0;
unknown(9,4) = t3;
unknown(10,1) = t16;
unknown(10,2) = t18;
unknown(10,3) = t2;
unknown(10,4) = (t21 + 0);
unknown(11,1) = -t23;
unknown(11,2) = -t24;
unknown(11,3) = t1;
unknown(11,4) = (-t25 + 0);
unknown(12,1) = t17;
unknown(12,2) = -t15;
unknown(12,3) = 0;
unknown(12,4) = (t27 + pkin(1) + 0);
unknown(13,1) = t33;
unknown(13,2) = t36;
unknown(13,3) = t18;
unknown(13,4) = t39;
unknown(14,1) = t42;
unknown(14,2) = t45;
unknown(14,3) = -t24;
unknown(14,4) = t48;
unknown(15,1) = t49;
unknown(15,2) = -t50;
unknown(15,3) = -t15;
unknown(15,4) = t53;
unknown(16,1) = t58;
unknown(16,2) = t61;
unknown(16,3) = -t36;
unknown(16,4) = t39;
unknown(17,1) = t64;
unknown(17,2) = t67;
unknown(17,3) = -t45;
unknown(17,4) = t48;
unknown(18,1) = t70;
unknown(18,2) = t73;
unknown(18,3) = t50;
unknown(18,4) = t53;
unknown(19,1) = (t36 * t78 - t58 * t76);
unknown(19,2) = (-t36 * t76 - t58 * t78);
unknown(19,3) = t61;
unknown(19,4) = (t61 * pkin(6) + t21 + t37 + t38 + 0);
unknown(20,1) = (t45 * t78 - t64 * t76);
unknown(20,2) = (-t45 * t76 - t64 * t78);
unknown(20,3) = t67;
unknown(20,4) = (t67 * pkin(6) - t25 - t46 - t47 + 0);
unknown(21,1) = (-t50 * t78 - t70 * t76);
unknown(21,2) = (t50 * t76 - t70 * t78);
unknown(21,3) = t73;
unknown(21,4) = (t73 * pkin(6) + pkin(1) + t27 + t51 - t52 + 0);
Tc_stack = unknown;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,6+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
