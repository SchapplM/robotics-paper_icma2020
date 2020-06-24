% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = CloosQRC350OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:55:31
% EndTime: 2020-06-23 21:55:31
% DurationCPUTime: 0.26s
% Computational Cost: add. (180->42), mult. (202->57), div. (0->0), fcn. (295->12), ass. (0->42)
t28 = qJ(2) + qJ(3);
t25 = cos(t28);
t31 = sin(qJ(4));
t49 = t25 * t31;
t36 = cos(qJ(4));
t48 = t25 * t36;
t24 = sin(t28);
t33 = sin(qJ(1));
t47 = t33 * t24;
t20 = t33 * t25;
t46 = t33 * t31;
t45 = t33 * t36;
t38 = cos(qJ(1));
t44 = t38 * t24;
t21 = t38 * t25;
t43 = t38 * t31;
t42 = t38 * t36;
t27 = pkin(1) + 0;
t32 = sin(qJ(2));
t23 = t32 * pkin(3) + pkin(2);
t41 = t33 * t23 + 0;
t40 = t38 * t23 + 0;
t37 = cos(qJ(2));
t39 = t37 * pkin(3) + t27;
t5 = pkin(4) * t47 + pkin(5) * t20 + t41;
t6 = pkin(4) * t44 + pkin(5) * t21 + t40;
t7 = t25 * pkin(4) - t24 * pkin(5) + t39;
t35 = cos(qJ(5));
t34 = cos(qJ(6));
t30 = sin(qJ(5));
t29 = sin(qJ(6));
t13 = t24 * t42 + t46;
t12 = t24 * t43 - t45;
t11 = t24 * t45 - t43;
t10 = t24 * t46 + t42;
t9 = -t24 * t30 + t35 * t48;
t8 = -t24 * t35 - t30 * t48;
t4 = t13 * t35 + t30 * t21;
t3 = -t13 * t30 + t35 * t21;
t2 = t11 * t35 + t30 * t20;
t1 = -t11 * t30 + t35 * t20;
t14 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t38, -t33, 0, 0; t33, t38, 0, 0; 0, 0, 1, t27; t38 * t32, t38 * t37, -t33, t38 * pkin(2) + 0; t33 * t32, t33 * t37, t38, t33 * pkin(2) + 0; t37, -t32, 0, t27; t44, t21, -t33, t40; t47, t20, t38, t41; t25, -t24, 0, t39; t13, -t12, t21, t6; t11, -t10, t20, t5; t48, -t49, -t24, t7; t4, t3, t12, t6; t2, t1, t10, t5; t9, t8, t49, t7; t12 * t29 - t4 * t34, t12 * t34 + t4 * t29, t3, t3 * pkin(6) + t6; t10 * t29 - t2 * t34, t10 * t34 + t2 * t29, t1, t1 * pkin(6) + t5; t29 * t49 - t9 * t34, t9 * t29 + t34 * t49, t8, t8 * pkin(6) + t7;];
Tc_stack = t14;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,6+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
