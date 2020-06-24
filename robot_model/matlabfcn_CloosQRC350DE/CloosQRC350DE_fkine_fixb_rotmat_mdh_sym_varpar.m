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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
% OptimizationMode: 2
% StartTime: 2020-06-23 20:58:47
% EndTime: 2020-06-23 20:58:47
% DurationCPUTime: 0.30s
% Computational Cost: add. (195->46), mult. (214->59), div. (0->0), fcn. (295->12), ass. (0->42)
t29 = qJ(2) + qJ(3);
t26 = cos(t29);
t31 = sin(qJ(4));
t48 = t26 * t31;
t35 = cos(qJ(4));
t47 = t26 * t35;
t33 = sin(qJ(1));
t46 = t33 * t26;
t45 = t33 * t31;
t44 = t33 * t35;
t25 = sin(t29);
t37 = cos(qJ(1));
t43 = t37 * t25;
t18 = t37 * t26;
t42 = t37 * t31;
t41 = t37 * t35;
t28 = pkin(1) + 0;
t32 = sin(qJ(2));
t23 = t32 * pkin(3) + pkin(2);
t40 = t37 * t23 + 0;
t39 = -t23 * t33 + 0;
t36 = cos(qJ(2));
t38 = t36 * pkin(3) + t28;
t6 = pkin(4) * t43 + pkin(5) * t18 + t40;
t7 = t26 * pkin(4) - t25 * pkin(5) + t38;
t5 = (-pkin(4) * t25 - pkin(5) * t26) * t33 + t39;
t34 = cos(qJ(5));
t30 = sin(qJ(5));
t24 = pkin(7) * qJ(5) - qJ(6);
t20 = cos(t24);
t19 = sin(t24);
t13 = t25 * t41 - t45;
t12 = t25 * t42 + t44;
t11 = -t25 * t44 - t42;
t10 = t25 * t45 - t41;
t9 = -t25 * t30 + t34 * t47;
t8 = -t25 * t34 - t30 * t47;
t4 = t13 * t34 + t30 * t18;
t3 = -t13 * t30 + t34 * t18;
t2 = t11 * t34 - t30 * t46;
t1 = -t11 * t30 - t34 * t46;
t14 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t37, t33, 0, 0; -t33, t37, 0, 0; 0, 0, 1, t28; t37 * t32, t37 * t36, t33, t37 * pkin(2) + 0; -t33 * t32, -t33 * t36, t37, -t33 * pkin(2) + 0; t36, -t32, 0, t28; t43, t18, t33, t40; -t33 * t25, -t46, t37, t39; t26, -t25, 0, t38; t13, -t12, t18, t6; t11, t10, -t46, t5; t47, -t48, -t25, t7; t4, t3, t12, t6; t2, t1, -t10, t5; t9, t8, t48, t7; -t12 * t19 - t4 * t20, t12 * t20 - t4 * t19, t3, t3 * pkin(6) + t6; t10 * t19 - t2 * t20, -t10 * t20 - t2 * t19, t1, t1 * pkin(6) + t5; -t19 * t48 - t9 * t20, -t9 * t19 + t20 * t48, t8, t8 * pkin(6) + t7;];
Tc_stack = t14;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,6+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
