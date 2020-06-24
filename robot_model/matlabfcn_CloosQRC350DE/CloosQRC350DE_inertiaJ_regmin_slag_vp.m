% Calculate minimal parameter regressor of joint inertia matrix for
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
% MM_reg [((6+1)*6/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = CloosQRC350DE_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:37
% EndTime: 2020-06-23 21:07:38
% DurationCPUTime: 0.70s
% Computational Cost: add. (202->43), mult. (476->95), div. (0->0), fcn. (550->10), ass. (0->51)
t33 = sin(qJ(3));
t34 = sin(qJ(2));
t37 = cos(qJ(3));
t38 = cos(qJ(2));
t12 = t34 * t33 - t38 * t37;
t51 = 0.2e1 * t12;
t50 = pkin(3) * t33;
t24 = pkin(7) * qJ(5) - qJ(6);
t17 = sin(t24);
t18 = cos(t24);
t36 = cos(qJ(4));
t32 = sin(qJ(4));
t35 = cos(qJ(5));
t44 = t32 * t35;
t8 = -t17 * t36 + t18 * t44;
t49 = t8 * t17;
t31 = sin(qJ(5));
t48 = t18 * t31;
t28 = t32 ^ 2;
t47 = t28 * t35;
t20 = t31 * t32;
t46 = t31 * t36;
t45 = t32 * t12;
t43 = t35 * t36;
t13 = t38 * t33 + t37 * t34;
t4 = t12 * t43 + t31 * t13;
t42 = t4 * t44;
t41 = t8 * t48;
t23 = -pkin(5) + t50;
t40 = t23 * t47;
t39 = pkin(7) * t20;
t22 = -t34 * pkin(3) - pkin(2);
t30 = t36 ^ 2;
t29 = t35 ^ 2;
t27 = t31 ^ 2;
t26 = pkin(3) * t37;
t21 = t29 * t28;
t19 = t27 * t28;
t16 = pkin(5) * t47;
t15 = t31 * t44;
t14 = -t31 * pkin(4) + pkin(5) * t43;
t11 = t12 ^ 2;
t10 = t36 * t45;
t9 = t23 * t43 + t31 * (t26 + pkin(4));
t7 = t8 ^ 2;
t6 = -t13 * pkin(4) + t12 * pkin(5) + t22;
t5 = -t12 * t46 + t35 * t13;
t3 = t5 * t20;
t2 = -t17 * t45 - t18 * t4;
t1 = t2 * t8;
t25 = [1, t38 ^ 2, 0, 0, 0.2e1 * pkin(2) * t34, t11, t13 * t51, 0, 0, 0, -0.2e1 * t22 * t13, t22 * t51, t30 * t11, t13 ^ 2, t4 ^ 2, t28 * t11, 0.2e1 * (t12 * t47 + t36 * t4) * t6, t2 ^ 2, t5 ^ 2; 0, 0, -t38, 0, 0, 0, 0, t12, t13, 0, 0, 0, -t10, 0, -t42, t10, (-t12 * t9 + t23 * t4) * t32, t1, t3; 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t50, t28, 0, t21, t30, -0.2e1 * t9 * t36 - 0.2e1 * t40, t7, t19; 0, 0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, -t10, 0, -t42, t10, (-pkin(5) * t4 + t12 * t14) * t32, t1, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, -t50, t28, 0, t21, t30, -t40 + t16 + (t14 - t9) * t36, t7, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t28, 0, t21, t30, 0.2e1 * t14 * t36 + 0.2e1 * t16, t7, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t4 * t31, 0, t6 * t46, -t2 * t48, t5 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, t23 * t20, -t41, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, -pkin(5) * t20, -t41, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, 0, 0, t18 ^ 2 * t27, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t6 * t44, -t2 * t17, -t5 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t9, -t49, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t14, -t49, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t48, -t35 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t17 ^ 2, pkin(7) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t25;
