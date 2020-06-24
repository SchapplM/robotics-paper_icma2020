% Calculate minimal parameter regressor of joint inertia matrix for
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
% MM_reg [((6+1)*6/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = CloosQRC350OL_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:18
% EndTime: 2020-06-23 22:03:20
% DurationCPUTime: 0.60s
% Computational Cost: add. (153->37), mult. (424->90), div. (0->0), fcn. (534->10), ass. (0->49)
t35 = sin(qJ(2));
t23 = t35 * pkin(3) + pkin(2);
t51 = 0.2e1 * t23;
t34 = sin(qJ(3));
t50 = t34 * pkin(3);
t33 = sin(qJ(4));
t28 = t33 ^ 2;
t37 = cos(qJ(5));
t49 = t28 * t37;
t32 = sin(qJ(5));
t20 = t32 * t33;
t36 = cos(qJ(6));
t48 = t32 * t36;
t38 = cos(qJ(4));
t47 = t32 * t38;
t39 = cos(qJ(3));
t40 = cos(qJ(2));
t15 = -t34 * t35 + t39 * t40;
t46 = t33 * t15;
t45 = t33 * t37;
t44 = t37 * t38;
t13 = t34 * t40 + t39 * t35;
t5 = -t32 * t13 + t15 * t44;
t43 = t5 * t45;
t31 = sin(qJ(6));
t11 = t31 * t38 + t36 * t45;
t42 = t11 * t48;
t22 = -pkin(5) + t50;
t41 = t22 * t49;
t30 = t38 ^ 2;
t29 = t37 ^ 2;
t27 = t32 ^ 2;
t26 = t39 * pkin(3);
t21 = t29 * t28;
t19 = t27 * t28;
t18 = pkin(5) * t49;
t17 = t32 * t45;
t16 = -t32 * pkin(4) + pkin(5) * t44;
t12 = t15 ^ 2;
t10 = t11 ^ 2;
t9 = t11 * t31;
t8 = t38 * t46;
t7 = t22 * t44 + t32 * (t26 + pkin(4));
t6 = t13 * pkin(4) + t15 * pkin(5) + t23;
t4 = -t37 * t13 - t15 * t47;
t3 = t4 * t20;
t2 = t31 * t46 - t36 * t5;
t1 = t2 * t11;
t14 = [1, t40 ^ 2, 0, 0, 0.2e1 * pkin(2) * t35, t12, -0.2e1 * t15 * t13, 0, 0, 0, t13 * t51, t15 * t51, t30 * t12, t13 ^ 2, t5 ^ 2, t28 * t12, 0.2e1 * (t15 * t49 + t38 * t5) * t6, t2 ^ 2, t4 ^ 2; 0, 0, t40, 0, 0, 0, 0, t15, -t13, 0, 0, 0, -t8, 0, -t43, t8, (-t15 * t7 + t22 * t5) * t33, t1, t3; 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t50, t28, 0, t21, t30, -0.2e1 * t7 * t38 - 0.2e1 * t41, t10, t19; 0, 0, 0, 0, 0, 0, 0, t15, -t13, 0, 0, 0, -t8, 0, -t43, t8, (-pkin(5) * t5 + t15 * t16) * t33, t1, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, -t50, t28, 0, t21, t30, -t41 + t18 + (t16 - t7) * t38, t10, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t28, 0, t21, t30, 0.2e1 * t16 * t38 + 0.2e1 * t18, t10, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t5 * t32, 0, t6 * t47, -t2 * t48, t4 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, t22 * t20, -t42, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, -pkin(5) * t20, -t42, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, 0, 0, t36 ^ 2 * t27, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t6 * t45, t2 * t31, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t7, t9, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t16, t9, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t48, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t31 ^ 2, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t14;
