% Calculate inertial parameters regressor of gravitation load for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = CloosQRC350OL_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:18:58
% EndTime: 2020-06-20 08:18:59
% DurationCPUTime: 0.44s
% Computational Cost: add. (151->43), mult. (175->60), div. (0->0), fcn. (194->10), ass. (0->40)
t17 = qJ(2) + qJ(3);
t15 = sin(t17);
t18 = sin(qJ(6));
t20 = sin(qJ(4));
t22 = cos(qJ(6));
t30 = t20 * t22;
t16 = cos(t17);
t19 = sin(qJ(5));
t23 = cos(qJ(5));
t24 = cos(qJ(4));
t28 = t23 * t24;
t5 = -t15 * t28 - t16 * t19;
t43 = g(3) * (-t15 * t30 + t5 * t18);
t32 = t18 * t20;
t42 = g(3) * (-t15 * t32 - t5 * t22);
t31 = t19 * t24;
t4 = t15 * t31 - t16 * t23;
t41 = g(3) * t4;
t40 = g(3) * t5;
t39 = g(3) * (-t15 * t23 - t16 * t31);
t7 = -t15 * t19 + t16 * t28;
t38 = g(3) * t7;
t21 = sin(qJ(2));
t33 = t21 * pkin(3);
t9 = -t15 * pkin(4) - t16 * pkin(5);
t8 = t9 - t33;
t37 = g(3) * t8;
t36 = g(3) * t9;
t14 = g(3) * t16;
t35 = g(3) * t20;
t34 = g(3) * t24;
t29 = t20 * t23;
t27 = t16 * t35;
t26 = t16 * t34;
t25 = t19 * t27;
t13 = g(3) * t15;
t11 = t15 * t34;
t10 = t15 * t35;
t3 = t4 * pkin(6);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t21, g(3) * cos(qJ(2)), 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, g(3) * t33, 0, 0, 0, 0, 0, 0, t11, -t10, t14, -t37, 0, 0, 0, 0, 0, 0, -t40, -t41, t10, -t37, 0, 0, 0, 0, 0, 0, -t42, -t43, -t41, -g(3) * (t3 + t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, t14, -t36, 0, 0, 0, 0, 0, 0, -t40, -t41, t10, -t36, 0, 0, 0, 0, 0, 0, -t42, -t43, -t41, -g(3) * (t3 + t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t27, -t25, -t26, 0, 0, 0, 0, 0, 0, 0, -(t18 * t24 + t22 * t29) * t14, -(-t18 * t29 + t22 * t24) * t14, -t25, -pkin(6) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t39, -t18 * t39, t38, pkin(6) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * (t16 * t30 + t7 * t18), -g(3) * (-t16 * t32 + t7 * t22), 0, 0;];
taug_reg = t1;
