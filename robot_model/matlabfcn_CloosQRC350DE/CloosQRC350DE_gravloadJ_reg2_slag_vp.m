% Calculate inertial parameters regressor of gravitation load for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = CloosQRC350DE_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:32
% EndTime: 2020-06-23 21:07:32
% DurationCPUTime: 0.28s
% Computational Cost: add. (66->29), mult. (148->44), div. (0->0), fcn. (169->8), ass. (0->24)
t14 = sin(qJ(5));
t19 = cos(qJ(4));
t18 = cos(qJ(5));
t16 = sin(qJ(3));
t17 = sin(qJ(2));
t20 = cos(qJ(3));
t21 = cos(qJ(2));
t8 = t17 * t16 - t21 * t20;
t25 = t18 * t8;
t9 = t21 * t16 + t20 * t17;
t28 = g(3) * (t14 * t9 + t19 * t25);
t24 = t14 * t19;
t27 = g(3) * (t9 * t24 + t25);
t6 = t8 * g(3);
t26 = t17 * g(3);
t5 = sin(qJ(4)) * t14 * t6;
t12 = pkin(6) * t24 - pkin(4);
t13 = pkin(6) * t18 + pkin(5);
t23 = -t12 * t20 - t13 * t16;
t22 = (-t12 * t16 + t13 * t20) * t21;
t7 = t9 * g(3);
t4 = g(3) * ((pkin(4) * t17 + pkin(5) * t21) * t20 + t16 * (pkin(4) * t21 - pkin(5) * t17));
t3 = g(3) * ((pkin(4) * t20 - pkin(5) * t16 + pkin(3)) * t17 + (pkin(4) * t16 + pkin(5) * t20) * t21);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, pkin(3) * t26, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t3, 0, 0, 0, 0, 0, 0, 0, -t27, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t27, ((pkin(3) + t23) * t17 + t22) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t4, 0, 0, 0, 0, 0, 0, 0, -t27, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, -t27, g(3) * (t23 * t17 + t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, pkin(6) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -pkin(6) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t1;
