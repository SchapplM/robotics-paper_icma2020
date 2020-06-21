% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [6x36]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = CloosQRC350OL_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:18:57
% EndTime: 2020-06-20 08:18:58
% DurationCPUTime: 0.28s
% Computational Cost: add. (90->28), mult. (112->49), div. (0->0), fcn. (134->10), ass. (0->29)
t12 = qJ(2) + qJ(3);
t10 = sin(t12);
t13 = sin(qJ(6));
t15 = sin(qJ(4));
t16 = cos(qJ(6));
t23 = t15 * t16;
t11 = cos(t12);
t14 = sin(qJ(5));
t17 = cos(qJ(5));
t18 = cos(qJ(4));
t21 = t17 * t18;
t4 = -t10 * t21 - t11 * t14;
t32 = g(3) * (-t10 * t23 + t4 * t13);
t25 = t13 * t15;
t31 = g(3) * (-t10 * t25 - t4 * t16);
t24 = t14 * t18;
t30 = g(3) * (t10 * t24 - t11 * t17);
t29 = g(3) * t4;
t28 = g(3) * (-t10 * t17 - t11 * t24);
t9 = g(3) * t11;
t27 = g(3) * t15;
t26 = g(3) * t18;
t22 = t15 * t17;
t20 = t10 * t27;
t19 = t11 * t27;
t8 = g(3) * t10;
t7 = t10 * t26;
t6 = -t10 * t14 + t11 * t21;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, g(3) * sin(qJ(2)), g(3) * cos(qJ(2)), 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t20, 0, 0, 0, 0, 0, -t29, -t30, 0, 0, 0, 0, 0, -t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t20, 0, 0, 0, 0, 0, -t29, -t30, 0, 0, 0, 0, 0, -t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t11 * t26, 0, 0, 0, 0, 0, t17 * t19, -t14 * t19, 0, 0, 0, 0, 0, -(t13 * t18 + t16 * t22) * t9, -(-t13 * t22 + t16 * t18) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, g(3) * t6, 0, 0, 0, 0, 0, t16 * t28, -t13 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * (t11 * t23 + t6 * t13), -g(3) * (-t11 * t25 + t6 * t16);];
taug_reg = t1;
