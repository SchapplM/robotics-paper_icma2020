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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 22:03:19
% EndTime: 2020-06-23 22:03:19
% DurationCPUTime: 0.27s
% Computational Cost: add. (75->23), mult. (77->22), div. (0->0), fcn. (77->7), ass. (0->21)
t15 = cos(qJ(4));
t10 = qJ(2) + qJ(3);
t9 = cos(t10);
t17 = t9 * cos(qJ(5));
t11 = sin(qJ(5));
t8 = sin(t10);
t18 = t8 * t11;
t2 = t15 * t18 - t17;
t23 = g(3) * t2;
t22 = g(3) * (-t15 * t17 + t18);
t13 = sin(qJ(2));
t19 = t13 * pkin(3);
t5 = -t8 * pkin(4) - t9 * pkin(5);
t4 = t5 - t19;
t21 = g(3) * t4;
t20 = g(3) * t5;
t7 = g(3) * t9;
t16 = sin(qJ(4)) * t11 * t7;
t6 = g(3) * t8;
t1 = t2 * pkin(6);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, g(3) * t19, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t21, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -g(3) * (t1 + t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t20, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -g(3) * (t1 + t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -pkin(6) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -pkin(6) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t3;
