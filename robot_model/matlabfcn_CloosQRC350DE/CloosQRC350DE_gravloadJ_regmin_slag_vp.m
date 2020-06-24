% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [6x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = CloosQRC350DE_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:36
% EndTime: 2020-06-23 21:07:36
% DurationCPUTime: 0.16s
% Computational Cost: add. (14->6), mult. (42->15), div. (0->0), fcn. (56->8), ass. (0->14)
t12 = cos(qJ(3));
t13 = cos(qJ(2));
t8 = sin(qJ(3));
t9 = sin(qJ(2));
t4 = -t13 * t12 + t9 * t8;
t2 = t4 * g(3);
t5 = t12 * t9 + t13 * t8;
t7 = sin(qJ(5));
t15 = t7 * t5;
t14 = cos(qJ(5)) * t4;
t11 = cos(qJ(4));
t3 = t5 * g(3);
t1 = (-t11 * t15 - t14) * g(3);
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t9 * g(3), 0, 0, 0, 0, 0, t3, -t2, 0, 0, 0, 0, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2, 0, 0, 0, 0, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, sin(qJ(4)) * t7 * t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t11 * t14 - t15) * g(3), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t6;
