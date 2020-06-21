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
% taug_reg [6x36]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
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
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:35
% EndTime: 2020-06-19 21:39:35
% DurationCPUTime: 0.28s
% Computational Cost: add. (125->36), mult. (150->60), div. (0->0), fcn. (152->10), ass. (0->244)
unknown=NaN(6,36);
t1 = sin(qJ(2));
t3 = cos(qJ(2));
t5 = qJ(2) + qJ(3);
t6 = sin(t5);
t7 = g(3) * t6;
t8 = cos(t5);
t9 = g(3) * t8;
t10 = cos(qJ(4));
t11 = t7 * t10;
t12 = sin(qJ(4));
t13 = t7 * t12;
t14 = t6 * t10;
t15 = cos(qJ(5));
t17 = sin(qJ(5));
t19 = -t14 * t15 - t8 * t17;
t20 = g(3) * t19;
t24 = g(3) * (t14 * t17 - t8 * t15);
t26 = pkin(7) * qJ(5) - qJ(6);
t27 = cos(t26);
t29 = t6 * t12;
t30 = sin(t26);
t33 = g(3) * (-t19 * t27 + t29 * t30);
t37 = g(3) * (-t19 * t30 - t29 * t27);
t44 = t8 * t12;
t47 = t8 * t10;
t58 = -t6 * t15 - t47 * t17;
t62 = -t47 * t15 + t6 * t17;
t65 = -t62 * pkin(7);
unknown(1,1) = 0;
unknown(1,2) = 0;
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = 0;
unknown(1,6) = 0;
unknown(1,7) = 0;
unknown(1,8) = 0;
unknown(1,9) = 0;
unknown(1,10) = 0;
unknown(1,11) = 0;
unknown(1,12) = 0;
unknown(1,13) = 0;
unknown(1,14) = 0;
unknown(1,15) = 0;
unknown(1,16) = 0;
unknown(1,17) = 0;
unknown(1,18) = 0;
unknown(1,19) = 0;
unknown(1,20) = 0;
unknown(1,21) = 0;
unknown(1,22) = 0;
unknown(1,23) = 0;
unknown(1,24) = 0;
unknown(1,25) = 0;
unknown(1,26) = 0;
unknown(1,27) = 0;
unknown(1,28) = 0;
unknown(1,29) = 0;
unknown(1,30) = 0;
unknown(1,31) = 0;
unknown(1,32) = 0;
unknown(1,33) = 0;
unknown(1,34) = 0;
unknown(1,35) = 0;
unknown(1,36) = 0;
unknown(2,1) = 0;
unknown(2,2) = 0;
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(2,5) = 0;
unknown(2,6) = 0;
unknown(2,7) = (g(3) * t1);
unknown(2,8) = (g(3) * t3);
unknown(2,9) = 0;
unknown(2,10) = 0;
unknown(2,11) = 0;
unknown(2,12) = 0;
unknown(2,13) = 0;
unknown(2,14) = t7;
unknown(2,15) = t9;
unknown(2,16) = 0;
unknown(2,17) = 0;
unknown(2,18) = 0;
unknown(2,19) = 0;
unknown(2,20) = 0;
unknown(2,21) = t11;
unknown(2,22) = -t13;
unknown(2,23) = 0;
unknown(2,24) = 0;
unknown(2,25) = 0;
unknown(2,26) = 0;
unknown(2,27) = 0;
unknown(2,28) = -t20;
unknown(2,29) = -t24;
unknown(2,30) = 0;
unknown(2,31) = 0;
unknown(2,32) = 0;
unknown(2,33) = 0;
unknown(2,34) = 0;
unknown(2,35) = -t33;
unknown(2,36) = -t37;
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 0;
unknown(3,4) = 0;
unknown(3,5) = 0;
unknown(3,6) = 0;
unknown(3,7) = 0;
unknown(3,8) = 0;
unknown(3,9) = 0;
unknown(3,10) = 0;
unknown(3,11) = 0;
unknown(3,12) = 0;
unknown(3,13) = 0;
unknown(3,14) = t7;
unknown(3,15) = t9;
unknown(3,16) = 0;
unknown(3,17) = 0;
unknown(3,18) = 0;
unknown(3,19) = 0;
unknown(3,20) = 0;
unknown(3,21) = t11;
unknown(3,22) = -t13;
unknown(3,23) = 0;
unknown(3,24) = 0;
unknown(3,25) = 0;
unknown(3,26) = 0;
unknown(3,27) = 0;
unknown(3,28) = -t20;
unknown(3,29) = -t24;
unknown(3,30) = 0;
unknown(3,31) = 0;
unknown(3,32) = 0;
unknown(3,33) = 0;
unknown(3,34) = 0;
unknown(3,35) = -t33;
unknown(3,36) = -t37;
unknown(4,1) = 0;
unknown(4,2) = 0;
unknown(4,3) = 0;
unknown(4,4) = 0;
unknown(4,5) = 0;
unknown(4,6) = 0;
unknown(4,7) = 0;
unknown(4,8) = 0;
unknown(4,9) = 0;
unknown(4,10) = 0;
unknown(4,11) = 0;
unknown(4,12) = 0;
unknown(4,13) = 0;
unknown(4,14) = 0;
unknown(4,15) = 0;
unknown(4,16) = 0;
unknown(4,17) = 0;
unknown(4,18) = 0;
unknown(4,19) = 0;
unknown(4,20) = 0;
unknown(4,21) = (t9 * t12);
unknown(4,22) = (t9 * t10);
unknown(4,23) = 0;
unknown(4,24) = 0;
unknown(4,25) = 0;
unknown(4,26) = 0;
unknown(4,27) = 0;
unknown(4,28) = (t9 * t12 * t15);
unknown(4,29) = -(t9 * t12 * t17);
unknown(4,30) = 0;
unknown(4,31) = 0;
unknown(4,32) = 0;
unknown(4,33) = 0;
unknown(4,34) = 0;
unknown(4,35) = -(g(3) * (t44 * t15 * t27 - t47 * t30));
unknown(4,36) = -(g(3) * (t44 * t15 * t30 + t47 * t27));
unknown(5,1) = 0;
unknown(5,2) = 0;
unknown(5,3) = 0;
unknown(5,4) = 0;
unknown(5,5) = 0;
unknown(5,6) = 0;
unknown(5,7) = 0;
unknown(5,8) = 0;
unknown(5,9) = 0;
unknown(5,10) = 0;
unknown(5,11) = 0;
unknown(5,12) = 0;
unknown(5,13) = 0;
unknown(5,14) = 0;
unknown(5,15) = 0;
unknown(5,16) = 0;
unknown(5,17) = 0;
unknown(5,18) = 0;
unknown(5,19) = 0;
unknown(5,20) = 0;
unknown(5,21) = 0;
unknown(5,22) = 0;
unknown(5,23) = 0;
unknown(5,24) = 0;
unknown(5,25) = 0;
unknown(5,26) = 0;
unknown(5,27) = 0;
unknown(5,28) = -(g(3) * t58);
unknown(5,29) = -(g(3) * t62);
unknown(5,30) = 0;
unknown(5,31) = 0;
unknown(5,32) = 0;
unknown(5,33) = 0;
unknown(5,34) = 0;
unknown(5,35) = -(g(3) * (-t44 * pkin(7) * t27 - t58 * t27 + t65 * t30));
unknown(5,36) = -(g(3) * (-t44 * pkin(7) * t30 - t65 * t27 - t58 * t30));
unknown(6,1) = 0;
unknown(6,2) = 0;
unknown(6,3) = 0;
unknown(6,4) = 0;
unknown(6,5) = 0;
unknown(6,6) = 0;
unknown(6,7) = 0;
unknown(6,8) = 0;
unknown(6,9) = 0;
unknown(6,10) = 0;
unknown(6,11) = 0;
unknown(6,12) = 0;
unknown(6,13) = 0;
unknown(6,14) = 0;
unknown(6,15) = 0;
unknown(6,16) = 0;
unknown(6,17) = 0;
unknown(6,18) = 0;
unknown(6,19) = 0;
unknown(6,20) = 0;
unknown(6,21) = 0;
unknown(6,22) = 0;
unknown(6,23) = 0;
unknown(6,24) = 0;
unknown(6,25) = 0;
unknown(6,26) = 0;
unknown(6,27) = 0;
unknown(6,28) = 0;
unknown(6,29) = 0;
unknown(6,30) = 0;
unknown(6,31) = 0;
unknown(6,32) = 0;
unknown(6,33) = 0;
unknown(6,34) = 0;
unknown(6,35) = -(g(3) * (t44 * t27 + t62 * t30));
unknown(6,36) = -(g(3) * (-t62 * t27 + t44 * t30));
taug_reg = unknown;
