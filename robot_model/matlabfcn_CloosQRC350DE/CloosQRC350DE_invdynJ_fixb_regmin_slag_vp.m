% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% tau_reg [6x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = CloosQRC350DE_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:31
% EndTime: 2020-06-23 21:07:35
% DurationCPUTime: 2.99s
% Computational Cost: add. (2247->200), mult. (4853->315), div. (0->0), fcn. (3892->10), ass. (0->119)
t83 = sin(qJ(3));
t84 = sin(qJ(2));
t132 = t84 * t83;
t107 = qJD(1) * t132;
t88 = cos(qJ(2));
t127 = qJD(1) * t88;
t87 = cos(qJ(3));
t50 = t87 * t127 - t107;
t82 = sin(qJ(4));
t149 = -qJD(5) * t82 + t50;
t85 = cos(qJ(5));
t86 = cos(qJ(4));
t131 = t85 * t86;
t130 = t87 * t84;
t55 = t88 * t83 + t130;
t51 = qJD(1) * t55;
t81 = sin(qJ(5));
t99 = t51 * t131 + t149 * t81;
t114 = qJD(1) * qJD(2);
t105 = t88 * t114;
t148 = qJDD(1) * t84 + t105;
t113 = qJD(1) * qJD(3);
t147 = t88 * t113 + t148;
t80 = qJD(2) + qJD(3);
t43 = t82 * t50 - t86 * t80;
t42 = qJD(5) - t43;
t146 = t42 * t81;
t44 = -t86 * t50 - t82 * t80;
t49 = qJD(4) + t51;
t30 = -t81 * t44 + t85 * t49;
t77 = pkin(7) * qJD(5) - qJD(6);
t26 = t30 - t77;
t145 = pkin(3) * t83;
t144 = pkin(3) * t87;
t106 = t84 * t114;
t115 = t88 * qJDD(1);
t28 = t113 * t130 + (t106 - t115) * t87 + t147 * t83;
t29 = -t80 * t107 + t83 * t115 + t147 * t87;
t73 = -t84 * pkin(3) - pkin(2);
t13 = -pkin(3) * t105 - t29 * pkin(4) + t28 * pkin(5) + t73 * qJDD(1);
t118 = t73 * qJD(1);
t33 = -t51 * pkin(4) - t50 * pkin(5) + t118;
t128 = pkin(3) * qJD(2);
t111 = t83 * t128;
t56 = -t80 * pkin(5) + t111;
t25 = -t82 * t33 + t86 * t56;
t116 = qJDD(2) * t83;
t125 = qJD(3) * t87;
t79 = qJDD(2) + qJDD(3);
t47 = -t79 * pkin(5) + (qJD(2) * t125 + t116) * pkin(3);
t7 = t25 * qJD(4) + t86 * t13 + t82 * t47;
t143 = t7 * t82;
t142 = t42 * t85;
t141 = t49 * t50;
t140 = t50 * t51;
t54 = -t88 * t87 + t132;
t139 = t54 * t82;
t138 = t54 * t86;
t31 = t85 * t44 + t81 * t49;
t78 = pkin(7) * qJ(5) - qJ(6);
t70 = sin(t78);
t137 = t70 * t31;
t71 = cos(t78);
t136 = t71 * t42;
t134 = t81 * t86;
t133 = t82 * t85;
t129 = t88 * qJD(1) ^ 2;
t126 = qJD(2) * t88;
t124 = qJD(4) * t82;
t123 = qJD(4) * t85;
t122 = qJD(4) * t86;
t120 = qJD(5) * t86;
t24 = t86 * t33 + t82 * t56;
t119 = t24 * qJD(4);
t112 = pkin(3) * t125;
t110 = t87 * t128;
t109 = pkin(3) * t126;
t38 = -t50 * pkin(4) + t51 * pkin(5);
t74 = -pkin(5) + t145;
t103 = pkin(3) * t127 + qJD(4) * t74 - t38;
t102 = qJD(2) * (-qJD(3) + t80);
t101 = qJD(3) * (-qJD(2) - t80);
t100 = -t54 * g(3) - t51 * t118;
t41 = -t88 * t125 - t87 * t126 + t80 * t132;
t98 = -t54 * t120 - t41;
t76 = qJDD(2) * t144;
t97 = t81 * (t79 * pkin(4) - qJD(3) * t111 + t76) + t85 * (-t82 * t13 + t86 * t47 - t119);
t57 = t80 * pkin(4) + t110;
t96 = t81 * t25 - t85 * t57;
t95 = t55 * g(3) + t50 * t118 + t76;
t40 = t80 * t55;
t94 = t54 * t122 + t40 * t82;
t93 = -t54 * t124 + t40 * t86;
t19 = t85 * t25 + t81 * t57;
t5 = -t96 * qJD(5) + t97;
t92 = -t5 * t86 + (-t55 * t134 - t85 * t54) * g(3) - t99 * t24 + (t82 * t51 + t124) * t19;
t16 = t44 * qJD(4) + t82 * t28 + t86 * t79 + qJDD(5);
t91 = -qJD(4) * t31 - qJD(5) * t146 + t85 * t16;
t90 = -qJD(5) * t55 - t93;
t75 = pkin(4) + t144;
t39 = -t55 * pkin(4) + t54 * pkin(5) + t73;
t37 = t54 * t131 + t81 * t55;
t32 = t50 ^ 2 - t51 ^ 2;
t27 = qJDD(4) + t29;
t22 = t41 * pkin(4) + t40 * pkin(5) - t109;
t21 = -t50 * t80 + t29;
t20 = -t51 * t80 + t28;
t17 = t43 * qJD(4) + t86 * t28 - t82 * t79;
t15 = -t71 * t31 - t70 * t42;
t12 = t98 * t81 - t90 * t85;
t11 = -t49 * t86 * t44 - t17 * t82;
t10 = -t49 * t82 * t42 + t16 * t86;
t9 = t30 * qJD(5) + t85 * t17 + t81 * t27;
t8 = -pkin(7) * qJDD(5) - t31 * qJD(5) - t81 * t17 + t85 * t27 + qJDD(6);
t4 = (-t42 * t77 - t9) * t71 + (t31 * t77 - t16) * t70;
t3 = -t9 * t133 + (-t85 * t122 - t99) * t31;
t2 = t8 * t81 * t82 + (t49 * t134 - t149 * t85) * t26;
t1 = t4 * (t71 * t133 - t70 * t86) + ((-t77 * t85 + t49) * t82 * t70 + ((-t77 + t123) * t86 + t99) * t71) * t15;
t6 = [qJDD(1), (-0.2e1 * t106 + t115) * t88, qJD(2) ^ 2 * t84 - t88 * qJDD(2), 0, 0.2e1 * t148 * pkin(2), t28 * t54 - t50 * t40, t28 * t55 + t54 * t29 + t40 * t51 + t50 * t41, t40 * t80 + t54 * t79, -t41 * t80 + t55 * t79, 0, 0.2e1 * t51 * t109 + (qJD(1) * t41 - qJDD(1) * t55 - t29) * t73, (-qJD(1) * t54 + t50) * t109 + (qJD(1) * t40 + qJDD(1) * t54 + t28) * t73, t17 * t138 + t93 * t44, t27 * t55 - t49 * t41, t31 * t12 + t9 * t37, t16 * t139 + t94 * t42, t24 * t12 + t7 * t37 + (t22 * t31 + t39 * t9 + (t39 * t142 - t19 * t54) * qJD(4)) * t86 + (t22 * t142 - t19 * t40 + t91 * t39 - t5 * t54) * t82, (-t4 * t37 + t15 * (-t77 * t139 - t12)) * t71 + (-t4 * t139 + t15 * (t37 * t77 - t94)) * t70, (t26 * t98 + t8 * t55) * t85 + (-t8 * t138 + t26 * t90) * t81; 0, t84 * t129, -t115, qJDD(2), -pkin(2) * t129 + t84 * g(3), t140, t32, t20, t21, t79, (t83 * t101 - t51 * t127 + t79 * t87) * pkin(3) + t95, (-t50 * t127 + (-qJDD(2) - t79) * t83 + t87 * t101) * pkin(3) + t100, t11, -t141, t3, t10, t82 * t74 * t9 + (-(-qJD(3) * t145 - t74 * t120) * t42 - t75 * t16) * t81 + (t103 * t86 + t82 * t112) * t31 + (-t143 + (-t74 * t16 - t119) * t86 + (-qJD(5) * t75 + t103 * t82 - t86 * t112) * t42) * t85 + t92, t1, t2; 0, 0, 0, 0, 0, t140, t32, t20, t21, t79, t102 * t145 + t95, (t87 * t102 - t116) * pkin(3) + t100, t11, -t141, t3, t10, -(t82 * t110 + t86 * t38) * t31 + (-pkin(4) * t16 - t42 * t111) * t81 + (-t86 * t119 - t143 + (-pkin(4) * qJD(5) + t86 * t110 - t82 * t38) * t42) * t85 + ((-t42 * t123 - t9) * t82 + t91 * t86) * pkin(5) + t92, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t43, t27, t31 * t142 + t9 * t81, -t42 * t44, t19 * t44 - t25 * t31 + (g(3) * t139 + t7) * t81, -t4 * t71 * t81 + (-t85 * t136 + (t77 * t81 + t44) * t70) * t15, -t26 * t146 + t8 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t30, t16, -t37 * g(3) - t24 * t30 - t97 + (qJD(5) - t42) * t96, -t4 * t70 + (-pkin(7) * t137 + (pkin(7) * t42 + t26) * t71) * t15, -t8 * pkin(7) + t26 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * (t136 - t137), t8;];
tau_reg = t6;
