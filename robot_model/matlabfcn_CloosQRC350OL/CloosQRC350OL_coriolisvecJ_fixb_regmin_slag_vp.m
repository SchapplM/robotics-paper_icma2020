% Calculate minimal parameter regressor of coriolis joint torque vector for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% tauc_reg [6x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = CloosQRC350OL_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:18
% EndTime: 2020-06-23 22:03:22
% DurationCPUTime: 1.95s
% Computational Cost: add. (1745->166), mult. (4232->283), div. (0->0), fcn. (3386->10), ass. (0->112)
t68 = sin(qJ(4));
t106 = qJD(5) * t68;
t71 = cos(qJ(5));
t72 = cos(qJ(4));
t118 = t71 * t72;
t69 = sin(qJ(3));
t129 = sin(qJ(2));
t94 = qJD(1) * t129;
t89 = t69 * t94;
t74 = cos(qJ(2));
t112 = qJD(1) * t74;
t73 = cos(qJ(3));
t98 = t73 * t112;
t46 = t89 - t98;
t80 = -t73 * t129 - t69 * t74;
t47 = t80 * qJD(1);
t67 = sin(qJ(5));
t31 = t47 * t118 + t67 * t46;
t131 = t67 * t106 - t31;
t64 = qJD(2) + qJD(3);
t42 = -t72 * t46 - t68 * t64;
t45 = qJD(4) + t47;
t27 = -t67 * t42 + t71 * t45;
t38 = t64 * t80;
t33 = t38 * qJD(1);
t93 = qJD(2) * t129;
t86 = qJD(1) * t93;
t115 = -qJD(3) * t89 - t69 * t86;
t117 = t73 * t74;
t85 = t64 * t117;
t34 = qJD(1) * t85 + t115;
t111 = qJD(2) * t74;
t95 = qJD(1) * t111;
t16 = pkin(3) * t95 + t34 * pkin(4) + t33 * pkin(5);
t53 = qJD(1) * pkin(2) + pkin(3) * t94;
t30 = -t47 * pkin(4) - t46 * pkin(5) + t53;
t114 = pkin(3) * qJD(2);
t102 = t69 * t114;
t51 = -t64 * pkin(5) + t102;
t23 = -t68 * t30 + t72 * t51;
t113 = pkin(3) * qJD(3);
t96 = qJD(2) * t113;
t90 = t73 * t96;
t9 = t23 * qJD(4) + t72 * t16 + t68 * t90;
t130 = t9 * t68;
t61 = t129 * pkin(3) + pkin(2);
t128 = t45 * t46;
t127 = t46 * t47;
t97 = t69 * t129;
t50 = -t97 + t117;
t126 = t50 * t68;
t125 = t50 * t72;
t124 = t53 * t46;
t123 = t53 * t47;
t121 = t67 * t45;
t120 = t68 * t71;
t41 = t68 * t46 - t72 * t64;
t40 = qJD(5) - t41;
t119 = t71 * t40;
t116 = t74 * qJD(1) ^ 2;
t110 = qJD(4) * t68;
t109 = qJD(4) * t71;
t108 = qJD(4) * t72;
t107 = qJD(5) * t67;
t105 = qJD(5) * t72;
t22 = t72 * t30 + t68 * t51;
t104 = t22 * qJD(4);
t26 = qJD(6) + t27;
t103 = t73 * t113;
t101 = t73 * t114;
t100 = pkin(3) * t111;
t36 = -t46 * pkin(4) + t47 * pkin(5);
t59 = t69 * pkin(3) - pkin(5);
t92 = -pkin(3) * t112 + qJD(4) * t59 - t36;
t91 = qJD(3) * (-qJD(2) - t64);
t39 = -qJD(3) * t97 - t69 * t93 + t85;
t88 = -t50 * t105 - t39;
t87 = (-qJD(3) + t64) * t114;
t52 = t64 * pkin(4) + t101;
t84 = t67 * t23 - t71 * t52;
t28 = t71 * t42 + t121;
t66 = sin(qJ(6));
t70 = cos(qJ(6));
t83 = t66 * t28 + t70 * t40;
t82 = t50 * t108 + t38 * t68;
t81 = -t50 * t110 + t38 * t72;
t17 = t71 * t23 + t67 * t52;
t78 = -t67 * t69 * t96 + t71 * (-t68 * t16 + t72 * t90 - t104);
t5 = -t84 * qJD(5) + t78;
t79 = -t5 * t72 + t131 * t22 + (t47 * t68 + t110) * t17;
t21 = t42 * qJD(4) + t68 * t33;
t77 = -qJD(4) * t28 - t40 * t107 + t71 * t21;
t76 = -qJD(5) * t80 - t81;
t60 = t73 * pkin(3) + pkin(4);
t37 = -pkin(4) * t80 + t50 * pkin(5) + t61;
t35 = t50 * t118 + t67 * t80;
t29 = t46 ^ 2 - t47 ^ 2;
t25 = -t115 + (-t46 - t98) * t64;
t24 = -t47 * t64 + t33;
t20 = t41 * qJD(4) + t72 * t33;
t18 = t39 * pkin(4) + t38 * pkin(5) + t100;
t14 = -t70 * t28 + t66 * t40;
t12 = t88 * t67 - t76 * t71;
t11 = -t45 * t72 * t42 - t20 * t68;
t10 = -t45 * t68 * t40 + t21 * t72;
t7 = -t28 * qJD(5) - t67 * t20 - t71 * t34;
t6 = t27 * qJD(5) + t71 * t20 - t67 * t34;
t4 = t83 * qJD(6) + t66 * t21 - t70 * t6;
t3 = -t6 * t120 + (-t71 * t108 + t131) * t28;
t2 = t7 * t67 * t68 + ((-t46 + t106) * t71 + t72 * t121) * t26;
t1 = t4 * (t70 * t120 + t66 * t72) + ((t31 + (qJD(6) + t109) * t72) * t70 + (-t70 * t107 + (-qJD(6) * t71 - t45) * t66) * t68) * t14;
t8 = [0, -0.2e1 * t74 * t86, -qJD(2) ^ 2 * t129, 0, 0.2e1 * pkin(2) * t95, t33 * t50 - t46 * t38, t33 * t80 - t50 * t34 + t38 * t47 + t46 * t39, t38 * t64, -t39 * t64, 0, -0.2e1 * t47 * t100 + t61 * t34 + t53 * t39, t61 * t33 + t53 * t38 + (qJD(1) * t50 - t46) * t100, t20 * t125 + t81 * t42, -t34 * t80 - t45 * t39, t28 * t12 + t6 * t35, t21 * t126 + t82 * t40, t22 * t12 + t9 * t35 + (t18 * t28 + t37 * t6 + (t37 * t119 - t17 * t50) * qJD(4)) * t72 + (t18 * t119 - t17 * t38 + t77 * t37 - t5 * t50) * t68, (-t4 * t35 + t14 * (qJD(6) * t126 - t12)) * t70 + (t4 * t126 + t14 * (qJD(6) * t35 + t82)) * t66, (t26 * t88 + t7 * t80) * t71 + (-t7 * t125 + t26 * t76) * t67; 0, t129 * t116, 0, 0, -pkin(2) * t116, t127, t29, t24, t25, 0, t124 + (t47 * t112 + t69 * t91) * pkin(3), -t123 + (t46 * t112 + t73 * t91) * pkin(3), t11, -t128, t3, t10, t68 * t59 * t6 + (-(-t59 * t105 - t69 * t113) * t40 - t60 * t21) * t67 + (t68 * t103 + t92 * t72) * t28 + (-t130 + (-t59 * t21 - t104) * t72 + (-qJD(5) * t60 - t72 * t103 + t92 * t68) * t40) * t71 + t79, t1, t2; 0, 0, 0, 0, 0, t127, t29, t24, t25, 0, t69 * t87 + t124, t73 * t87 - t123, t11, -t128, t3, t10, -(t68 * t101 + t72 * t36) * t28 + (-pkin(4) * t21 - t40 * t102) * t67 + (-t72 * t104 - t130 + (-pkin(4) * qJD(5) + t72 * t101 - t68 * t36) * t40) * t71 + ((-t40 * t109 - t6) * t68 + t77 * t72) * pkin(5) + t79, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * t41, -t34, t28 * t119 + t6 * t67, -t40 * t42, t17 * t42 - t23 * t28 + t9 * t67, -t4 * t70 * t67 + (-t70 * t119 + (qJD(6) * t67 - t42) * t66) * t14, -t40 * t67 * t26 + t7 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t27, t21, -t22 * t27 - t78 + (qJD(5) - t40) * t84, t26 * t70 * t14 + t4 * t66, t26 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t83, t7;];
tauc_reg = t8;
