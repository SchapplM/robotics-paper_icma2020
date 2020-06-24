% Calculate minimal parameter regressor of coriolis joint torque vector for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% tauc_reg [6x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = CloosQRC350DE_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:31
% EndTime: 2020-06-23 21:07:35
% DurationCPUTime: 2.49s
% Computational Cost: add. (1881->167), mult. (4373->284), div. (0->0), fcn. (3443->10), ass. (0->108)
t74 = cos(qJ(2));
t108 = qJD(1) * t74;
t70 = sin(qJ(2));
t109 = qJD(1) * t70;
t69 = sin(qJ(3));
t73 = cos(qJ(3));
t46 = t73 * t108 - t69 * t109;
t68 = sin(qJ(4));
t127 = -qJD(5) * t68 + t46;
t71 = cos(qJ(5));
t72 = cos(qJ(4));
t113 = t71 * t72;
t49 = t74 * t69 + t73 * t70;
t47 = qJD(1) * t49;
t67 = sin(qJ(5));
t85 = t47 * t113 + t127 * t67;
t66 = qJD(2) + qJD(3);
t41 = t68 * t46 - t72 * t66;
t40 = qJD(5) - t41;
t126 = t40 * t67;
t42 = -t72 * t46 - t68 * t66;
t45 = qJD(4) + t47;
t27 = -t67 * t42 + t71 * t45;
t64 = pkin(7) * qJD(5) - qJD(6);
t26 = t27 - t64;
t115 = t70 * t69;
t48 = -t74 * t73 + t115;
t125 = qJD(1) * t48;
t106 = qJD(3) * t74;
t107 = qJD(2) * t74;
t89 = qJD(1) * t107;
t33 = t66 * t73 * t109 + (qJD(1) * t106 + t89) * t69;
t34 = t66 * t125;
t16 = -pkin(3) * t89 + t34 * pkin(4) + t33 * pkin(5);
t61 = -t70 * pkin(3) - pkin(2);
t99 = t61 * qJD(1);
t30 = -t47 * pkin(4) - t46 * pkin(5) + t99;
t111 = pkin(3) * qJD(2);
t97 = t69 * t111;
t50 = -t66 * pkin(5) + t97;
t23 = -t68 * t30 + t72 * t50;
t110 = pkin(3) * qJD(3);
t90 = qJD(2) * t110;
t86 = t73 * t90;
t9 = t23 * qJD(4) + t72 * t16 + t68 * t86;
t124 = t9 * t68;
t123 = t45 * t46;
t122 = t46 * t47;
t121 = t48 * t68;
t120 = t48 * t72;
t117 = t67 * t45;
t28 = t71 * t42 + t117;
t65 = pkin(7) * qJ(5) - qJ(6);
t59 = sin(t65);
t119 = t59 * t28;
t116 = t68 * t71;
t114 = t71 * t40;
t112 = t74 * qJD(1) ^ 2;
t105 = qJD(4) * t68;
t104 = qJD(4) * t71;
t103 = qJD(4) * t72;
t101 = qJD(5) * t72;
t22 = t72 * t30 + t68 * t50;
t100 = t22 * qJD(4);
t98 = t73 * t110;
t96 = t73 * t111;
t95 = pkin(3) * t107;
t92 = t46 * t99;
t91 = t47 * t99;
t36 = -t46 * pkin(4) + t47 * pkin(5);
t62 = pkin(3) * t69 - pkin(5);
t88 = pkin(3) * t108 + qJD(4) * t62 - t36;
t87 = qJD(3) * (-qJD(2) - t66);
t39 = t66 * t115 + (-t106 - t107) * t73;
t84 = -t48 * t101 - t39;
t83 = (-qJD(3) + t66) * t111;
t51 = t66 * pkin(4) + t96;
t82 = t67 * t23 - t71 * t51;
t38 = t66 * t49;
t81 = t48 * t103 + t38 * t68;
t80 = -t48 * t105 + t38 * t72;
t17 = t71 * t23 + t67 * t51;
t78 = -t67 * t69 * t90 + t71 * (-t68 * t16 + t72 * t86 - t100);
t5 = -t82 * qJD(5) + t78;
t79 = -t5 * t72 - t85 * t22 + (t47 * t68 + t105) * t17;
t21 = t42 * qJD(4) + t68 * t33;
t77 = -qJD(4) * t28 - qJD(5) * t126 + t71 * t21;
t76 = -qJD(5) * t49 - t80;
t63 = pkin(3) * t73 + pkin(4);
t60 = cos(t65);
t37 = -t49 * pkin(4) + t48 * pkin(5) + t61;
t35 = t48 * t113 + t67 * t49;
t29 = t46 ^ 2 - t47 ^ 2;
t25 = -t46 * t66 - t34;
t24 = -t47 * t66 + t33;
t20 = t41 * qJD(4) + t72 * t33;
t18 = t39 * pkin(4) + t38 * pkin(5) - t95;
t14 = -t60 * t28 - t59 * t40;
t12 = t84 * t67 - t76 * t71;
t11 = -t45 * t72 * t42 - t20 * t68;
t10 = -t45 * t68 * t40 + t21 * t72;
t7 = -t28 * qJD(5) - t67 * t20 - t71 * t34;
t6 = t27 * qJD(5) + t71 * t20 - t67 * t34;
t4 = (-t40 * t64 - t6) * t60 + (t28 * t64 - t21) * t59;
t3 = -t6 * t116 + (-t71 * t103 - t85) * t28;
t2 = t7 * t67 * t68 + (t72 * t117 - t127 * t71) * t26;
t1 = t4 * (t60 * t116 - t59 * t72) + ((-t64 * t71 + t45) * t68 * t59 + ((-t64 + t104) * t72 + t85) * t60) * t14;
t8 = [0, -0.2e1 * t70 * t89, qJD(2) ^ 2 * t70, 0, 0.2e1 * pkin(2) * t89, t33 * t48 - t46 * t38, t33 * t49 - t48 * t34 + t38 * t47 + t46 * t39, t38 * t66, -t39 * t66, 0, (qJD(1) * t39 + t34) * t61 + 0.2e1 * t47 * t95, (qJD(1) * t38 + t33) * t61 + (t46 - t125) * t95, t20 * t120 + t80 * t42, -t34 * t49 - t45 * t39, t28 * t12 + t6 * t35, t21 * t121 + t81 * t40, t22 * t12 + t9 * t35 + (t18 * t28 + t37 * t6 + (t37 * t114 - t17 * t48) * qJD(4)) * t72 + (t18 * t114 - t17 * t38 + t77 * t37 - t48 * t5) * t68, (-t4 * t35 + t14 * (-t64 * t121 - t12)) * t60 + (-t4 * t121 + t14 * (t35 * t64 - t81)) * t59, (t26 * t84 + t7 * t49) * t71 + (-t7 * t120 + t26 * t76) * t67; 0, t70 * t112, 0, 0, -pkin(2) * t112, t122, t29, t24, t25, 0, t92 + (-t47 * t108 + t69 * t87) * pkin(3), -t91 + (-t46 * t108 + t73 * t87) * pkin(3), t11, -t123, t3, t10, t68 * t62 * t6 + (-(-t62 * t101 - t69 * t110) * t40 - t63 * t21) * t67 + (t68 * t98 + t88 * t72) * t28 + (-t124 + (-t62 * t21 - t100) * t72 + (-qJD(5) * t63 + t88 * t68 - t72 * t98) * t40) * t71 + t79, t1, t2; 0, 0, 0, 0, 0, t122, t29, t24, t25, 0, t69 * t83 + t92, t73 * t83 - t91, t11, -t123, t3, t10, -(t72 * t36 + t68 * t96) * t28 + (-pkin(4) * t21 - t40 * t97) * t67 + (-t72 * t100 - t124 + (-pkin(4) * qJD(5) - t68 * t36 + t72 * t96) * t40) * t71 + ((-t40 * t104 - t6) * t68 + t77 * t72) * pkin(5) + t79, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * t41, -t34, t28 * t114 + t6 * t67, -t40 * t42, t17 * t42 - t23 * t28 + t9 * t67, -t4 * t60 * t67 + (-t60 * t114 + (t64 * t67 + t42) * t59) * t14, -t26 * t126 + t7 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t27, t21, -t22 * t27 - t78 + (qJD(5) - t40) * t82, -t4 * t59 + (-pkin(7) * t119 + (pkin(7) * t40 + t26) * t60) * t14, -t7 * pkin(7) + t26 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * (t60 * t40 - t119), t7;];
tauc_reg = t8;
