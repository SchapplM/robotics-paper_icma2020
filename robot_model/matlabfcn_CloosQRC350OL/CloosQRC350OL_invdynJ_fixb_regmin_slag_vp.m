% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% CloosQRC350OL
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% tau_reg [6x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = CloosQRC350OL_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:18
% EndTime: 2020-06-23 22:03:23
% DurationCPUTime: 2.36s
% Computational Cost: add. (2093->200), mult. (4679->315), div. (0->0), fcn. (3785->12), ass. (0->124)
t84 = sin(qJ(4));
t124 = qJD(5) * t84;
t88 = cos(qJ(5));
t89 = cos(qJ(4));
t136 = t88 * t89;
t86 = sin(qJ(2));
t133 = qJD(1) * t86;
t85 = sin(qJ(3));
t112 = t85 * t133;
t91 = cos(qJ(2));
t132 = qJD(1) * t91;
t90 = cos(qJ(3));
t51 = -t90 * t132 + t112;
t54 = t85 * t91 + t90 * t86;
t52 = t54 * qJD(1);
t83 = sin(qJ(5));
t34 = -t52 * t136 + t83 * t51;
t150 = t83 * t124 - t34;
t118 = qJD(1) * qJD(2);
t121 = qJDD(1) * t86;
t149 = t91 * t118 + t121;
t78 = qJD(2) + qJD(3);
t43 = -t89 * t51 - t84 * t78;
t49 = qJD(4) - t52;
t30 = -t83 * t43 + t88 * t49;
t79 = qJ(2) + qJ(3);
t75 = cos(t79);
t67 = g(3) * t75;
t119 = t91 * qJDD(1);
t39 = t78 * t54;
t28 = -t39 * qJD(1) + t90 * t119 - t85 * t121;
t105 = t78 * t91;
t111 = t86 * t118;
t99 = -qJD(3) * t112 + (-t111 + t119) * t85;
t29 = (qJD(1) * t105 + t121) * t90 + t99;
t50 = qJDD(1) * pkin(2) + t149 * pkin(3);
t13 = t29 * pkin(4) + t28 * pkin(5) + t50;
t58 = qJD(1) * pkin(2) + pkin(3) * t133;
t33 = t52 * pkin(4) - t51 * pkin(5) + t58;
t131 = qJD(2) * t85;
t116 = pkin(3) * t131;
t56 = -t78 * pkin(5) + t116;
t25 = -t84 * t33 + t89 * t56;
t120 = qJDD(2) * t85;
t129 = qJD(3) * t90;
t77 = qJDD(2) + qJDD(3);
t47 = -t77 * pkin(5) + (qJD(2) * t129 + t120) * pkin(3);
t7 = t25 * qJD(4) + t89 * t13 + t84 * t47;
t148 = t7 * t84;
t42 = t84 * t51 - t89 * t78;
t16 = t42 * qJD(4) + t89 * t28 - t84 * t77;
t93 = (-t78 * t132 - t121) * t90 - t99;
t26 = qJDD(4) + t93;
t9 = t30 * qJD(5) + t88 * t16 + t83 * t26;
t147 = t84 * t9;
t146 = t85 * pkin(3);
t145 = t90 * pkin(3);
t70 = t86 * pkin(3) + pkin(2);
t41 = qJD(5) - t42;
t144 = t41 * t88;
t143 = t49 * t51;
t142 = t51 * t52;
t55 = -t85 * t86 + t90 * t91;
t141 = t55 * t84;
t140 = t55 * t89;
t138 = t83 * t89;
t87 = cos(qJ(6));
t137 = t87 * t88;
t135 = t91 * qJD(1) ^ 2;
t134 = pkin(3) * qJD(2);
t130 = qJD(3) * t85;
t128 = qJD(4) * t84;
t127 = qJD(4) * t88;
t126 = qJD(4) * t89;
t125 = qJD(5) * t83;
t123 = qJD(5) * t89;
t24 = t89 * t33 + t84 * t56;
t122 = t24 * qJD(4);
t27 = qJD(6) + t30;
t117 = pkin(3) * t129;
t115 = t90 * t134;
t114 = t91 * t134;
t37 = -t51 * pkin(4) - t52 * pkin(5);
t109 = t58 * t52 + t67;
t68 = -pkin(5) + t146;
t108 = -pkin(3) * t132 + qJD(4) * t68 - t37;
t107 = qJD(2) * (-qJD(3) + t78);
t106 = qJD(3) * (-qJD(2) - t78);
t72 = qJDD(2) * t145;
t74 = sin(t79);
t104 = g(3) * t74 + t58 * t51 + t72;
t40 = t90 * t105 + (-t130 - t131) * t86;
t103 = -t55 * t123 - t40;
t102 = t83 * (t77 * pkin(4) - qJD(3) * t116 + t72) + t88 * (-t84 * t13 + t89 * t47 - t122);
t57 = t78 * pkin(4) + t115;
t101 = t83 * t25 - t88 * t57;
t31 = t88 * t43 + t83 * t49;
t82 = sin(qJ(6));
t100 = t82 * t31 + t87 * t41;
t98 = t55 * t126 - t39 * t84;
t97 = -t55 * t128 - t39 * t89;
t15 = t43 * qJD(4) + t84 * t28 + t89 * t77 + qJDD(5);
t96 = -qJD(4) * t31 - t41 * t125 + t88 * t15;
t95 = qJD(5) * t54 - t97;
t19 = t88 * t25 + t83 * t57;
t5 = -t101 * qJD(5) + t102;
t94 = -g(3) * (t74 * t138 - t75 * t88) - t5 * t89 + t150 * t24 + (-t84 * t52 + t128) * t19;
t69 = pkin(4) + t145;
t38 = t54 * pkin(4) + t55 * pkin(5) + t70;
t36 = t55 * t136 - t83 * t54;
t32 = t51 ^ 2 - t52 ^ 2;
t22 = t40 * pkin(4) - t39 * pkin(5) + t114;
t21 = -t51 * t78 + t93;
t20 = t52 * t78 + t28;
t17 = -t87 * t31 + t82 * t41;
t12 = t103 * t83 - t95 * t88;
t11 = -t49 * t89 * t43 - t16 * t84;
t10 = -t49 * t84 * t41 + t15 * t89;
t8 = -t31 * qJD(5) - t83 * t16 + t88 * t26 + qJDD(6);
t4 = t100 * qJD(6) + t82 * t15 - t87 * t9;
t3 = -t88 * t147 + (-t88 * t126 + t150) * t31;
t2 = t8 * t83 * t84 + ((-t51 + t124) * t88 + t49 * t138) * t27;
t1 = t4 * (t84 * t137 + t82 * t89) + ((t34 + (qJD(6) + t127) * t89) * t87 + (-t87 * t125 + (-qJD(6) * t88 - t49) * t82) * t84) * t17;
t6 = [qJDD(1), (-0.2e1 * t111 + t119) * t91, -qJD(2) ^ 2 * t86 + t91 * qJDD(2), 0, 0.2e1 * t149 * pkin(2), t28 * t55 + t51 * t39, -t28 * t54 - t55 * t29 + t39 * t52 + t51 * t40, -t39 * t78 + t55 * t77, -t40 * t78 - t54 * t77, 0, t52 * t114 + t70 * t29 + t58 * t40 + t50 * t54, -t51 * t114 + t70 * t28 - t58 * t39 + t50 * t55, t16 * t140 + t97 * t43, -t26 * t54 - t49 * t40, t31 * t12 + t9 * t36, t15 * t141 + t98 * t41, t24 * t12 + t7 * t36 + (t22 * t31 + t38 * t9 + (t38 * t144 - t19 * t55) * qJD(4)) * t89 + (t22 * t144 + t19 * t39 + t96 * t38 - t5 * t55) * t84, (-t4 * t36 + t17 * (qJD(6) * t141 - t12)) * t87 + (t4 * t141 + t17 * (qJD(6) * t36 + t98)) * t82, (t27 * t103 - t8 * t54) * t88 + (-t8 * t140 + t27 * t95) * t83; 0, t86 * t135, t119, qJDD(2), -pkin(2) * t135 + g(3) * t86, -t142, t32, t20, t21, t77, (t85 * t106 - t52 * t132 + t77 * t90) * pkin(3) + t104, (t51 * t132 + (-qJDD(2) - t77) * t85 + t90 * t106) * pkin(3) + t109, t11, -t143, t3, t10, t68 * t147 + (-(-pkin(3) * t130 - t68 * t123) * t41 - t69 * t15) * t83 + (t108 * t89 + t84 * t117) * t31 + (-t148 + (-t68 * t15 - t122) * t89 + (-qJD(5) * t69 + t108 * t84 - t89 * t117) * t41) * t88 + t94, t1, t2; 0, 0, 0, 0, 0, -t142, t32, t20, t21, t77, t107 * t146 + t104, (t90 * t107 - t120) * pkin(3) + t109, t11, -t143, t3, t10, -(t84 * t115 + t89 * t37) * t31 + (-pkin(4) * t15 - t41 * t116) * t83 + (-t89 * t122 - t148 + (-pkin(4) * qJD(5) + t89 * t115 - t84 * t37) * t41) * t88 + ((-t41 * t127 - t9) * t84 + t96 * t89) * pkin(5) + t94, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t42, t26, t31 * t144 + t9 * t83, -t41 * t43, t19 * t43 - t25 * t31 + (-t84 * t67 + t7) * t83, -t4 * t87 * t83 + (-t41 * t137 + (qJD(6) * t83 - t43) * t82) * t17, -t41 * t83 * t27 + t8 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t30, t15, -t24 * t30 - g(3) * (-t75 * t136 + t74 * t83) - t102 + (-t41 + qJD(5)) * t101, t27 * t87 * t17 + t4 * t82, t27 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t100, t8;];
tau_reg = t6;
