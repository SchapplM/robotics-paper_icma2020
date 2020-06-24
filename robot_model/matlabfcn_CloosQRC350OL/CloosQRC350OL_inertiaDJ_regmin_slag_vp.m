% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((6+1)*6/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = CloosQRC350OL_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:18
% EndTime: 2020-06-23 22:03:21
% DurationCPUTime: 1.49s
% Computational Cost: add. (734->108), mult. (2048->225), div. (0->0), fcn. (2062->10), ass. (0->93)
t105 = sin(qJ(3));
t106 = sin(qJ(2));
t58 = cos(qJ(3));
t59 = cos(qJ(2));
t35 = -t105 * t106 + t58 * t59;
t108 = qJD(2) + qJD(3);
t107 = t58 * pkin(3);
t47 = t106 * pkin(3) + pkin(2);
t34 = t105 * t59 + t58 * t106;
t19 = t108 * t34;
t104 = t35 * t19;
t54 = sin(qJ(4));
t50 = t54 ^ 2;
t103 = t35 * t50;
t102 = t35 * t54;
t53 = sin(qJ(5));
t101 = t53 * t54;
t55 = cos(qJ(6));
t100 = t53 * t55;
t57 = cos(qJ(4));
t99 = t53 * t57;
t98 = t54 * t19;
t56 = cos(qJ(5));
t97 = t54 * t56;
t96 = t56 * t57;
t36 = t53 * pkin(4) - pkin(5) * t96;
t89 = qJD(4) * t57;
t72 = t54 * t89;
t68 = t56 * t72;
t91 = qJD(4) * t54;
t94 = 0.2e1 * pkin(5) * t68 + t36 * t91;
t93 = qJD(2) * t59;
t92 = qJD(4) * t35;
t90 = qJD(4) * t56;
t88 = qJD(5) * t53;
t87 = qJD(5) * t54;
t86 = qJD(5) * t56;
t85 = qJD(5) * t57;
t84 = qJD(6) * t53;
t83 = qJD(6) * t55;
t82 = t35 * t96;
t78 = t105 * pkin(3);
t46 = t78 - pkin(5);
t71 = pkin(4) + t107;
t26 = t46 * t96 + t53 * t71;
t77 = t50 * t88;
t81 = t26 * t91 + (-0.2e1 * t68 + t77) * t46;
t80 = qJD(3) * t107;
t79 = pkin(3) * t93;
t76 = t55 * t86;
t52 = sin(qJ(6));
t75 = t52 * t83;
t74 = t53 * t89;
t73 = t53 * t86;
t70 = qJD(2) * t106;
t69 = t56 * t80;
t67 = qJD(3) * t78;
t20 = t35 * t108;
t65 = -t35 * t85 - t20;
t64 = t35 * t89 - t98;
t63 = -t53 * t87 + t56 * t89;
t62 = t53 * t85 + t54 * t90;
t30 = t54 * t86 + t74;
t61 = t52 * t84 - t76;
t60 = qJD(5) * t34 + t19 * t57 + t35 * t91;
t51 = t55 ^ 2;
t49 = t53 ^ 2;
t42 = -0.2e1 * t72;
t41 = 0.2e1 * t72;
t32 = t35 ^ 2;
t31 = t52 * t57 + t55 * t97;
t25 = pkin(4) * t86 + t62 * pkin(5);
t24 = 0.2e1 * (t68 - t77) * t56;
t23 = 0.2e1 * t49 * t72 + 0.2e1 * t50 * t73;
t22 = t49 * t87 + (-t87 * t56 - t74) * t56;
t18 = t34 * pkin(4) + t35 * pkin(5) + t47;
t17 = -t53 * t34 + t82;
t16 = -t56 * t34 - t35 * t99;
t14 = (qJD(6) + t90) * t57 * t55 + (-t55 * t88 + (-qJD(6) * t56 - qJD(4)) * t52) * t54;
t13 = t62 * t46 + t53 * t67 - t57 * t69 - t71 * t86;
t12 = t52 * t102 - t55 * t17;
t11 = 0.2e1 * t31 * t14;
t10 = t14 * t52 + t31 * t83;
t9 = t20 * pkin(4) - t19 * pkin(5) + t79;
t8 = t50 * t92 + (-t92 * t57 + t98) * t57;
t7 = -t14 * t100 + t61 * t31;
t6 = t60 * t53 + t65 * t56;
t5 = t65 * t53 - t60 * t56;
t4 = t6 * t101 + t30 * t16;
t3 = -t63 * t17 - t5 * t97;
t2 = (qJD(6) * t102 - t5) * t55 + (qJD(6) * t17 + t64) * t52;
t1 = t12 * t14 + t2 * t31;
t15 = [0, -0.2e1 * t59 * t70, 0, 0, 0.2e1 * pkin(2) * t93, -0.2e1 * t104, 0.2e1 * t19 * t34 - 0.2e1 * t35 * t20, 0, 0, 0, 0.2e1 * t47 * t20 + 0.2e1 * t34 * t79, -0.2e1 * t47 * t19 + 0.2e1 * t35 * t79, 0.2e1 * (-t57 * t104 - t32 * t91) * t57, 0.2e1 * t34 * t20, 0.2e1 * t17 * t5, -0.2e1 * t19 * t103 + 0.2e1 * t32 * t72, 0.2e1 * (t56 * t103 + t17 * t57) * t9 + 0.2e1 * (t5 * t57 + (-t19 * t56 - t35 * t88) * t50 + (-t17 + 0.2e1 * t82) * t91) * t18, 0.2e1 * t12 * t2, 0.2e1 * t16 * t6; 0, 0, -t70, 0, 0, 0, 0, -t19, -t20, 0, 0, 0, t8, 0, t3, -t8, (t17 * t46 - t26 * t35) * t89 + (t13 * t35 + t17 * t80 + t19 * t26 + t46 * t5) * t54, t1, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t67, -0.2e1 * t80, t41, 0, t24, t42, 0.2e1 * t13 * t57 - 0.2e1 * t50 * t69 + 0.2e1 * t81, t11, t23; 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, 0, 0, t8, 0, t3, -t8, (-pkin(5) * t17 - t35 * t36) * t89 + (-pkin(5) * t5 + t19 * t36 - t25 * t35) * t54, t1, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t80, t41, 0, t24, t42, (t13 - t25) * t57 + (-pkin(5) * t88 - t69) * t50 + t81 + t94, t11, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t24, t42, -0.2e1 * pkin(5) * t77 - 0.2e1 * t25 * t57 + 0.2e1 * t94, t11, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t17 * t86 + t5 * t53, 0, t9 * t99 + (-t53 * t91 + t56 * t85) * t18, -t2 * t100 + t61 * t12, -t16 * t88 + t6 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t80 * t101 + t30 * t46, t7, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t30 * pkin(5), t7, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t73, 0, 0, -0.2e1 * t49 * t75 + 0.2e1 * t51 * t73, -0.2e1 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t63 * t18 + t9 * t97, t12 * t83 + t2 * t52, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t13, t10, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t25, t10, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t76 + (t52 ^ 2 - t51) * t84, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t75, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
