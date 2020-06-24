% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((6+1)*6/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = CloosQRC350DE_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:32
% EndTime: 2020-06-23 21:07:35
% DurationCPUTime: 1.62s
% Computational Cost: add. (915->109), mult. (2262->228), div. (0->0), fcn. (2107->10), ass. (0->95)
t106 = cos(qJ(3));
t107 = cos(qJ(2));
t58 = sin(qJ(3));
t59 = sin(qJ(2));
t35 = -t107 * t106 + t59 * t58;
t109 = qJD(2) + qJD(3);
t108 = pkin(3) * t58;
t36 = t106 * t59 + t107 * t58;
t19 = t109 * t36;
t105 = t35 * t19;
t57 = sin(qJ(4));
t55 = t57 ^ 2;
t104 = t35 * t55;
t53 = pkin(7) * qJ(5) - qJ(6);
t48 = sin(t53);
t103 = t48 * t57;
t49 = cos(t53);
t52 = pkin(7) * qJD(5) - qJD(6);
t102 = t49 * t52;
t56 = sin(qJ(5));
t101 = t49 * t56;
t100 = t52 * t56;
t99 = t56 * t57;
t61 = cos(qJ(4));
t98 = t56 * t61;
t97 = t57 * t19;
t60 = cos(qJ(5));
t96 = t57 * t60;
t94 = t60 * t61;
t37 = t56 * pkin(4) - pkin(5) * t94;
t89 = qJD(4) * t61;
t75 = t57 * t89;
t73 = t60 * t75;
t91 = qJD(4) * t57;
t93 = 0.2e1 * pkin(5) * t73 + t37 * t91;
t92 = qJD(4) * t35;
t90 = qJD(4) * t60;
t88 = qJD(5) * t56;
t87 = qJD(5) * t57;
t86 = qJD(5) * t60;
t85 = qJD(5) * t61;
t84 = t48 * t102;
t83 = t35 * t94;
t51 = -pkin(5) + t108;
t80 = pkin(3) * t106;
t70 = t80 + pkin(4);
t27 = t51 * t94 + t56 * t70;
t79 = t55 * t88;
t82 = t27 * t91 + (-0.2e1 * t73 + t79) * t51;
t81 = qJD(3) * t108;
t78 = t56 * t87;
t77 = t56 * t89;
t76 = t56 * t86;
t50 = -t59 * pkin(3) - pkin(2);
t74 = t107 * qJD(2);
t72 = qJD(3) * t80;
t71 = pkin(3) * t74;
t20 = t35 * t109;
t68 = -t35 * t85 - t20;
t67 = t60 * t72;
t66 = t35 * t89 + t97;
t65 = -t48 * t100 + t49 * t86;
t64 = t60 * t89 - t78;
t63 = t56 * t85 + t57 * t90;
t32 = t57 * t86 + t77;
t62 = -qJD(5) * t36 - t19 * t61 + t35 * t91;
t54 = t56 ^ 2;
t47 = t49 ^ 2;
t43 = -0.2e1 * t75;
t42 = 0.2e1 * t75;
t33 = t35 ^ 2;
t28 = t32 * pkin(7);
t26 = -t48 * t61 + t49 * t96;
t25 = pkin(4) * t86 + t63 * pkin(5);
t24 = 0.2e1 * (t73 - t79) * t60;
t23 = 0.2e1 * t54 * t75 + 0.2e1 * t55 * t76;
t22 = t54 * t87 + (-t87 * t60 - t77) * t60;
t18 = -t36 * pkin(4) + t35 * pkin(5) + t50;
t17 = -t35 * t98 + t60 * t36;
t16 = t56 * t36 + t83;
t14 = t63 * t51 + t56 * t81 - t61 * t67 - t70 * t86;
t13 = -t35 * t103 - t49 * t16;
t12 = t20 * pkin(4) + t19 * pkin(5) - t71;
t11 = t55 * t92 + (-t92 * t61 - t97) * t61;
t10 = (-t52 * t60 + qJD(4)) * t103 + (-t78 + (-t52 + t90) * t61) * t49;
t9 = 0.2e1 * t26 * t10;
t8 = -t10 * t48 - t26 * t102;
t7 = t62 * t56 + t68 * t60;
t6 = t68 * t56 - t62 * t60;
t5 = -t10 * t101 - t65 * t26;
t4 = t32 * t17 + t7 * t99;
t3 = -t64 * t16 - t6 * t96;
t2 = (-t35 * t52 * t57 - t6) * t49 + (t16 * t52 - t66) * t48;
t1 = t13 * t10 + t2 * t26;
t15 = [0, -0.2e1 * t59 * t74, 0, 0, 0.2e1 * pkin(2) * t74, 0.2e1 * t105, 0.2e1 * t19 * t36 - 0.2e1 * t35 * t20, 0, 0, 0, 0.2e1 * t50 * t20 + 0.2e1 * t36 * t71, 0.2e1 * t50 * t19 - 0.2e1 * t35 * t71, 0.2e1 * (t61 * t105 - t33 * t91) * t61, -0.2e1 * t36 * t20, 0.2e1 * t16 * t6, 0.2e1 * t19 * t104 + 0.2e1 * t33 * t75, 0.2e1 * (t60 * t104 + t16 * t61) * t12 + 0.2e1 * (t6 * t61 + (t19 * t60 - t35 * t88) * t55 + (-t16 + 0.2e1 * t83) * t91) * t18, 0.2e1 * t13 * t2, 0.2e1 * t17 * t7; 0, 0, qJD(2) * t59, 0, 0, 0, 0, t19, -t20, 0, 0, 0, t11, 0, t3, -t11, (t16 * t51 - t27 * t35) * t89 + (t14 * t35 + t16 * t72 - t19 * t27 + t51 * t6) * t57, t1, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t81, -0.2e1 * t72, t42, 0, t24, t43, 0.2e1 * t14 * t61 - 0.2e1 * t55 * t67 + 0.2e1 * t82, t9, t23; 0, 0, 0, 0, 0, 0, 0, t19, -t20, 0, 0, 0, t11, 0, t3, -t11, (-pkin(5) * t16 - t35 * t37) * t89 + (-pkin(5) * t6 - t19 * t37 - t25 * t35) * t57, t1, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t72, t42, 0, t24, t43, (t14 - t25) * t61 + (-pkin(5) * t88 - t67) * t55 + t82 + t93, t9, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t24, t43, -0.2e1 * pkin(5) * t79 - 0.2e1 * t25 * t61 + 0.2e1 * t93, t9, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t16 * t86 + t6 * t56, 0, t12 * t98 + (-t56 * t91 + t60 * t85) * t18, -t2 * t101 - t65 * t13, -t17 * t88 + t7 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t32 * t51 + t72 * t99, t5, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t32 * pkin(5), t5, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t76, 0, 0, 0.2e1 * t47 * t76 - 0.2e1 * t54 * t84, -0.2e1 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t12 * t96 + t64 * t18, -t13 * t102 - t2 * t48, -t7 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t14, t8, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t25, t8, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t100 + t48 * t65, pkin(7) * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t84, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
