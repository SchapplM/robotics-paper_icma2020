% Calculate inertial parameters regressor of joint inertia matrix for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = CloosQRC350OL_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:19
% EndTime: 2020-06-23 22:03:25
% DurationCPUTime: 3.09s
% Computational Cost: add. (1046->121), mult. (2299->246), div. (0->0), fcn. (2611->10), ass. (0->120)
t81 = sin(qJ(3));
t82 = sin(qJ(2));
t86 = cos(qJ(3));
t87 = cos(qJ(2));
t48 = t81 * t87 + t86 * t82;
t50 = -t81 * t82 + t86 * t87;
t65 = t82 * pkin(3) + pkin(2);
t26 = t48 * pkin(4) + t50 * pkin(5) + t65;
t128 = 0.2e1 * t26;
t127 = 0.2e1 * t65;
t80 = sin(qJ(4));
t84 = cos(qJ(5));
t110 = t80 * t84;
t78 = sin(qJ(6));
t83 = cos(qJ(6));
t85 = cos(qJ(4));
t45 = t83 * t110 + t78 * t85;
t126 = pkin(6) * t45;
t79 = sin(qJ(5));
t125 = pkin(6) * t79;
t124 = pkin(6) * t84;
t74 = t80 ^ 2;
t68 = t74 * pkin(5);
t123 = t81 * pkin(3);
t111 = t80 * t50;
t109 = t84 * t85;
t22 = t50 * t109 - t79 * t48;
t9 = t83 * t111 + t78 * t22;
t122 = t9 * t83;
t10 = t78 * t111 - t83 * t22;
t121 = t10 * t78;
t71 = t86 * pkin(3);
t64 = t71 + pkin(4);
t114 = t79 * t64;
t63 = -pkin(5) + t123;
t28 = t114 + (t63 * t84 - pkin(6)) * t85;
t39 = (t63 - t124) * t80;
t15 = -t83 * t28 + t78 * t39;
t120 = t15 * t78;
t66 = t79 * pkin(4);
t40 = t66 + (-pkin(5) * t84 - pkin(6)) * t85;
t55 = (-pkin(5) - t124) * t80;
t24 = -t83 * t40 + t78 * t55;
t119 = t24 * t78;
t118 = t26 * t80;
t44 = -t78 * t110 + t83 * t85;
t37 = t44 * t83;
t38 = t45 * t78;
t25 = t26 ^ 2;
t117 = t74 * t25;
t116 = t74 * t84;
t76 = t84 ^ 2;
t62 = t76 * t74;
t115 = t78 * t79;
t61 = t79 * t80;
t113 = t79 * t83;
t112 = t79 * t85;
t108 = t85 * t26;
t107 = (t37 + t38) * pkin(6);
t77 = t85 ^ 2;
t106 = t77 * pkin(5) + t68;
t72 = t78 ^ 2;
t75 = t83 ^ 2;
t105 = t72 + t75;
t104 = t74 + t77;
t103 = t22 * t110;
t102 = t26 * t61;
t33 = t44 * t115;
t101 = t45 * t113;
t100 = t63 * t116;
t99 = t78 * t113;
t11 = (-pkin(6) * t50 - t26 * t84) * t80;
t4 = t22 * pkin(6) + t108;
t2 = t78 * t11 + t83 * t4;
t3 = -t83 * t11 + t78 * t4;
t98 = -t2 * t45 + t3 * t44;
t97 = t105 * pkin(6) ^ 2;
t96 = t104 * t63;
t95 = 0.2e1 * t105 * pkin(6);
t94 = t2 * t83 + t3 * t78;
t93 = -t2 * t78 + t3 * t83;
t14 = t78 * t28 + t83 * t39;
t92 = -t14 * t78 + t15 * t83;
t23 = t78 * t40 + t83 * t55;
t91 = -t23 * t78 + t24 * t83;
t89 = pkin(5) ^ 2;
t73 = t79 ^ 2;
t67 = t74 * t89;
t60 = t73 * t74;
t59 = t63 ^ 2;
t58 = pkin(5) * t116;
t57 = t79 * t110;
t56 = t74 * t59;
t53 = pkin(5) * t109 - t66;
t52 = t84 * pkin(4) + pkin(5) * t112;
t51 = t52 ^ 2;
t47 = t50 ^ 2;
t46 = t48 ^ 2;
t43 = t45 ^ 2;
t42 = t44 ^ 2;
t41 = t74 * t47;
t36 = t85 * t111;
t32 = t63 * t109 + t114;
t31 = -t63 * t112 + t84 * t64;
t30 = pkin(6) * t33;
t29 = t31 ^ 2;
t27 = t31 * t52;
t21 = -t50 * t112 - t84 * t48;
t20 = t21 ^ 2;
t19 = t77 * t25;
t18 = t21 * t84;
t17 = t73 * t117;
t16 = t21 * t61;
t13 = t24 * t44;
t12 = t52 * t102;
t8 = t31 * t102;
t7 = t15 * t44;
t6 = t10 * t45;
t5 = t9 * t44;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t87 ^ 2, 0, 0, t82 ^ 2, 0, 0, 0.2e1 * pkin(2) * t82, 0, 0, pkin(2) ^ 2, t47, -0.2e1 * t50 * t48, 0, t46, 0, 0, t48 * t127, t50 * t127, 0, t65 ^ 2, t77 * t47, 0, 0, t41, 0, t46, 0, 0, t104 * t50 * t128, t19 + t117, t22 ^ 2, 0, 0, t20, 0, t41, 0, (t50 * t116 + t22 * t85) * t128, 0, t25 * t62 + t17 + t19, t10 ^ 2, 0, 0, t9 ^ 2, 0, t20, 0, 0, -0.2e1 * t2 * t10 + 0.2e1 * t3 * t9, t2 ^ 2 + t3 ^ 2 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t48, 0, 0, 0, (-t48 * t81 - t50 * t86) * pkin(3), 0, -t36, 0, 0, t36, 0, 0, 0, 0, 0, 0, -t103, 0, 0, t16, 0, t36, 0, (t22 * t63 - t32 * t50) * t80, 0, t8 + (-t32 * t84 + t63 * t85) * t118, t6, 0, 0, t5, 0, t16, 0, 0, -t14 * t10 + t15 * t9 + t98, t2 * t14 + t3 * t15 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t123, 0, (t81 ^ 2 + t86 ^ 2) * pkin(3) ^ 2, t74, 0, 0, t77, 0, 0, 0, 0, -0.2e1 * t96, t77 * t59 + t64 ^ 2 + t56, t62, 0, 0, t60, 0, t77, 0, -0.2e1 * t32 * t85 - 0.2e1 * t100, 0, t32 ^ 2 + t29 + t56, t43, 0, 0, t42, 0, t60, 0, 0, -0.2e1 * t14 * t45 + 0.2e1 * t7, t14 ^ 2 + t15 ^ 2 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t48, 0, 0, 0, 0, 0, -t36, 0, 0, t36, 0, 0, 0, 0, 0, 0, -t103, 0, 0, t16, 0, t36, 0, (-pkin(5) * t22 + t50 * t53) * t80, 0, t12 + (-pkin(5) * t85 + t53 * t84) * t118, t6, 0, 0, t5, 0, t16, 0, 0, -t23 * t10 + t24 * t9 + t98, t2 * t23 + t3 * t24 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t71, -t123, 0, 0, t74, 0, 0, t77, 0, 0, 0, 0, -t96 + t106, t64 * pkin(4) - pkin(5) * t96, t62, 0, 0, t60, 0, t77, 0, -t100 + t58 + (-t32 + t53) * t85, 0, -t32 * t53 - t63 * t68 + t27, t43, 0, 0, t42, 0, t60, 0, 0, t13 + t7 + (-t14 - t23) * t45, t14 * t23 + t15 * t24 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t74, 0, 0, t77, 0, 0, 0, 0, 0.2e1 * t106, pkin(4) ^ 2 + t77 * t89 + t67, t62, 0, 0, t60, 0, t77, 0, 0.2e1 * t53 * t85 + 0.2e1 * t58, 0, t53 ^ 2 + t51 + t67, t43, 0, 0, t42, 0, t60, 0, 0, -0.2e1 * t23 * t45 + 0.2e1 * t13, t23 ^ 2 + t24 ^ 2 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, 0, 0, 0, t22 * t79, 0, 0, t18, 0, 0, 0, t79 * t108, 0, 0, -t10 * t113, 0, 0, t9 * t115, 0, t18, 0, 0, ((-t10 * t83 + t78 * t9) * pkin(6) + t94) * t79, t94 * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, t57, 0, 0, 0, t63 * t61, 0, 0, -t101, 0, 0, t33, 0, t57, 0, 0, t30 + (t120 + (t14 - t126) * t83) * t79, (t14 * t83 + t120) * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, t57, 0, 0, 0, -pkin(5) * t61, 0, 0, -t101, 0, 0, t33, 0, t57, 0, 0, t30 + (t119 + (t23 - t126) * t83) * t79, (t23 * t83 + t119) * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t73, 0, 0, t76, 0, 0, 0, 0, 0, 0, t75 * t73, 0, 0, t72 * t73, 0, t76, 0, 0, t73 * t95, t73 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, 0, t26 * t110, 0, 0, t121, 0, 0, t122, 0, 0, 0, 0, (t121 + t122) * pkin(6) + t93, t93 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, -t32, 0, 0, t38, 0, 0, t37, 0, 0, 0, 0, t92 + t107, t92 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, t53, 0, 0, t38, 0, 0, t37, 0, 0, 0, 0, t91 + t107, t91 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, 0, 0, t99, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t72, 0, 0, t75, 0, 0, 0, 0, t95, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
