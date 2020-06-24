% Calculate joint inertia matrix for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350DE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350DE_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:02:06
% EndTime: 2020-06-23 21:02:08
% DurationCPUTime: 1.47s
% Computational Cost: add. (652->187), mult. (864->285), div. (0->0), fcn. (462->14), ass. (0->122)
t119 = Icges(7,3) + Icges(6,2);
t80 = pkin(6) + rSges(7,3);
t124 = m(6) * rSges(6,2) ^ 2 + t80 ^ 2 * m(7);
t100 = -t119 + t124;
t151 = Icges(6,1) + Icges(7,2);
t30 = t100 + t151;
t75 = cos(qJ(5));
t65 = t75 ^ 2;
t143 = t30 * t65;
t154 = Icges(5,1) + t143;
t53 = pkin(7) * qJ(5) - qJ(6);
t43 = cos(t53);
t41 = t43 ^ 2;
t71 = sin(qJ(4));
t132 = t71 * t41;
t42 = sin(t53);
t139 = t42 * t43;
t69 = Icges(7,2) - Icges(7,1);
t76 = cos(qJ(4));
t153 = (t75 * t132 - t76 * t139) * t69;
t120 = Icges(5,2) + Icges(6,3);
t14 = -Icges(7,1) - t100 - t120 + t154;
t152 = t143 - Icges(5,3) - t124 - t151;
t73 = sin(qJ(2));
t44 = t73 * pkin(3) + pkin(2);
t68 = qJ(2) + qJ(3);
t55 = sin(t68);
t107 = -t55 * pkin(4) - t44;
t70 = sin(qJ(5));
t127 = t76 * t70;
t56 = cos(t68);
t150 = (t55 * t127 - t56 * t75) * rSges(6,2) - t56 * pkin(5) + t107;
t148 = m(5) + m(6);
t147 = m(5) * rSges(5,3);
t146 = rSges(4,3) * m(4);
t72 = sin(qJ(3));
t145 = t72 * pkin(4);
t77 = cos(qJ(3));
t78 = cos(qJ(2));
t29 = -t73 * t72 + t78 * t77;
t144 = t29 * t76;
t142 = t30 * t75;
t39 = m(6) * rSges(6,2) + t80 * m(7);
t141 = t39 * t75;
t140 = t39 * t77;
t59 = pkin(3) * t72;
t45 = t59 - pkin(5);
t138 = t45 * t75;
t137 = t56 * t71;
t136 = (t75 + 0.1e1) * (t75 - 0.1e1);
t135 = t69 * t41;
t134 = t70 * t71;
t133 = t70 * t75;
t131 = t72 * t75;
t130 = t72 * t76;
t47 = rSges(4,2) * t146 - Icges(4,6);
t48 = rSges(4,1) * t146 - Icges(4,5);
t129 = t73 * (t47 * t77 + t72 * t48);
t128 = t75 * t76;
t126 = t77 * t76;
t122 = Icges(7,3) * t71;
t125 = pkin(4) * t39 + pkin(7) * t122;
t86 = pkin(5) ^ 2;
t87 = pkin(4) ^ 2;
t63 = t86 + t87;
t123 = Icges(7,3) * pkin(7);
t121 = t75 * Icges(7,3);
t118 = m(7) + t148;
t117 = 0.2e1 * t39;
t116 = t29 * t128;
t115 = t69 * t139;
t114 = t71 * t139;
t74 = sin(qJ(1));
t113 = t74 * t134;
t112 = t69 * t132;
t79 = cos(qJ(1));
t111 = t79 * t134;
t110 = -t127 / 0.2e1;
t109 = -pkin(5) * t72 + pkin(3);
t108 = -rSges(3,1) * t73 - pkin(2);
t106 = t80 * t75 + pkin(5);
t105 = 0.2e1 * t115;
t104 = -t47 * t72 + t48 * t77;
t66 = t76 ^ 2;
t103 = (t66 - 0.1e1 / 0.2e1) * t75;
t102 = t75 * t114;
t101 = Icges(7,1) + Icges(6,3) + t124;
t96 = -rSges(4,1) * t55 - rSges(4,2) * t56 - t44;
t95 = (-(t118 * pkin(5) + m(4) * rSges(4,2) + t147) * t72 + (t118 * pkin(4) + m(4) * rSges(4,1)) * t77) * pkin(3);
t93 = (-pkin(5) - rSges(5,3)) * t56 + t107;
t92 = (t76 * t41 + t102) * t69;
t54 = t65 + 0.1e1;
t90 = (t54 * t66 - t65) * t135 + Icges(4,3) - t14 * t66 + t119 + 0.2e1 * t76 * t69 * t102 + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + 0.2e1 * pkin(5) * t147 + (rSges(5,3) ^ 2 + t87) * m(5) + t154;
t89 = t86 * m(5) + t63 * m(6) + t90;
t88 = pkin(3) ^ 2;
t60 = pkin(3) * t77;
t50 = t63 * m(7);
t40 = t70 * t122;
t36 = -t79 * rSges(2,1) - t74 * rSges(2,2);
t35 = -t74 * rSges(2,1) + t79 * rSges(2,2);
t34 = t80 * t127 - pkin(4);
t32 = pkin(4) * t140;
t28 = t78 * t72 + t73 * t77;
t24 = -t74 * rSges(3,3) + t108 * t79;
t23 = t79 * rSges(3,3) + t108 * t74;
t21 = t56 * t128 - t55 * t70;
t19 = t56 * t127 + t55 * t75;
t16 = pkin(5) * t39 + t142;
t15 = (t77 * pkin(4) + t109) * t73 + (t77 * pkin(5) + t145) * t78 + pkin(2);
t13 = -t74 * rSges(4,3) + t96 * t79;
t12 = t79 * rSges(4,3) + t96 * t74;
t10 = t93 * t79;
t9 = t93 * t74;
t8 = -t42 * t137 - t21 * t43;
t7 = t43 * t137 - t21 * t42;
t6 = t39 * t145 + t16 * t77;
t5 = (-t106 * t72 - t34 * t77 + pkin(3)) * t73 + (t106 * t77 - t34 * t72) * t78 + pkin(2);
t4 = -rSges(6,2) * t113 + t150 * t79;
t3 = rSges(6,2) * t111 + t150 * t74;
t2 = t80 * t111 - t74 * t5;
t1 = -t80 * t113 - t5 * t79;
t11 = [Icges(3,1) * t78 ^ 2 + Icges(6,1) * t21 ^ 2 + Icges(7,1) * t8 ^ 2 + t73 ^ 2 * Icges(3,2) + Icges(7,2) * t7 ^ 2 + Icges(2,3) + t119 * t19 ^ 2 + (-0.2e1 * Icges(4,4) * t56 + (Icges(5,3) + Icges(4,2)) * t55) * t55 + m(7) * (t1 ^ 2 + t2 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(3) * (t23 ^ 2 + t24 ^ 2) + m(2) * (t35 ^ 2 + t36 ^ 2) + (t120 * t71 ^ 2 + t66 * Icges(5,1) + Icges(4,1)) * t56 ^ 2; -((t54 * t126 - t70 * t131) * t78 - t73 * (t54 * t130 + t77 * t133)) * t112 + ((t77 * t103 + t72 * t110) * t78 - t73 * (t70 * t126 + (0.2e1 * t66 - 0.1e1) * t131) / 0.2e1) * t105 + ((t14 * t126 + t70 * (t39 * t109 - t30 * t131 + t32)) * t78 - t73 * (t14 * t130 + t70 * t6)) * t71 + (rSges(3,3) * rSges(3,1) * m(3) + pkin(3) * t146 - Icges(3,5) + t104) * t78 - t129; (t88 + t63) * m(7) + (m(4) + t148) * t88 + (-(t60 + pkin(4)) * t127 - t138) * t117 + 0.2e1 * t95 + m(3) * rSges(3,1) ^ 2 + t89 + Icges(3,3); -(-t28 * t133 + t54 * t144) * t112 + (t29 * t103 + t28 * t110) * t105 + (t14 * t144 + t70 * ((-t16 * t72 + t32) * t78 - t73 * t6)) * t71 + t104 * t78 - t129; t50 + t89 + t95 + (-(t60 + 0.2e1 * pkin(4)) * t127 - (t59 - 0.2e1 * pkin(5)) * t75) * t39; t148 * t86 + t87 * m(6) + (-pkin(4) * t127 + pkin(5) * t75) * t117 + t90 + t50; -(t29 * t142 + t15 * t39) * t127 - t152 * t28 + ((t70 * t116 + t28 * t136) * t41 + t70 * t29 * t114) * t69; ((t39 * t45 - t142) * t71 + t153) * t70; -(t71 * t16 - t153) * t70; t135 * t136 - t152; (-t101 * t29 - t15 * t141) * t71 - (t29 * t127 + t28 * t75) * t123 + (-t29 * t132 + (-t28 * t70 + t116) * t139) * t69; (-t39 * t138 + t101) * t76 - (pkin(3) * t140 + t125) * t70 + t92; (pkin(5) * t141 + t101) * t76 - t125 * t70 + t92; -pkin(7) * t121 - t70 * t115; Icges(7,3) * pkin(7) ^ 2 + t101 + t135; Icges(7,3) * t19; t40; t40; t121; -t123; Icges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11), t11(16); t11(2), t11(3), t11(5), t11(8), t11(12), t11(17); t11(4), t11(5), t11(6), t11(9), t11(13), t11(18); t11(7), t11(8), t11(9), t11(10), t11(14), t11(19); t11(11), t11(12), t11(13), t11(14), t11(15), t11(20); t11(16), t11(17), t11(18), t11(19), t11(20), t11(21);];
Mq = res;
