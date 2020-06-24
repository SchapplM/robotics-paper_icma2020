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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350DE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:01:52
% EndTime: 2020-06-23 21:01:55
% DurationCPUTime: 3.17s
% Computational Cost: add. (826->232), mult. (1047->289), div. (0->0), fcn. (473->10), ass. (0->141)
t73 = m(5) + m(6);
t168 = m(7) + t73;
t159 = pkin(5) * t168;
t104 = mrSges(5,3) + t159;
t27 = mrSges(4,2) + t104;
t128 = (Ifges(7,1) + Ifges(6,3));
t171 = Ifges(5,2) + t128;
t127 = (-Ifges(6,2) - Ifges(7,3));
t129 = (Ifges(6,1) + Ifges(7,2));
t157 = (mrSges(7,3) * pkin(6));
t56 = 2 * t157;
t76 = pkin(6) ^ 2;
t74 = (m(7) * t76);
t94 = t56 + t74 + t129;
t20 = t94 + t127;
t69 = cos(qJ(5));
t60 = t69 ^ 2;
t148 = t20 * t60;
t170 = Ifges(5,3) - t148;
t169 = -Ifges(5,1) + t127 - t148;
t40 = pkin(7) * qJ(5) - qJ(6);
t32 = cos(t40);
t145 = sin(t40) * t32;
t29 = t32 ^ 2;
t66 = sin(qJ(4));
t146 = t29 * t66;
t64 = Ifges(7,2) - Ifges(7,1);
t70 = cos(qJ(4));
t167 = (-t70 * t145 + t69 * t146) * t64;
t55 = m(4) + t73;
t79 = pkin(3) ^ 2;
t166 = -t55 * t79 - Ifges(3,2);
t11 = t20 * t69;
t39 = pkin(6) * m(7) + mrSges(6,2) + mrSges(7,3);
t6 = t39 * pkin(5) + t11;
t77 = pkin(5) ^ 2;
t140 = t168 * t77;
t165 = -t140 - Ifges(4,1) - t171;
t160 = pkin(5) * mrSges(5,3);
t78 = pkin(4) ^ 2;
t53 = t73 * t78;
t164 = Ifges(4,2) + t53 - 0.2e1 * t160 + t129 + t165 + t170;
t109 = t74 + t128;
t163 = Ifges(5,2) + t109 + t169;
t71 = cos(qJ(3));
t62 = t71 ^ 2;
t162 = 0.4e1 * t62;
t67 = sin(qJ(3));
t48 = pkin(3) * t67;
t161 = pkin(4) * t67;
t158 = pkin(5) * t67;
t3 = -(2 * t157) - t163;
t61 = t70 ^ 2;
t156 = t3 * t61;
t68 = sin(qJ(2));
t47 = t68 * pkin(2);
t72 = cos(qJ(2));
t153 = (pkin(5) * t71 + t161) * t72 + pkin(2);
t152 = Ifges(7,3) * t66;
t151 = Ifges(7,3) * t69;
t130 = t71 * t72;
t16 = -t67 * t68 + t130;
t150 = t16 * t70;
t19 = t104 * pkin(4) - Ifges(4,4);
t149 = t19 * t67;
t147 = t29 * t64;
t23 = (t60 / 0.2e1 + 0.1e1 / 0.2e1) * t61;
t35 = t47 + pkin(3);
t144 = t35 * t67;
t143 = t39 * t69;
t142 = t39 * t71;
t141 = (t69 + 0.1e1) * (t69 - 0.1e1);
t138 = t67 * t70;
t137 = t68 * t72;
t15 = t67 * t72 + t68 * t71;
t7 = t69 * t15;
t136 = t69 * t67;
t65 = sin(qJ(5));
t135 = t70 * t65;
t134 = t70 * t69;
t133 = t70 * t71;
t132 = t71 * t67;
t131 = t71 * t69;
t126 = pkin(4) * t39 + pkin(7) * t152;
t125 = t62 - 0.1e1 / 0.2e1;
t124 = -0.2e1 * t143;
t123 = -t23 - t60 / 0.2e1 + 0.1e1;
t122 = 0.2e1 * t39;
t121 = 0.2e1 * t71;
t120 = t39 * t161;
t119 = pkin(5) * t143;
t118 = t6 * t135;
t117 = t16 * t134;
t116 = t64 * t146;
t115 = t64 * t145;
t114 = t66 * t145;
t36 = t48 - pkin(5);
t113 = t36 * t143;
t112 = t39 * t135;
t111 = t20 * t136;
t108 = 0.2e1 * t171;
t107 = -t135 / 0.2e1;
t106 = t132 / 0.2e1;
t105 = pkin(3) - t158;
t103 = 0.2e1 * t115;
t102 = mrSges(5,3) + t143;
t101 = -Ifges(4,5) * t71 + Ifges(4,6) * t67;
t100 = pkin(4) * t112;
t99 = t69 * t120;
t98 = (t61 - 0.1e1 / 0.2e1) * t69;
t97 = t69 * t114;
t96 = t131 * t138;
t93 = t56 + t109;
t92 = t125 * t134;
t91 = t65 * t96;
t90 = t94 + t170;
t89 = Ifges(7,3) * (t16 * t135 + t7);
t25 = pkin(4) * t168 + mrSges(4,1);
t88 = (t25 * t71 - t27 * t67) * pkin(3);
t85 = (t29 * t70 + t97) * t64;
t4 = t56 + t163;
t2 = t4 * t61;
t43 = t60 + 0.1e1;
t26 = t43 * t61;
t58 = 0.2e1 * t160;
t83 = 0.2e1 * t70 * t64 * t97 + (t26 - t60) * t147 + Ifges(4,3) + t2 + t53 + t58 - t169;
t13 = -0.2e1 * t100;
t82 = t13 + t2 + t164;
t81 = (t77 + t78) * m(7) + t73 * t77 + t83;
t80 = pkin(2) ^ 2;
t75 = t78 * m(7);
t63 = t72 ^ 2;
t49 = pkin(3) * t71;
t28 = t65 * t152;
t22 = pkin(4) * t142;
t21 = pkin(4) * t71 + t105;
t14 = t68 * (Ifges(4,5) * t67 + Ifges(4,6) * t71);
t12 = t26 + t60 - 0.2e1;
t8 = 0.2e1 * t148;
t1 = t6 * t71 + t120;
t5 = [(((-Ifges(4,4) - t118) * t162 + 0.2e1 * t118 - 0.2e1 * t19 + ((t102 + t159) * t162 + t124) * pkin(4) + 0.2e1 * ((-pkin(4) * t135 - pkin(5) * t69) * t122 + t75 - t156 + t164) * t132) * t68 + 0.2e1 * ((t143 + t27) * t71 + (t25 - t112) * t67) * (t68 * pkin(3) + pkin(2))) * t72 + (t82 + 0.2e1 * t27 * t48 + 0.2e1 * t113 + ((0.2e1 * Ifges(5,1) + (2 * Ifges(6,2)) + 0.2e1 * Ifges(7,3) - t108 - (2 * t74) + t8 - (4 * t157)) * t61 + 0.4e1 * t100 + t8 - 0.2e1 * t75 + 0.2e1 * Ifges(4,1) - (2 * Ifges(6,1)) - 0.2e1 * Ifges(4,2) - (2 * Ifges(7,2)) - 0.2e1 * Ifges(5,3) + 0.4e1 * t102 * pkin(5) + 0.2e1 * t140 - 0.2e1 * t53 + t108) * t62 + ((-0.2e1 * t111 + (pkin(3) - 0.2e1 * t158) * t39) * t135 + 0.2e1 * t99 + 0.2e1 * t149 - t25 * pkin(3)) * t121 + (t78 - t79) * m(7) + Ifges(3,1) + t166) * t63 + (t75 + t82 - 0.2e1 * t119) * t62 + (-pkin(5) + t144) * t124 - 0.2e1 * t27 * t144 - 0.2e1 * ((t12 * t62 + t123 - 0.2e1 * t91) * t63 - 0.2e1 * (t12 * t106 + t65 * t92) * t137 + t123 * t62 + t91 - 0.1e1 / 0.2e1 + t23) * t147 - 0.4e1 * ((-t65 * t132 + t92) * t63 - (t125 * t65 + t96) * t137 + t65 * t106 + (-t62 / 0.2e1 + 0.1e1 / 0.2e1) * t134) * t64 * t114 + t56 + (-(-t111 + (t105 + t47) * t39) * t135 - t99 - t149 + t25 * t35) * t121 + (t76 + t79 + t80) * m(7) + (m(3) + t55) * t80 + 0.2e1 * (mrSges(3,1) + (m(7) + t55) * pkin(3)) * t47 + t156 + t58 + Ifges(2,3) - t165 - t166; -((t43 * t133 - t65 * t136) * t72 - (t65 * t131 + t43 * t138) * t68) * t116 + ((t67 * t107 + t71 * t98) * t72 - t68 * (t65 * t133 + (0.2e1 * t61 - 0.1e1) * t136) / 0.2e1) * t103 + ((t3 * t133 + t65 * (t39 * t105 - t111 + t22)) * t72 - t68 * (t1 * t65 + t3 * t138)) * t66 + (mrSges(4,3) * pkin(3) - Ifges(3,5) + t101) * t72 + t14; t81 + (-(t49 + pkin(4)) * t135 - t36 * t69) * t122 + 0.2e1 * t88 + (m(4) + t168) * t79 + Ifges(3,3); -(t43 * t150 - t65 * t7) * t116 + (t15 * t107 + t16 * t98) * t103 + (-t4 * t150 + ((-t6 * t67 + t22) * t72 - t68 * t1) * t65) * t66 + t101 * t72 + t14; t81 + (-(t49 + 0.2e1 * pkin(4)) * t135 - (t48 - 0.2e1 * pkin(5)) * t69) * t39 + t88; t13 + t140 + t75 + t83 + 0.2e1 * t119; -((t21 * t39 - t111) * t68 + t130 * t11 + t153 * t39) * t135 + t90 * t15 + ((t65 * t117 + t15 * t141) * t29 + t65 * t16 * t114) * t64; ((t36 * t39 - t11) * t66 + t167) * t65; -t65 * (t6 * t66 - t167); t141 * t147 + t90; (-(t21 * t68 + t153) * t143 - t16 * t93) * t66 - pkin(7) * t89 + (-t16 * t146 + (-t15 * t65 + t117) * t145) * t64; (t93 - t113) * t70 - t65 * (pkin(3) * t142 + t126) + t85; (t93 + t119) * t70 - t65 * t126 + t85; -pkin(7) * t151 - t65 * t115; pkin(7) ^ 2 * Ifges(7,3) + t147 + t93; t89; t28; t28; t151; -Ifges(7,3) * pkin(7); Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11), t5(16); t5(2), t5(3), t5(5), t5(8), t5(12), t5(17); t5(4), t5(5), t5(6), t5(9), t5(13), t5(18); t5(7), t5(8), t5(9), t5(10), t5(14), t5(19); t5(11), t5(12), t5(13), t5(14), t5(15), t5(20); t5(16), t5(17), t5(18), t5(19), t5(20), t5(21);];
Mq = res;
