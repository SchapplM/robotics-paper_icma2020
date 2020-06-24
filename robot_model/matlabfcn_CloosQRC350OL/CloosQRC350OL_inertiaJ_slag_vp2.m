% Calculate joint inertia matrix for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350OL_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:55:58
% EndTime: 2020-06-23 21:56:00
% DurationCPUTime: 2.02s
% Computational Cost: add. (1186->209), mult. (2598->321), div. (0->0), fcn. (2618->10), ass. (0->109)
t120 = Ifges(7,3) + Ifges(6,2);
t80 = sin(qJ(5));
t86 = cos(qJ(4));
t125 = t80 * t86;
t82 = sin(qJ(3));
t83 = sin(qJ(2));
t87 = cos(qJ(3));
t88 = cos(qJ(2));
t47 = -t82 * t88 - t83 * t87;
t48 = -t82 * t83 + t87 * t88;
t85 = cos(qJ(5));
t19 = -t48 * t125 + t47 * t85;
t152 = t120 * t19;
t151 = pkin(6) * m(7);
t150 = Ifges(5,2) + Ifges(6,3);
t122 = t85 * t86;
t20 = t48 * t122 + t47 * t80;
t149 = Ifges(6,1) * t20;
t121 = t86 * mrSges(6,2);
t65 = t83 * pkin(3) + pkin(2);
t24 = -pkin(4) * t47 + pkin(5) * t48 + t65;
t148 = t24 * t121;
t147 = mrSges(7,3) + t151;
t81 = sin(qJ(4));
t75 = t81 ^ 2;
t78 = t86 ^ 2;
t118 = t75 + t78;
t128 = t75 * t85;
t93 = 0.2e1 * mrSges(6,2) * t128 + 0.2e1 * t118 * mrSges(5,3);
t145 = pkin(3) * t82;
t144 = pkin(3) * t87;
t143 = pkin(6) * mrSges(7,3);
t142 = pkin(6) * t85;
t135 = t24 * t85;
t10 = (-pkin(6) * t48 - t135) * t81;
t4 = pkin(6) * t20 + t24 * t86;
t79 = sin(qJ(6));
t84 = cos(qJ(6));
t2 = t10 * t79 + t4 * t84;
t124 = t81 * t85;
t43 = t84 * t124 + t79 * t86;
t141 = t2 * t43;
t132 = t48 * t81;
t8 = t84 * t132 + t20 * t79;
t140 = t79 * t8;
t9 = t79 * t132 - t20 * t84;
t139 = t79 * t9;
t138 = Ifges(6,1) * t85;
t137 = Ifges(7,1) * t84;
t136 = Ifges(7,2) * t84;
t63 = -pkin(5) + t145;
t64 = pkin(4) + t144;
t31 = t63 * t122 + t80 * t64;
t134 = t31 * mrSges(6,2);
t131 = t48 * t86;
t74 = t80 ^ 2;
t130 = t74 * t75;
t77 = t85 ^ 2;
t129 = t75 * t77;
t127 = t79 * t80;
t126 = t80 * t81;
t123 = t84 * mrSges(7,3);
t117 = -0.2e1 * t121;
t115 = 0.2e1 * mrSges(7,3);
t113 = t24 * t126;
t108 = t43 * t137;
t107 = t43 * t123;
t106 = t118 * t63;
t42 = -t79 * t124 + t84 * t86;
t105 = (Ifges(7,2) + t143) * t42 * t127 + t120 * t80 * t124;
t104 = Ifges(6,3) * t86 + (pkin(6) * t123 + t136) * t42 + (Ifges(7,1) + t143) * t43 * t79;
t3 = -t10 * t84 + t4 * t79;
t103 = t2 * t84 + t3 * t79;
t102 = -t2 * t79 + t3 * t84;
t26 = -pkin(6) * t86 + t31;
t40 = (t63 - t142) * t81;
t12 = t26 * t79 + t40 * t84;
t13 = -t26 * t84 + t40 * t79;
t101 = t12 * t84 + t13 * t79;
t70 = t80 * pkin(4);
t41 = t70 + (-pkin(5) * t85 - pkin(6)) * t86;
t52 = (-pkin(5) - t142) * t81;
t21 = t41 * t79 + t52 * t84;
t22 = -t41 * t84 + t52 * t79;
t99 = t21 * t84 + t22 * t79;
t97 = (mrSges(4,1) * t87 - mrSges(4,2) * t82) * pkin(3);
t96 = (mrSges(6,2) * t63 - t138) * t81;
t95 = t43 * Ifges(7,1) * t9 + Ifges(4,5) * t48 + Ifges(4,6) * t47 + (t3 * mrSges(7,3) + Ifges(7,2) * t8) * t42 + t126 * t152 + t150 * t81 * t131;
t94 = t43 ^ 2 * Ifges(7,1) + t42 ^ 2 * Ifges(7,2) + t75 * Ifges(5,1) + Ifges(6,1) * t129 + t120 * t130 + t150 * t78 + Ifges(4,3);
t73 = t79 ^ 2;
t76 = t84 ^ 2;
t92 = (t115 + t151) * pkin(6) * (t73 + t76);
t90 = pkin(5) ^ 2;
t71 = t75 * t90;
t62 = t63 ^ 2;
t61 = Ifges(7,3) * t126;
t57 = t75 * t62;
t51 = -pkin(5) * t122 + t70;
t50 = pkin(4) * t85 + pkin(5) * t125;
t49 = t50 ^ 2;
t30 = -t63 * t125 + t64 * t85;
t29 = t30 ^ 2;
t25 = t50 * t30;
t23 = t24 ^ 2;
t17 = t78 * t23;
t16 = t23 * t130;
t11 = t50 * t113;
t7 = t30 * t113;
t1 = [m(4) * t65 ^ 2 + t88 ^ 2 * Ifges(3,1) + t9 ^ 2 * Ifges(7,1) + t83 ^ 2 * Ifges(3,2) + t8 ^ 2 * Ifges(7,2) + Ifges(2,3) + t120 * t19 ^ 2 + (-t2 * t9 + t3 * t8) * t115 + (m(3) * pkin(2) + 0.2e1 * t83 * mrSges(3,1)) * pkin(2) + (-0.2e1 * t65 * mrSges(4,1) + (Ifges(5,3) + Ifges(4,2)) * t47) * t47 + m(7) * (t2 ^ 2 + t3 ^ 2 + t16) + m(5) * (t23 * t75 + t17) + m(6) * (t23 * t129 + t16 + t17) + (0.2e1 * t148 + t149) * t20 + (0.2e1 * t65 * mrSges(4,2) + 0.2e1 * Ifges(4,4) * t47 + t93 * t24 + (t78 * Ifges(5,1) + t150 * t75 + Ifges(4,1)) * t48) * t48; (-mrSges(4,3) * t144 + (-Ifges(5,1) * t86 - t134) * t81) * t48 + m(7) * (t12 * t2 + t13 * t3 + t7) + t20 * t96 + (-t12 * t9 + t13 * t8 - t141) * mrSges(7,3) + m(6) * (t7 + (-t31 * t85 + t63 * t86) * t81 * t24) + t95 + Ifges(3,5) * t88 + t47 * mrSges(4,3) * t145; m(6) * (t31 ^ 2 + t29 + t57) + m(7) * (t12 ^ 2 + t13 ^ 2 + t29) + m(5) * (t62 * t78 + t64 ^ 2 + t57) - t93 * t63 + t94 + m(4) * (t82 ^ 2 + t87 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t97 + (-t12 * t43 + t13 * t42) * t115 + Ifges(3,3) + t31 * t117; m(7) * (t2 * t21 + t22 * t3 + t11) + m(6) * t11 + (-Ifges(5,1) * t131 - t20 * t138 + m(6) * (-pkin(5) * t86 - t51 * t85) * t24 + (-pkin(5) * t20 - t48 * t51) * mrSges(6,2)) * t81 + (-t21 * t9 + t22 * t8 - t141) * mrSges(7,3) + t95; ((-t31 - t51) * t86 + (pkin(5) - t63) * t128) * mrSges(6,2) + (t118 * pkin(5) - t106) * mrSges(5,3) + ((-t12 - t21) * t43 + (t13 + t22) * t42) * mrSges(7,3) + m(7) * (t12 * t21 + t13 * t22 + t25) + m(6) * (-pkin(5) * t63 * t75 + t31 * t51 + t25) + t97 + m(5) * (pkin(4) * t64 - pkin(5) * t106) + t94; m(7) * (t21 ^ 2 + t22 ^ 2 + t49) + m(6) * (t51 ^ 2 + t49 + t71) + m(5) * (pkin(4) ^ 2 + t78 * t90 + t71) + t93 * pkin(5) + t94 + (-t21 * t43 + t22 * t42) * t115 + t51 * t117; Ifges(5,3) * t47 + t85 * t152 + (t148 - t9 * t137 + Ifges(7,2) * t140 + t149 + t103 * mrSges(7,3) + (m(7) * t103 + (-t84 * t9 + t140) * mrSges(7,3)) * pkin(6)) * t80; (-t108 + t96 + t101 * mrSges(7,3) + (m(7) * t101 - t107) * pkin(6)) * t80 + t105; (-t108 + (-pkin(5) * mrSges(6,2) - t138) * t81 + t99 * mrSges(7,3) + (m(7) * t99 - t107) * pkin(6)) * t80 + t105; Ifges(5,3) + t120 * t77 + (t76 * Ifges(7,1) + t73 * Ifges(7,2) + Ifges(6,1) + t92) * t74; Ifges(7,1) * t139 + t8 * t136 + (mrSges(6,2) * t135 + Ifges(6,3) * t48) * t81 + t102 * mrSges(7,3) + (m(7) * t102 + (t8 * t84 + t139) * mrSges(7,3)) * pkin(6); t104 - t134 + t147 * (-t12 * t79 + t13 * t84); -t51 * mrSges(6,2) + t104 + t147 * (-t21 * t79 + t22 * t84); (-Ifges(7,1) + Ifges(7,2)) * t84 * t127; t73 * Ifges(7,1) + t76 * Ifges(7,2) + Ifges(6,3) + t92; Ifges(7,3) * t19; t61; t61; Ifges(7,3) * t85; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
