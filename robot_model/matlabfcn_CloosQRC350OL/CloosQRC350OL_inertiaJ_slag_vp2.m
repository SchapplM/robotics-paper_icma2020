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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 08:00:14
% EndTime: 2020-06-20 08:00:18
% DurationCPUTime: 3.93s
% Computational Cost: add. (2306->366), mult. (4922->524), div. (0->0), fcn. (5371->10), ass. (0->173)
t166 = sin(qJ(3));
t167 = sin(qJ(2));
t171 = cos(qJ(3));
t172 = cos(qJ(2));
t117 = -t166 * t172 - t167 * t171;
t118 = -t166 * t167 + t171 * t172;
t147 = t167 * pkin(3) + pkin(2);
t71 = -pkin(4) * t117 + pkin(5) * t118 + t147;
t259 = 0.2e1 * t71;
t165 = sin(qJ(4));
t210 = t118 * t165;
t169 = cos(qJ(5));
t164 = sin(qJ(5));
t170 = cos(qJ(4));
t204 = t164 * t170;
t64 = t117 * t169 - t118 * t204;
t201 = t169 * t170;
t65 = t117 * t164 + t118 * t201;
t20 = Ifges(6,4) * t65 + Ifges(6,2) * t64 + Ifges(6,6) * t210;
t163 = sin(qJ(6));
t168 = cos(qJ(6));
t29 = t163 * t65 + t168 * t210;
t30 = t163 * t210 - t168 * t65;
t6 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t64;
t258 = t20 + t6;
t224 = Ifges(6,4) * t169;
t202 = t165 * t169;
t109 = -t163 * t202 + t168 * t170;
t110 = t163 * t170 + t168 * t202;
t206 = t164 * t165;
t53 = Ifges(7,5) * t110 + Ifges(7,6) * t109 + Ifges(7,3) * t206;
t257 = t53 + Ifges(6,6) * t170 + (Ifges(6,2) * t164 - t224) * t165;
t255 = mrSges(6,1) * t170 - mrSges(7,1) * t109 + mrSges(7,2) * t110 + mrSges(6,3) * t202;
t225 = Ifges(6,4) * t164;
t134 = Ifges(6,2) * t169 + t225;
t207 = t163 * t164;
t199 = Ifges(7,6) * t207 + Ifges(7,3) * t169;
t205 = t164 * t168;
t94 = -Ifges(7,5) * t205 + t199;
t256 = t94 + t134;
t213 = mrSges(6,1) * t169 - mrSges(6,2) * t164 + mrSges(5,1);
t216 = t170 * mrSges(5,2);
t254 = t213 * t165 + t216;
t209 = t118 * t170;
t230 = -mrSges(5,1) * t117 - mrSges(6,1) * t64 + mrSges(6,2) * t65 + mrSges(5,3) * t209;
t73 = -mrSges(5,2) * t117 - mrSges(5,3) * t210;
t253 = t230 * t165 + t170 * t73;
t228 = mrSges(6,2) * t169;
t112 = (-mrSges(6,1) * t164 - t228) * t165;
t159 = t165 ^ 2;
t162 = t170 ^ 2;
t197 = t159 + t162;
t252 = 0.2e1 * t197 * mrSges(5,3) - 0.2e1 * t165 * t112;
t85 = -mrSges(7,2) * t206 + mrSges(7,3) * t109;
t251 = 0.2e1 * t85;
t86 = mrSges(7,1) * t206 - mrSges(7,3) * t110;
t250 = 0.2e1 * t86;
t121 = -mrSges(6,2) * t170 + mrSges(6,3) * t206;
t249 = 0.2e1 * t121;
t130 = t170 * mrSges(5,1) - t165 * mrSges(5,2);
t248 = 0.2e1 * t130;
t247 = pkin(6) * m(7);
t54 = Ifges(7,4) * t110 + Ifges(7,2) * t109 + Ifges(7,6) * t206;
t246 = t54 / 0.2e1;
t245 = t64 / 0.2e1;
t244 = t65 / 0.2e1;
t222 = Ifges(7,4) * t168;
t96 = Ifges(7,6) * t169 + (Ifges(7,2) * t163 - t222) * t164;
t243 = t96 / 0.2e1;
t223 = Ifges(7,4) * t163;
t98 = Ifges(7,5) * t169 + (-Ifges(7,1) * t168 + t223) * t164;
t242 = t98 / 0.2e1;
t241 = pkin(6) * t85;
t240 = t109 / 0.2e1;
t239 = t110 / 0.2e1;
t133 = Ifges(7,2) * t168 + t223;
t238 = t133 / 0.2e1;
t136 = Ifges(7,1) * t163 + t222;
t237 = t136 / 0.2e1;
t236 = t163 / 0.2e1;
t235 = t164 / 0.2e1;
t234 = t168 / 0.2e1;
t233 = t169 / 0.2e1;
t232 = t170 / 0.2e1;
t231 = pkin(6) * t169;
t229 = mrSges(6,1) * t210 - mrSges(7,1) * t29 + mrSges(7,2) * t30 - mrSges(6,3) * t65;
t227 = Ifges(5,4) * t165;
t226 = Ifges(5,4) * t170;
t221 = Ifges(5,6) * t117;
t68 = t71 ^ 2;
t220 = t159 * t68;
t219 = t163 * t86;
t218 = t165 * t71;
t145 = pkin(3) * t166 - pkin(5);
t146 = pkin(3) * t171 + pkin(4);
t93 = t145 * t201 + t164 * t146;
t217 = t169 * t93;
t215 = t170 * t71;
t186 = mrSges(7,1) * t168 - mrSges(7,2) * t163;
t212 = -t186 + mrSges(6,1);
t154 = t164 * pkin(4);
t125 = -pkin(5) * t201 + t154;
t208 = t125 * t169;
t200 = Ifges(5,5) * t209 + Ifges(5,3) * t117;
t131 = Ifges(7,5) * t163 + Ifges(7,6) * t168;
t132 = Ifges(6,5) * t164 + Ifges(6,6) * t169;
t198 = t163 ^ 2 + t168 ^ 2;
t196 = 0.2e1 * pkin(6);
t19 = Ifges(6,5) * t65 + Ifges(6,6) * t64 + Ifges(6,3) * t210;
t194 = t71 * t206;
t193 = t206 / 0.2e1;
t192 = -t202 / 0.2e1;
t138 = -Ifges(5,1) * t165 - t226;
t99 = Ifges(6,5) * t170 + (-Ifges(6,1) * t169 + t225) * t165;
t191 = -t169 * t99 - t138;
t190 = t197 * t145;
t189 = 0.2e1 * t255;
t188 = m(7) * t198 * pkin(6) ^ 2;
t23 = pkin(6) * t65 + t215;
t31 = (-pkin(6) * t118 - t169 * t71) * t165;
t10 = t163 * t23 - t168 * t31;
t9 = t163 * t31 + t168 * t23;
t187 = t10 * t168 - t163 * t9;
t185 = -mrSges(7,1) * t163 - mrSges(7,2) * t168;
t184 = -Ifges(5,5) * t165 - Ifges(5,6) * t170;
t107 = (t145 - t231) * t165;
t83 = -pkin(6) * t170 + t93;
t43 = t107 * t168 + t163 * t83;
t44 = t107 * t163 - t168 * t83;
t183 = -t163 * t43 + t168 * t44;
t108 = t154 + (-pkin(5) * t169 - pkin(6)) * t170;
t126 = (-pkin(5) - t231) * t165;
t66 = t108 * t163 + t126 * t168;
t67 = -t108 * t168 + t126 * t163;
t182 = -t163 * t66 + t168 * t67;
t95 = -Ifges(6,5) * t202 + Ifges(6,6) * t206 + Ifges(6,3) * t170;
t181 = (t171 * mrSges(4,1) - t166 * mrSges(4,2)) * pkin(3);
t135 = -Ifges(5,2) * t170 - t227;
t55 = Ifges(7,1) * t110 + Ifges(7,4) * t109 + Ifges(7,5) * t206;
t180 = t109 * t54 + t110 * t55 + Ifges(4,3) + t257 * t206 + (-t135 + t95) * t170;
t179 = t109 * t238 + t110 * t237 + t131 * t193 + t168 * t241 + t54 * t234 + t55 * t236 + t95;
t178 = t165 * t191 + t180;
t137 = Ifges(6,1) * t164 + t224;
t177 = t132 * t232 + t137 * t192 + t99 * t235 + t98 * t239 + t96 * t240 + t184 + t257 * t233 + (t246 + t241) * t207 + (-t55 / 0.2e1 + pkin(6) * t86) * t205 + t256 * t193;
t21 = Ifges(6,1) * t65 + Ifges(6,4) * t64 + Ifges(6,5) * t210;
t45 = t221 + (-Ifges(5,2) * t165 + t226) * t118;
t46 = Ifges(5,5) * t117 + (Ifges(5,1) * t170 - t227) * t118;
t7 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t64;
t8 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t64;
t176 = t138 * t209 / 0.2e1 - t170 * t45 / 0.2e1 - t165 * t46 / 0.2e1 + Ifges(4,5) * t118 + t10 * t85 + t9 * t86 + t30 * t55 / 0.2e1 + t29 * t246 + t8 * t239 + t7 * t240 + t99 * t244 + t21 * t192 + t112 * t215 + t19 * t232 - t121 * t202 * t71 + t257 * t245 + (-t135 / 0.2e1 + t95 / 0.2e1) * t210 + t255 * t194 + t258 * t193 + (Ifges(4,6) + t184 / 0.2e1) * t117;
t174 = pkin(5) ^ 2;
t161 = t169 ^ 2;
t158 = t164 ^ 2;
t155 = t159 * t174;
t144 = t145 ^ 2;
t140 = t159 * t144;
t124 = pkin(4) * t169 + pkin(5) * t204;
t122 = mrSges(7,1) * t169 + mrSges(7,3) * t205;
t120 = -mrSges(7,2) * t169 + mrSges(7,3) * t207;
t119 = t124 ^ 2;
t111 = t185 * t164;
t92 = -t145 * t204 + t146 * t169;
t91 = t92 ^ 2;
t72 = t124 * t92;
t69 = (mrSges(5,1) * t165 + t216) * t118;
t57 = t162 * t68;
t56 = t158 * t220;
t42 = t124 * t194;
t36 = -mrSges(6,2) * t210 + mrSges(6,3) * t64;
t26 = t92 * t194;
t17 = mrSges(7,1) * t64 - mrSges(7,3) * t30;
t16 = -mrSges(7,2) * t64 + mrSges(7,3) * t29;
t1 = [0.2e1 * pkin(2) * (mrSges(3,1) * t167 + mrSges(3,2) * t172) + t172 * (Ifges(3,1) * t172 - Ifges(3,4) * t167) - t167 * (Ifges(3,4) * t172 - Ifges(3,2) * t167) + m(4) * t147 ^ 2 + t65 * t21 + t29 * t7 + t30 * t8 + 0.2e1 * t10 * t16 + 0.2e1 * t9 * t17 + Ifges(2,3) + m(3) * pkin(2) ^ 2 + (0.2e1 * t147 * mrSges(4,2) + Ifges(4,1) * t118) * t118 + (-0.2e1 * t147 * mrSges(4,1) + 0.2e1 * Ifges(4,4) * t118 + Ifges(4,2) * t117 + t200) * t117 + t258 * t64 + (t118 * t46 + t230 * t259) * t170 + m(5) * (t57 + t220) + m(6) * (t161 * t220 + t56 + t57) + m(7) * (t10 ^ 2 + t9 ^ 2 + t56) + ((t19 - t45 - t221) * t118 + (t164 * t229 - t169 * t36 - t73) * t259) * t165; t229 * t92 + (t117 * t166 - t118 * t171) * mrSges(4,3) * pkin(3) + Ifges(3,5) * t172 - Ifges(3,6) * t167 + t146 * t69 + t93 * t36 + t43 * t17 + t44 * t16 + m(6) * (t26 + (t145 * t170 - t217) * t218) + t176 + t253 * t145 + m(7) * (t10 * t44 + t43 * t9 + t26); t146 * t248 + t93 * t249 + t44 * t251 + t43 * t250 + t92 * t189 + 0.2e1 * t181 + m(4) * (t166 ^ 2 + t171 ^ 2) * pkin(3) ^ 2 - t252 * t145 + t178 + m(5) * (t144 * t162 + t146 ^ 2 + t140) + m(6) * (t93 ^ 2 + t140 + t91) + m(7) * (t43 ^ 2 + t44 ^ 2 + t91) + Ifges(3,3); -t253 * pkin(5) + m(7) * (t10 * t67 + t66 * t9 + t42) + t229 * t124 + t125 * t36 + t66 * t17 + t67 * t16 + pkin(4) * t69 + m(6) * (t42 + (-pkin(5) * t170 - t208) * t218) + t176; ((-pkin(5) + t145) * t112 + t191) * t165 + (pkin(5) * t197 - t190) * mrSges(5,3) + m(5) * (pkin(4) * t146 - pkin(5) * t190) + m(6) * (-pkin(5) * t145 * t159 + t125 * t93 + t72) + m(7) * (t43 * t66 + t44 * t67 + t72) + (t66 + t43) * t86 + (t67 + t44) * t85 + (t146 + pkin(4)) * t130 + (t93 + t125) * t121 + t181 + t180 + t255 * (t92 + t124); m(7) * (t66 ^ 2 + t67 ^ 2 + t119) + m(6) * (t125 ^ 2 + t119 + t155) + m(5) * (pkin(4) ^ 2 + t162 * t174 + t155) + t125 * t249 + pkin(4) * t248 + t67 * t251 + t66 * t250 + t124 * t189 + t252 * pkin(5) + t178; t9 * t122 + t137 * t244 + t10 * t120 + t29 * t243 + t30 * t242 + (t134 / 0.2e1 + t94 / 0.2e1) * t64 - t213 * t215 + (t20 / 0.2e1 + t6 / 0.2e1) * t169 + ((t132 / 0.2e1 - Ifges(5,6)) * t118 + (mrSges(5,2) + (-t158 - t161) * mrSges(6,3)) * t71) * t165 + (t21 / 0.2e1 - t168 * t8 / 0.2e1 + t7 * t236 + t111 * t218 + (m(7) * (t10 * t163 + t168 * t9) + t163 * t16 + t168 * t17) * pkin(6)) * t164 + t200; -t254 * t145 + t43 * t122 + t44 * t120 + t92 * t111 + mrSges(6,3) * t217 + (-t92 * mrSges(6,3) + (t163 * t44 + t168 * t43) * t247) * t164 + t177; t254 * pkin(5) + t66 * t122 + t124 * t111 + t67 * t120 + mrSges(6,3) * t208 + (-t124 * mrSges(6,3) + (t163 * t67 + t168 * t66) * t247) * t164 + t177; Ifges(5,3) + t256 * t169 + t158 * t188 + (t163 * t96 - t168 * t98 + t137 + (t120 * t163 + t122 * t168) * t196) * t164; t29 * t238 + t131 * t245 + t8 * t236 + t7 * t234 + t30 * t237 + t187 * mrSges(7,3) + (t212 * t164 + t228) * t218 + (m(7) * t187 + t168 * t16 - t163 * t17) * pkin(6) + t19; -t93 * mrSges(6,2) + t212 * t92 + t183 * mrSges(7,3) + (m(7) * t183 - t219) * pkin(6) + t179; -t125 * mrSges(6,2) + t212 * t124 + t182 * mrSges(7,3) + (m(7) * t182 - t219) * pkin(6) + t179; t131 * t233 + (pkin(6) * t120 + t243 - t164 * t136 / 0.2e1) * t168 + (-pkin(6) * t122 + t133 * t235 + t242) * t163 + t132; mrSges(7,3) * t196 * t198 + t168 * t133 + t163 * t136 + Ifges(6,3) + t188; mrSges(7,1) * t9 - mrSges(7,2) * t10 + t6; mrSges(7,1) * t43 - mrSges(7,2) * t44 + t53; mrSges(7,1) * t66 - mrSges(7,2) * t67 + t53; (pkin(6) * t186 - Ifges(7,5) * t168) * t164 + t199; pkin(6) * t185 + t131; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
