% Calculate kinetic energy for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = CloosQRC350DE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:01:18
% EndTime: 2020-06-23 21:01:29
% DurationCPUTime: 11.15s
% Computational Cost: add. (1620->309), mult. (2233->546), div. (0->0), fcn. (2083->14), ass. (0->193)
t163 = V_base(6) - qJD(1);
t181 = sin(qJ(1));
t186 = cos(qJ(1));
t195 = t181 * V_base(4) + V_base(5) * t186;
t267 = -t195 * pkin(2) + V_base(3);
t182 = cos(qJ(5));
t187 = pkin(6) + rSges(7,3);
t154 = t187 * t182 + pkin(5);
t184 = cos(qJ(3));
t265 = t154 * t184;
t180 = sin(qJ(2));
t255 = t180 * pkin(3);
t264 = t163 * (pkin(2) + t255);
t167 = qJD(2) * t186;
t152 = V_base(5) + t167;
t185 = cos(qJ(2));
t183 = cos(qJ(4));
t177 = sin(qJ(5));
t230 = t177 * t187;
t139 = t183 * t230 - pkin(4);
t179 = sin(qJ(3));
t89 = -t139 * t179 + t265;
t249 = t89 * t185;
t125 = t139 * t184;
t88 = t154 * t179 + t125;
t87 = -pkin(3) + t88;
t263 = t87 * t180 - t249;
t176 = qJ(2) + qJ(3);
t169 = cos(t176);
t161 = t169 ^ 2;
t262 = t161 / 0.2e1;
t261 = pkin(3) - t125;
t124 = t154 * V_base(6);
t166 = qJD(2) * t181;
t153 = V_base(4) + t166;
t168 = sin(t176);
t244 = Icges(4,4) * t169;
t121 = -Icges(4,2) * t168 + t244;
t245 = Icges(4,4) * t168;
t122 = Icges(4,1) * t169 - t245;
t165 = qJD(3) * t186;
t128 = t165 + t152;
t164 = qJD(3) * t181;
t129 = t164 + t153;
t199 = Icges(4,2) * t169 + t245;
t93 = Icges(4,6) * t186 - t199 * t181;
t94 = Icges(4,6) * t181 + t199 * t186;
t200 = Icges(4,1) * t168 + t244;
t95 = Icges(4,5) * t186 - t200 * t181;
t96 = Icges(4,5) * t181 + t200 * t186;
t260 = (t168 * t96 + t169 * t94) * t129 + (t168 * t95 + t169 * t93) * t128 + (t121 * t169 + t122 * t168) * t163;
t158 = t186 * t183;
t178 = sin(qJ(4));
t227 = t181 * t178;
t110 = t168 * t227 - t158;
t259 = t110 ^ 2;
t223 = t186 * t178;
t225 = t181 * t183;
t113 = t168 * t223 + t225;
t258 = t113 ^ 2;
t256 = pkin(3) * t185;
t115 = t168 * t158 - t227;
t222 = t186 * t182;
t207 = -t115 * t177 + t169 * t222;
t112 = -t168 * t225 - t223;
t226 = t181 * t182;
t77 = -t112 * t177 - t169 * t226;
t254 = t77 * t207;
t253 = rSges(3,1) * t180;
t252 = rSges(3,1) * t185;
t224 = t182 * t183;
t106 = -t168 * t177 + t169 * t224;
t162 = pkin(7) * qJ(5) - qJ(6);
t156 = sin(t162);
t157 = cos(t162);
t234 = t169 * t178;
t69 = -t106 * t157 - t156 * t234;
t251 = Icges(7,1) * t69;
t68 = -t106 * t156 + t157 * t234;
t250 = Icges(7,2) * t68;
t248 = Icges(5,1) * t183;
t247 = Icges(6,1) * t106;
t246 = Icges(2,4) * t181;
t243 = Icges(3,5) * t181;
t174 = t185 ^ 2;
t242 = Icges(3,2) * t174;
t241 = Icges(5,2) * t178;
t231 = t177 * t183;
t105 = -t168 * t182 - t169 * t231;
t240 = Icges(6,2) * t105;
t239 = Icges(5,3) * t161;
t238 = Icges(5,3) * t186;
t237 = Icges(7,3) * t105;
t236 = t110 * t113;
t235 = t184 * t124;
t233 = t169 * t181;
t232 = t169 * t186;
t229 = t180 * t181;
t228 = t180 * t186;
t219 = qJD(4) * t169;
t218 = qJD(5) * t187;
t217 = V_base(5) * pkin(1) + V_base(1);
t216 = t186 * t242;
t215 = Icges(6,3) * t234;
t214 = t178 * t230;
t212 = t181 * V_base(6);
t211 = t186 * V_base(6);
t209 = qJD(4) * t230;
t208 = pkin(2) + t253;
t206 = pkin(4) * t168 + pkin(5) * t169;
t205 = rSges(4,1) * t168 + rSges(4,2) * t169;
t202 = -V_base(4) * pkin(1) + V_base(2);
t201 = V_base(6) * t214;
t100 = t186 * t219 + t129;
t198 = Icges(4,5) * t168 + Icges(4,6) * t169;
t136 = -qJD(4) * t168 + t163;
t67 = qJD(5) * t113 + t100;
t133 = -Icges(3,1) * t229 + Icges(3,5) * t186;
t196 = t139 * V_base(6);
t101 = qJD(5) * t234 + t136;
t194 = (Icges(4,5) * t169 - Icges(4,6) * t168) * t163 + (Icges(4,3) * t186 - t198 * t181) * t128 + (Icges(4,3) * t181 + t198 * t186) * t129;
t99 = -t181 * t219 + t128;
t66 = -qJD(5) * t110 + t99;
t193 = t152 * t256 + t264 * t181 + t217;
t107 = t206 * t181;
t127 = t169 * pkin(4) - t168 * pkin(5);
t192 = t163 * t107 + t128 * t127 + t193;
t191 = -t153 * t256 + t264 * t186 + t202;
t108 = t206 * t186;
t190 = t163 * t108 - t129 * t127 + t191;
t189 = (-qJD(2) - t195) * t255 + t267;
t188 = -t129 * t107 - t128 * t108 + t189;
t175 = t186 ^ 2;
t173 = t181 ^ 2;
t172 = t178 ^ 2;
t170 = Icges(2,4) * t186;
t160 = -pkin(7) * qJD(5) + qJD(6);
t148 = t186 * rSges(2,1) + t181 * rSges(2,2);
t147 = t181 * rSges(2,1) - t186 * rSges(2,2);
t146 = t185 * Icges(3,2) * t229;
t145 = Icges(2,1) * t186 + t246;
t144 = -Icges(2,1) * t181 + t170;
t143 = Icges(2,2) * t181 + t170;
t142 = Icges(2,2) * t186 - t246;
t134 = Icges(3,1) * t228 + t243;
t132 = Icges(3,5) * t228 + Icges(3,3) * t181;
t131 = -Icges(3,5) * t229 + Icges(3,3) * t186;
t130 = t185 * t179 + t180 * t184;
t126 = t169 * rSges(4,1) - t168 * rSges(4,2);
t117 = -t179 * t177 + t184 * t224;
t116 = t184 * t177 + t179 * t224;
t109 = t195 * t187;
t102 = t105 ^ 2;
t98 = t181 * rSges(4,3) + t205 * t186;
t97 = t186 * rSges(4,3) - t205 * t181;
t86 = V_base(5) * rSges(2,3) + t163 * t147 + t217;
t85 = t163 * t148 + V_base(2) + (-pkin(1) - rSges(2,3)) * V_base(4);
t83 = -V_base(4) * t147 - V_base(5) * t148 + V_base(3);
t82 = t89 * t180;
t80 = t115 * t182 + t177 * t232;
t78 = t112 * t182 - t177 * t233;
t75 = t207 ^ 2;
t74 = t77 ^ 2;
t73 = t154 * V_base(4) - t186 * t196;
t72 = -t154 * V_base(5) - t181 * t196;
t71 = t116 * t185 + t180 * t117;
t70 = -t195 * pkin(4) + t109 * t231;
t65 = t160 * t105 + t101;
t64 = t152 * t252 + (-t186 * rSges(3,3) + t208 * t181) * t163 + t217;
t63 = -t153 * t252 + (t181 * rSges(3,3) + t208 * t186) * t163 + t202;
t62 = V_base(3) - qJD(2) * t253 + (V_base(4) * rSges(3,3) - t208 * V_base(5)) * t186 + (-V_base(5) * rSges(3,3) - t208 * V_base(4)) * t181;
t61 = t88 * t185 + t82;
t60 = t87 * t185 + t82;
t59 = pkin(2) - t263;
t58 = -t113 * t156 - t80 * t157;
t57 = t113 * t157 - t80 * t156;
t56 = t110 * t156 - t78 * t157;
t55 = -t110 * t157 - t78 * t156;
t52 = t160 * t207 + t67;
t51 = t160 * t77 + t66;
t49 = -t129 * t126 + t163 * t98 + t191;
t48 = t128 * t126 - t163 * t97 + t193;
t47 = -t128 * t98 + t129 * t97 + t189;
t44 = (t100 * t168 + t136 * t232) * rSges(5,3) + t190;
t43 = (t136 * t233 - t168 * t99) * rSges(5,3) + t192;
t42 = (-t100 * t181 - t186 * t99) * t169 * rSges(5,3) + t188;
t40 = (t101 * t207 - t105 * t67) * rSges(6,2) + t190;
t39 = (-t101 * t77 + t105 * t66) * rSges(6,2) + t192;
t38 = (-t207 * t66 + t67 * t77) * rSges(6,2) + t188;
t34 = t263 * qJD(2) + (t88 * t180 - t249) * qJD(3) - (-t116 * t180 + t185 * t117) * t218 + t178 * (-t180 * t179 + t185 * t184) * t209 + (t70 * t184 + (t195 * pkin(5) + t109 * t182) * t179 - t195 * pkin(3)) * t180 + (t70 * t179 - t195 * t265) * t185 + (-V_base(5) * t181 + V_base(4) * t186) * t214 + t267;
t33 = (-t181 * t59 + t186 * t214) * qJD(1) - t60 * t167 - t61 * t165 - (-t178 * t226 + t71 * t186) * t218 + (t130 * t223 + t225) * t209 - t186 * t201 + (t72 * t179 + t181 * t235 + t261 * V_base(5)) * t185 + (t72 * t184 + (-t124 * t181 + t139 * V_base(5)) * t179 + pkin(3) * t212) * t180 + pkin(2) * t212 + t217;
t32 = (-t181 * t214 - t59 * t186) * qJD(1) + t60 * t166 + t61 * t164 + (t178 * t222 + t71 * t181) * t218 + (-t130 * t227 + t158) * t209 + t181 * t201 + (t73 * t179 + t186 * t235 - t261 * V_base(4)) * t185 + (t73 * t184 + (-t124 * t186 - t139 * V_base(4)) * t179 + pkin(3) * t211) * t180 + pkin(2) * t211 + t202;
t1 = m(3) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(2) * (t83 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + (Icges(6,1) * t80 ^ 2 + Icges(6,2) * t75 + Icges(6,3) * t258) * t67 ^ 2 / 0.2e1 + (Icges(6,1) * t78 ^ 2 + Icges(6,2) * t74 + Icges(6,3) * t259) * t66 ^ 2 / 0.2e1 + t100 * t99 * (t115 * Icges(5,1) * t112 - t161 * t181 * t238 - Icges(5,2) * t236) + (Icges(5,1) * t115 ^ 2 + Icges(5,2) * t258 + t175 * t239) * t100 ^ 2 / 0.2e1 + t152 * ((t186 * t132 + (-t134 * t180 - t216) * t181) * t153 + (t186 * t131 - t133 * t229 + t173 * t242) * t152 + t163 * (t133 * t185 + t146)) / 0.2e1 + t153 * ((t181 * t132 + t134 * t228 + t175 * t242) * t153 + (t133 * t228 + (t131 - t216) * t181) * t152 + (t243 + (Icges(3,1) - Icges(3,2)) * t228) * t163 * t185) / 0.2e1 + t129 * (t194 * t181 + t260 * t186) / 0.2e1 + t128 * (-t260 * t181 + t194 * t186) / 0.2e1 + (Icges(7,1) * t69 ^ 2 + Icges(7,2) * t68 ^ 2 + Icges(7,3) * t102) * t65 ^ 2 / 0.2e1 + (Icges(7,1) * t58 ^ 2 + Icges(7,2) * t57 ^ 2 + Icges(7,3) * t75) * t52 ^ 2 / 0.2e1 + ((t168 ^ 2 * Icges(5,3) / 0.2e1 + (Icges(5,1) * t183 ^ 2 + Icges(5,2) * t172) * t262) * t136 + ((t113 * t241 + t115 * t248 - t168 * t238) * t100 + (Icges(5,3) * t168 * t181 - t110 * t241 + t112 * t248) * t99) * t169) * t136 + V_base(5) * ((t186 * t143 - t181 * t145) * V_base(4) + (t186 * t142 - t181 * t144) * V_base(5)) / 0.2e1 + t65 * t52 * (t207 * t237 + t57 * t250 + t58 * t251) + t66 * t67 * (t80 * Icges(6,1) * t78 + Icges(6,2) * t254 - Icges(6,3) * t236) + ((t58 * Icges(7,1) * t56 + t57 * Icges(7,2) * t55 + Icges(7,3) * t254) * t52 + (t77 * t237 + t55 * t250 + t56 * t251) * t65 + (Icges(7,1) * t56 ^ 2 / 0.2e1 + Icges(7,2) * t55 ^ 2 / 0.2e1 + Icges(7,3) * t74 / 0.2e1) * t51) * t51 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(7) * (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + t163 * V_base(5) * (-Icges(2,5) * t181 + Icges(2,6) * t186) + t163 * V_base(4) * (Icges(2,5) * t186 + Icges(2,6) * t181) + (Icges(5,1) * t112 ^ 2 + Icges(5,2) * t259 + t173 * t239) * t99 ^ 2 / 0.2e1 + ((-t168 * t94 + t169 * t96) * t129 + (-t168 * t93 + t169 * t95) * t128 + t146 * t152 + ((-Icges(3,2) * t228 + t134) * t153 + t133 * t152) * t185 + (t180 ^ 2 * Icges(3,2) + Icges(3,1) * t174 - t168 * t121 + t169 * t122 + Icges(2,3)) * t163) * t163 / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + V_base(4) * ((t181 * t143 + t186 * t145) * V_base(4) + (t181 * t142 + t186 * t144) * V_base(5)) / 0.2e1 + ((-t110 * t215 + t77 * t240 + t78 * t247) * t66 + (t113 * t215 + t207 * t240 + t80 * t247) * t67 + (Icges(6,1) * t106 ^ 2 / 0.2e1 + t172 * Icges(6,3) * t262 + Icges(6,2) * t102 / 0.2e1) * t101) * t101 + m(6) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(5) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1;
T = t1;
