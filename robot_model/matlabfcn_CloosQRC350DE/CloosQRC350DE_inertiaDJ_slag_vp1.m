% Calculate time derivative of joint inertia matrix for
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350DE_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:02:09
% EndTime: 2020-06-23 21:02:31
% DurationCPUTime: 14.83s
% Computational Cost: add. (1895->348), mult. (2972->606), div. (0->0), fcn. (1758->14), ass. (0->234)
t306 = -Icges(4,1) + Icges(4,2);
t305 = Icges(7,3) + Icges(6,2);
t303 = m(3) * rSges(3,1);
t103 = sin(qJ(2));
t102 = sin(qJ(3));
t107 = cos(qJ(3));
t106 = cos(qJ(4));
t100 = sin(qJ(5));
t110 = pkin(6) + rSges(7,3);
t235 = t100 * t110;
t60 = t106 * t235 - pkin(4);
t105 = cos(qJ(5));
t68 = t105 * t110 + pkin(5);
t31 = -t102 * t60 + t107 * t68;
t302 = t103 * t31;
t108 = cos(qJ(2));
t134 = t102 * t103 - t107 * t108;
t101 = sin(qJ(4));
t222 = qJD(4) * t101;
t188 = t134 * t222;
t84 = pkin(7) * qJ(5) - qJ(6);
t73 = sin(t84);
t74 = cos(t84);
t272 = t73 * t74;
t99 = Icges(7,2) - Icges(7,1);
t208 = t99 * t272;
t242 = pkin(7) * qJD(5);
t82 = -qJD(6) + t242;
t164 = t82 * t208;
t220 = qJD(4) * t106;
t72 = t74 ^ 2;
t192 = t72 * t220;
t280 = 0.2e1 * t101;
t301 = t164 * t280 - t99 * t192;
t293 = 0.2e1 * t72;
t212 = t99 * t293;
t271 = t82 * t99;
t71 = t73 ^ 2;
t294 = -0.2e1 * t71;
t300 = t82 * t212 + t271 * t294;
t299 = -t102 * t68 - t107 * t60;
t92 = qJD(2) + qJD(3);
t158 = t106 * t92 + qJD(5);
t98 = qJ(2) + qJ(3);
t86 = sin(t98);
t248 = t105 * t86;
t87 = cos(t98);
t298 = t100 * (t101 * (qJD(4) * t86 + qJD(1)) - t158 * t87) - t92 * t248;
t204 = t105 * t272;
t157 = t101 * t204;
t219 = qJD(5) * t100;
t178 = t101 * t219;
t203 = t106 * t272;
t217 = qJD(5) * t105;
t267 = t71 - t72;
t297 = ((t101 * t105 * t72 - t203) * t217 + ((t106 * t267 - 0.2e1 * t157) * t82 + t105 * t192 - t178 * t72 + t222 * t272) * t100) * t99;
t295 = 0.2e1 * pkin(4);
t95 = t105 ^ 2;
t265 = t100 ^ 2 - t95;
t144 = m(6) * rSges(6,2) ^ 2 + m(7) * t110 ^ 2;
t130 = t144 - t305;
t291 = Icges(7,1) + Icges(6,3);
t292 = Icges(6,1) + Icges(7,2);
t54 = t130 + t292;
t50 = t54 * t95;
t21 = Icges(5,1) - Icges(5,2) - t130 - t291 + t50;
t251 = t105 * t54;
t69 = m(6) * rSges(6,2) + m(7) * t110;
t279 = pkin(5) * t69;
t29 = t251 + t279;
t243 = t107 * t69;
t57 = pkin(4) * t243;
t15 = -t102 * t29 + t57;
t278 = rSges(4,3) * m(4);
t79 = rSges(4,2) * t278 - Icges(4,6);
t80 = rSges(4,1) * t278 - Icges(4,5);
t165 = -t102 * t79 + t80 * t107;
t52 = t102 * t108 + t103 * t107;
t23 = t92 * t52;
t131 = -t106 * t23 + t188;
t96 = t106 ^ 2;
t229 = t96 - 0.1e1 / 0.2e1;
t290 = t158 * t86 + t87 * t222;
t236 = t100 * t106;
t240 = qJD(1) * t86;
t286 = (qJD(4) + t240) * t236 - qJD(1) * t87 * t105;
t283 = -0.2e1 * pkin(5);
t282 = -0.2e1 * t82;
t281 = 0.2e1 * t96;
t277 = -t103 / 0.2e1;
t276 = -t106 / 0.2e1;
t275 = pkin(4) * t102;
t274 = pkin(5) * t102;
t273 = t71 * t82;
t270 = t86 * t92;
t269 = t87 * t92;
t90 = t105 + 0.1e1;
t91 = t105 - 0.1e1;
t268 = t90 * t91;
t260 = Icges(7,3) * pkin(7);
t266 = pkin(4) * t69 + t101 * t260;
t264 = rSges(6,2) * t100;
t259 = t100 * t52;
t257 = t101 * t87;
t256 = t101 * t99;
t65 = pkin(5) * t107 + t275;
t252 = t103 * t65;
t250 = t105 * t69;
t249 = t105 * t82;
t246 = t106 * t134;
t245 = t106 * t86;
t239 = qJD(3) * t96;
t238 = qJD(5) * t52;
t237 = t100 * t105;
t234 = t102 * t105;
t233 = t102 * t106;
t232 = t105 * t106;
t231 = t106 * t107;
t230 = t107 * t105;
t104 = sin(qJ(1));
t228 = qJD(1) * t104;
t226 = qJD(2) * t103;
t225 = qJD(2) * t108;
t224 = qJD(3) * t102;
t223 = qJD(3) * t107;
t221 = qJD(4) * t105;
t218 = qJD(5) * t102;
t216 = t106 * qJD(5);
t215 = t107 * qJD(5);
t213 = m(5) + m(6) + m(7);
t211 = pkin(5) * t270;
t210 = pkin(4) * t269;
t209 = t69 * t275;
t207 = t86 * t269;
t206 = pkin(4) * t240;
t205 = t72 * t256;
t201 = pkin(3) * t223;
t200 = pkin(3) * t225;
t199 = t106 * t260;
t197 = t21 * t231;
t196 = t21 * t233;
t195 = t87 * t228;
t190 = t100 * t232;
t189 = t100 * t231;
t75 = pkin(3) * t103 + pkin(2);
t187 = -t86 * pkin(4) - t75;
t186 = t72 * (-t90 - t91);
t184 = t104 * t225;
t109 = cos(qJ(1));
t183 = t109 * t225;
t182 = t106 * t224;
t181 = t100 * t220;
t180 = t101 * t220;
t179 = t107 * t222;
t177 = t100 * t218;
t176 = t100 * t217;
t175 = t101 * t217;
t174 = t105 * t216;
t173 = t100 * t215;
t172 = -t236 / 0.2e1;
t169 = 0.2e1 * t208;
t168 = -0.2e1 * t207;
t167 = rSges(3,1) * t103 + pkin(2);
t166 = t222 / 0.2e1;
t59 = t75 * t228;
t163 = pkin(5) * t195 + t104 * t206 + t109 * t211 + t59;
t162 = t105 * t229;
t161 = t99 * t204;
t160 = -0.2e1 * t189;
t159 = -0.2e1 * t100 * t233;
t156 = -0.2e1 * t180;
t153 = pkin(3) * t182;
t151 = t219 * t272;
t150 = t134 * t190;
t149 = t87 ^ 2 * t180;
t148 = t54 * t176;
t146 = -0.2e1 * t164;
t136 = (-t92 - t216) * t87;
t9 = t100 * t136 - t105 * t290;
t145 = -t257 * t82 - t9;
t24 = t92 * t134;
t143 = -t134 * t216 - t24;
t40 = -0.2e1 * t148;
t141 = pkin(4) * t107 - t274;
t140 = rSges(4,1) * t86 + rSges(4,2) * t87;
t132 = -t140 - t75;
t16 = t107 * t29 + t209;
t63 = t144 + t291;
t129 = (-pkin(5) - rSges(5,3)) * t87 + t187;
t85 = t95 + 0.1e1;
t41 = t100 * t230 + t233 * t85;
t42 = -t100 * t234 + t231 * t85;
t127 = (-t102 * (pkin(4) * t213 + m(4) * rSges(4,1)) - t107 * (pkin(5) * t213 + m(4) * rSges(4,2) + m(5) * rSges(5,3))) * qJD(3) * pkin(3);
t126 = t175 + t181;
t125 = t100 * t222 - t174;
t123 = t109 * t129;
t122 = (-rSges(6,2) * t105 - pkin(5)) * t87 + t187;
t121 = t131 - t238;
t45 = -t100 * t86 + t232 * t87;
t120 = -t101 * t270 + t220 * t87 - t45 * t82;
t119 = -t200 + (-pkin(4) * t87 + rSges(5,3) * t86) * t92;
t118 = -rSges(4,3) * t104 + t109 * t132;
t117 = -t200 + (-rSges(4,1) * t87 + rSges(4,2) * t86) * t92;
t94 = t101 ^ 2;
t115 = t212 * (-t85 * t180 + (-t96 + 0.1e1) * t176) - t21 * t156 + (t85 * t96 - t95) * t146 + t148 * t281 + t40 + (t281 - 0.2e1 * t94) * qJD(4) * t161 + (-0.2e1 * t151 + (t294 + t293) * t249) * t106 * t256;
t113 = -t103 * t165 * qJD(3) + (-qJD(3) * t108 - t225) * (t102 * t80 + t107 * t79);
t112 = ((t282 + t221) * t203 + (-qJD(4) * t72 - t249 * t267 - t151) * t101) * t99;
t89 = pkin(3) * t107;
t88 = pkin(3) * t102;
t78 = t89 + pkin(4);
t77 = t89 + t295;
t76 = t88 - pkin(5);
t64 = pkin(3) + t141;
t55 = t104 * t211;
t49 = t101 * t104 - t109 * t245;
t48 = -t101 * t109 - t104 * t245;
t47 = t126 * Icges(7,3);
t43 = t54 * t177;
t39 = t125 * t110;
t34 = t189 + (t281 - 0.1e1) * t234;
t30 = t102 * t172 + t107 * t162;
t28 = pkin(3) + t299;
t26 = t103 * t64 + t108 * t65 + pkin(2);
t22 = t69 * t76 - t251;
t19 = -t134 * t232 - t259;
t14 = t57 - t54 * t234 + t69 * (pkin(3) - t274);
t8 = t100 * t290 + t105 * t136;
t7 = (t108 * t64 - t252) * qJD(2) + (t108 * t141 - t252) * qJD(3);
t6 = qJD(3) * t15 - t54 * t173;
t5 = -t134 * t251 + t26 * t69;
t4 = t100 * t16 + t196;
t3 = t100 * t14 + t197;
t2 = -t103 * t16 + t108 * t15;
t1 = [0.2e1 * m(6) * ((t104 * t122 - t264 * t48) * (-pkin(3) * t184 - t104 * t210 + t55 + (-t104 * t298 - t48 * t217) * rSges(6,2) + (-t206 + t286 * rSges(6,2) + (-t87 * pkin(5) - t75) * qJD(1)) * t109) + (t109 * t122 - t264 * t49) * (-pkin(3) * t183 - t109 * t210 - (t286 * t104 + t109 * t298 + t49 * t217) * rSges(6,2) + t163)) + 0.2e1 * m(5) * (t129 * t104 * (qJD(1) * t123 + t104 * t119 + t55) + (rSges(5,3) * t195 + t109 * t119 + t163) * t123) + 0.2e1 * m(4) * ((rSges(4,3) * t109 + t104 * t132) * (qJD(1) * t118 + t104 * t117) + t118 * (t59 + t140 * t228 + (-rSges(4,3) * qJD(1) + t117) * t109)) + 0.2e1 * (-(-rSges(3,3) * t104 - t109 * t167) * t183 + t184 * (-rSges(3,3) * t109 + t104 * t167)) * t303 + 0.2e1 * Icges(5,3) * t207 + 0.2e1 * Icges(7,1) * (-t257 * t73 - t45 * t74) * (-t120 * t73 + t145 * t74) + 0.2e1 * Icges(7,2) * (t257 * t74 - t45 * t73) * (t120 * t74 + t145 * t73) + 0.2e1 * Icges(6,1) * t45 * t9 + (0.2e1 * Icges(4,4) * t86 + t306 * t87) * t270 + (-0.2e1 * Icges(4,4) * t87 + t306 * t86) * t269 + 0.2e1 * m(7) * (t104 ^ 2 + t109 ^ 2) * (-(t103 * t28 + t108 * t31 + pkin(2)) * (-(t107 * t39 + t110 * t177) * t103 - (t102 * t39 - t110 * t173) * t108 - (t108 * t299 - t302) * qJD(3) - (t108 * t28 - t302) * qJD(2)) + t101 * t235 * t110 * t126) + 0.2e1 * t305 * (-t236 * t87 - t248) * t8 + (-2 * Icges(3,1) + 2 * Icges(3,2)) * t103 * t225 + (t168 * t96 - 0.2e1 * t149) * Icges(5,1) + (Icges(5,2) + Icges(6,3)) * (t94 * t168 + 0.2e1 * t149); -(t85 * t188 + ((t102 * t265 + t105 * t160) * t108 - t103 * (t105 * t159 - t107 * t265)) * qJD(5) + t92 * (-t103 * t42 - t108 * t41)) * t205 + ((-t30 * t103 - t108 * t34 / 0.2e1) * qJD(2) + ((t108 * t156 + (-qJD(3) + t216 + 0.2e1 * t239) * t277) * t107 + ((-t239 - t216 / 0.2e1 + qJD(3) / 0.2e1) * t108 + 0.2e1 * t103 * t180) * t102) * t105 + ((t102 * t166 - t215 * t229 + t223 * t276) * t108 + (-0.2e1 * t218 * t96 - t179 - t182 + t218) * t277) * t100) * t169 + ((t100 * t43 - t21 * t179) * t108 - t103 * (-t102 * t21 * t222 + t100 * t6) + ((t160 * t54 + t14) * t108 - t103 * (t159 * t54 + t16)) * t217 + ((-t196 + t100 * (-pkin(5) * t243 - t230 * t54 - t209)) * t108 - t103 * t197) * qJD(3) + (-t103 * t3 - t108 * t4) * qJD(2)) * t101 + (-t103 * t4 + t108 * t3) * t220 - (pkin(3) * t278 + rSges(3,3) * t303 - Icges(3,5) + t165) * t226 + t113 + t301 * (-t103 * t41 + t108 * t42) + t300 * (t30 * t108 + t277 * t34); 0.2e1 * t127 + 0.2e1 * ((-t216 * t78 - t201) * t105 + (qJD(5) * t76 + t222 * t78 + t153) * t100) * t69 + t115; -(t24 * t237 + t131 * t85 + (t265 * t52 + 0.2e1 * t150) * qJD(5)) * t205 + ((-t229 * t23 + (0.2e1 * t188 - t238 / 0.2e1) * t106) * t105 + (qJD(5) * t134 * t229 + t166 * t52 - t24 * t276) * t100) * t169 + (-t40 * t246 + t2 * t217 + t100 * ((-qJD(3) * t16 + t43) * t108 - t103 * t6 + (-t103 * t15 - t108 * t16) * qJD(2)) + t131 * t21) * t101 + (t100 * t2 - t21 * t246) * t220 - t165 * t226 + t113 + t301 * (-t237 * t52 - t246 * t85) + t300 * (-t134 * t162 + t172 * t52); t127 + ((-t216 * t77 - t201) * t105 + (t153 + t77 * t222 + (t88 + t283) * qJD(5)) * t100) * t69 + t115; (t125 * t295 + t219 * t283) * t69 + t115; -t5 * t174 + (-Icges(5,3) - t144 + t50 - t292) * t24 + (t5 * t222 - t106 * (-t23 * t251 + t69 * t7) + (0.2e1 * t105 * t52 - t134 * t236) * t54 * qJD(5)) * t100 + ((-t190 * t23 - t24 * t268) * t72 + ((t268 * t52 - t150) * t282 - t134 * t181) * t272 + (-t23 * t272 - (-t273 + (t82 - t221) * t72) * t134) * t100 * t101 + (t186 * t259 - (-t106 * t265 * t72 + t157) * t134) * qJD(5)) * t99; ((t201 * t69 + t219 * t54) * t101 + t22 * t220) * t100 + t22 * t175 + t297; -(-t178 * t54 + t220 * t29) * t100 - t29 * t175 + t297; t146 * t268 + (t186 * t99 + 0.2e1 * t251) * t219; (t23 * t63 + (-t105 * t7 + t219 * t26) * t69) * t101 + (t134 * t63 - t250 * t26) * t220 - (t100 * t121 + t105 * t143) * t260 + (-t19 * t273 + (-t134 * t280 * t82 - t100 * t143 + t105 * t121) * t272 + (t101 * t23 + t134 * t220 + t19 * t82) * t72) * t99; (-t105 * t201 + t219 * t76) * t69 * t106 - (-t250 * t76 + t63) * t222 - (-pkin(3) * t224 * t69 + qJD(4) * t199) * t100 - (pkin(3) * t243 + t266) * t217 + t112; (-t105 * t266 - t236 * t279) * qJD(5) + (-(pkin(5) * t250 + t63) * t101 - t100 * t199) * qJD(4) + t112; -qJD(5) * t161 + (Icges(7,3) * t242 + t267 * t271) * t100; t146; -Icges(7,3) * t8; t47; t47; -Icges(7,3) * t219; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
