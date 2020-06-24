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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350DE_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:02:09
% EndTime: 2020-06-23 21:02:34
% DurationCPUTime: 19.54s
% Computational Cost: add. (2397->452), mult. (3970->698), div. (0->0), fcn. (1974->10), ass. (0->296)
t132 = sin(qJ(5));
t106 = pkin(6) * m(7) + mrSges(6,2) + mrSges(7,3);
t95 = pkin(4) * t106;
t267 = t132 * t95;
t133 = sin(qJ(4));
t287 = qJD(4) * t133;
t137 = cos(qJ(4));
t331 = (m(7) * pkin(6) ^ 2);
t343 = (mrSges(7,3) * pkin(6));
t177 = t331 + 2 * t343;
t160 = -Ifges(6,2) - Ifges(7,3) + t177;
t326 = (Ifges(7,1) + Ifges(6,3));
t272 = Ifges(5,2) + t326;
t136 = cos(qJ(5));
t127 = t136 ^ 2;
t327 = (Ifges(6,1) + Ifges(7,2));
t81 = t160 + t327;
t70 = t81 * t127;
t38 = -Ifges(5,1) + t160 + t272 - t70;
t319 = t137 * t38;
t370 = (t267 - t319) * t287;
t140 = m(5) + m(6);
t358 = m(7) + t140;
t333 = pkin(5) * t358;
t220 = mrSges(5,3) + t333;
t91 = (mrSges(4,2) + t220);
t134 = sin(qJ(3));
t334 = pkin(5) * t106;
t71 = t81 * t136;
t50 = t71 + t334;
t138 = cos(qJ(3));
t314 = t106 * t138;
t84 = pkin(4) * t314;
t20 = -t134 * t50 + t84;
t109 = t358 * pkin(5) ^ 2;
t144 = pkin(4) ^ 2;
t119 = t140 * t144;
t128 = t137 ^ 2;
t349 = -2 * pkin(5);
t367 = Ifges(5,3) + t327 - t70;
t368 = (mrSges(5,3) * t349) + t38 * t128 - Ifges(4,1) + Ifges(4,2) - t109 + t119 - t272 + t367;
t335 = pkin(4) * t134;
t266 = t106 * t335;
t21 = t138 * t50 + t266;
t286 = qJD(4) * t137;
t232 = t133 * t286;
t198 = -0.2e1 * t232;
t283 = qJD(5) * t136;
t92 = t132 * t283;
t83 = t128 * t92;
t325 = t38 * t198 + 0.2e1 * t81 * t83;
t108 = pkin(7) * qJ(5) - qJ(6);
t96 = sin(t108);
t93 = t96 ^ 2;
t97 = cos(t108);
t94 = t97 ^ 2;
t324 = t93 - t94;
t139 = cos(qJ(2));
t135 = sin(qJ(2));
t305 = t135 * t134;
t165 = -t138 * t139 + t305;
t241 = t165 * t287;
t124 = t132 ^ 2;
t297 = t124 - t127;
t315 = t106 * t136;
t211 = mrSges(5,3) + t315;
t303 = t137 * t132;
t13 = -t50 * t303 - Ifges(4,4) + (t211 + t333) * pkin(4);
t179 = -Ifges(4,4) + (t220 + t315) * pkin(4);
t129 = t138 ^ 2;
t345 = 0.2e1 * t129;
t366 = t13 * t345 - t179;
t107 = pkin(7) * qJD(5) - qJD(6);
t131 = Ifges(7,2) - Ifges(7,1);
t329 = t96 * t97;
t262 = t131 * t329;
t203 = t107 * t262;
t164 = 0.2e1 * t133 * t203;
t249 = t94 * t286;
t188 = t131 * t249;
t365 = -t188 + t164;
t313 = t107 * t131;
t256 = t94 * t313;
t257 = t93 * t313;
t364 = 0.2e1 * t256 - 0.2e1 * t257;
t130 = t139 ^ 2;
t363 = (t135 ^ 2 - t130) * qJD(2);
t115 = pkin(3) * t134;
t102 = t115 - pkin(5);
t36 = t102 * t106 - t71;
t285 = qJD(5) * t132;
t229 = t133 * t285;
t258 = qJD(4) * t329;
t261 = t136 * t329;
t312 = t107 * t132;
t362 = (t312 * (-0.2e1 * t133 * t261 + t137 * t324) + (t133 * t258 + t136 * t249 - t229 * t94) * t132 + (t133 * t136 * t94 - t137 * t329) * t283) * t131;
t308 = t134 * t136;
t253 = t81 * t308;
t332 = pkin(5) * t134;
t278 = -0.2e1 * t332;
t32 = -0.2e1 * t253 + (pkin(3) + t278) * t106;
t360 = 0.2e1 * t32;
t359 = 2 * t91;
t344 = 0.2e1 * t135;
t282 = t137 * qJD(5);
t225 = t136 * t282;
t196 = -0.2e1 * t225;
t233 = t132 * t287;
t199 = 0.2e1 * t233;
t357 = (t196 + t199) * t95;
t98 = t127 / 0.2e1 + 0.1e1 / 0.2e1;
t355 = 0.2e1 * t98 * t232 + t83;
t210 = -Ifges(4,5) * t138 + Ifges(4,6) * t134;
t354 = t233 - t225;
t279 = t129 - 0.1e1 / 0.2e1;
t353 = qJD(5) * t279;
t352 = qJD(3) + qJD(2);
t88 = pkin(4) * t358 + mrSges(4,1);
t348 = 0.2e1 * t88;
t346 = 0.2e1 * t128;
t342 = t88 * pkin(3);
t125 = t134 ^ 2;
t341 = -t125 / 0.2e1;
t340 = t129 / 0.2e1;
t339 = -t135 / 0.2e1;
t338 = -t137 / 0.2e1;
t337 = t128 + 0.1e1;
t87 = pkin(5) * t138 + t335;
t336 = t87 * t139 + pkin(2);
t330 = Ifges(7,3) * pkin(7);
t100 = t135 * pkin(3) + pkin(2);
t322 = t131 * t94;
t320 = t132 * t81;
t318 = t137 * t165;
t85 = t98 * t128;
t317 = t133 * t330 + t95;
t316 = t100 * t134;
t117 = t136 + 0.1e1;
t118 = t136 - 0.1e1;
t311 = t117 * t118;
t310 = t131 * t133;
t309 = t132 * t136;
t307 = t134 * t137;
t306 = t134 * t138;
t304 = t135 * t139;
t302 = t137 * t136;
t301 = t137 * t138;
t300 = t138 * t136;
t110 = t127 + 0.1e1;
t89 = t110 * t128;
t72 = t89 + t127 - 0.2e1;
t299 = t72 * qJD(3);
t75 = t134 * t139 + t135 * t138;
t298 = t75 * qJD(5);
t296 = t125 - t129;
t293 = qJD(2) * t135;
t292 = qJD(2) * t139;
t291 = qJD(3) * t128;
t290 = qJD(3) * t134;
t289 = qJD(3) * t137;
t288 = qJD(3) * t138;
t284 = qJD(5) * t134;
t281 = t138 * qJD(5);
t280 = t128 - 0.1e1 / 0.2e1;
t64 = -t85 - t127 / 0.2e1 + 0.1e1;
t276 = 0.2e1 * t132;
t275 = -0.2e1 * t304;
t274 = -0.2e1 * t303;
t271 = 0.2e1 * t272;
t157 = t139 * (Ifges(4,5) * t134 + Ifges(4,6) * t138);
t270 = qJD(2) * t157 + (-t135 * t210 + t157) * qJD(3);
t269 = t92 + t355;
t268 = 0.2e1 * qJD(5);
t264 = t137 * t95;
t263 = pkin(5) * t315;
t260 = pkin(2) * t292;
t259 = pkin(3) * t288;
t255 = t94 * t310;
t254 = t38 * t307;
t252 = t81 * t300;
t251 = t38 * t301;
t250 = Ifges(7,3) * t285;
t248 = t102 * t315;
t247 = t106 * t303;
t246 = t107 * t310;
t245 = t137 * t300;
t244 = t132 * t302;
t243 = t132 * t306;
t242 = t132 * t301;
t240 = pkin(3) - t332;
t239 = t134 * t292;
t237 = t106 * t288;
t236 = t134 * t289;
t235 = t135 * t288;
t234 = t136 * t288;
t231 = t138 * t287;
t230 = t106 * t285;
t228 = t132 * t284;
t227 = t135 * t285;
t226 = t133 * t283;
t224 = t132 * t281;
t222 = t306 / 0.2e1;
t221 = -t303 / 0.2e1;
t219 = t287 / 0.2e1;
t218 = -t282 / 0.2e1;
t217 = t281 / 0.2e1;
t216 = 0.2e1 * t262;
t215 = 0.2e1 * t260;
t178 = -0.2e1 * t203;
t187 = t136 * t246;
t204 = t133 * t262;
t171 = qJD(5) * t204;
t48 = t132 * t171;
t212 = -qJD(4) * t255 + t137 * t178 - t93 * t187 - t48;
t209 = t132 * (-t117 - t118);
t208 = t134 * t274;
t207 = -0.2e1 * t242;
t206 = t136 * t280;
t205 = pkin(4) * t247;
t201 = qJD(2) * t275;
t200 = -0.2e1 * t134 * t288;
t197 = -0.2e1 * t230;
t194 = pkin(3) * t236;
t193 = pkin(3) * t237;
t192 = qJD(5) * t266;
t191 = t131 * t258;
t190 = t165 * t244;
t189 = t94 * t246;
t59 = t81 * t228;
t186 = t134 * t245;
t185 = t81 * t92;
t183 = qJD(3) * t244;
t182 = t134 * t234;
t180 = t134 * t231;
t73 = -0.2e1 * t205;
t56 = -0.2e1 * t185;
t57 = 0.2e1 * t185;
t176 = pkin(4) * t138 - t332;
t175 = t279 * t302;
t174 = t279 * t287;
t172 = t137 * t191;
t167 = t132 * t186;
t163 = -t110 * t232 - t83;
t141 = t144 * m(7);
t162 = t141 + t368;
t28 = t352 * t75;
t161 = -t137 * t28 + t241;
t90 = t177 + t326;
t159 = 0.2e1 * t179;
t158 = 0.4e1 * t179;
t54 = t110 * t301 - t132 * t308;
t55 = t110 * t307 + t132 * t300;
t155 = (-t88 * t134 - t138 * t91) * qJD(3) * pkin(3);
t154 = t161 - t298;
t153 = (t131 * t261 - t132 * t330) * t137;
t150 = 0.2e1 * (t163 + t92) * t322 + (t89 - t127) * t178 + t325 + t56 + (-0.2e1 * t133 ^ 2 + t346) * t136 * t191 + (-0.2e1 * t187 * t324 - 0.2e1 * t48) * t137;
t149 = t368 + t73;
t147 = 0.2e1 * (-pkin(4) * t303 - pkin(5) * t136) * t106 + t162;
t27 = t352 * t165;
t146 = Ifges(7,3) * (t132 * t154 - t136 * t27 - t165 * t225);
t145 = pkin(3) ^ 2;
t121 = m(4) + t140;
t116 = pkin(3) * t138;
t114 = t135 * pkin(2);
t104 = t116 + pkin(4);
t103 = t116 + 0.2e1 * pkin(4);
t101 = t114 + pkin(3);
t82 = pkin(3) + t176;
t79 = pkin(5) * t197;
t78 = t87 * qJD(3);
t74 = t135 * t278 + t100;
t63 = t176 * qJD(3) * t139;
t62 = 0.2e1 * t70;
t58 = (t132 * t286 + t226) * Ifges(7,3);
t49 = t279 * t132 + t186;
t42 = t242 + (t346 - 0.1e1) * t308;
t41 = t175 - t243;
t39 = t134 * t221 + t138 * t206;
t34 = t135 * t82 + t336;
t33 = 0.2e1 * Ifges(5,1) + (2 * Ifges(6,2)) + (2 * Ifges(7,3)) - t271 - (2 * t331) + t62 - (4 * t343);
t30 = t163 - t92;
t26 = t106 * t82 - t253;
t25 = -t253 + (t114 + t240) * t106;
t18 = -t135 * t71 + t106 * (-pkin(5) * t135 + t316);
t17 = t106 * t240 - t253 + t84;
t15 = t132 * t175 + t222 * t72;
t14 = t72 * t129 - 0.2e1 * t167 + t64;
t11 = qJD(3) * t20 - t81 * t224;
t8 = t132 * t21 - t254;
t7 = t134 * t158 + t303 * t360 - 0.2e1 * t342;
t6 = t132 * t17 - t251;
t4 = -t135 * t21 + t139 * t20;
t3 = t41 * t130 - t49 * t304 + t132 * t222 + (-t129 / 0.2e1 + 0.1e1 / 0.2e1) * t302;
t2 = 0.4e1 * pkin(5) * t211 + t128 * t33 + 0.2e1 * Ifges(4,1) - (2 * Ifges(6,1)) - 0.2e1 * Ifges(4,2) - (2 * Ifges(7,2)) - 0.2e1 * Ifges(5,3) + 0.2e1 * t109 - 0.2e1 * t119 - 0.2e1 * t141 + 0.4e1 * t205 + t271 + t62;
t1 = t74 * t315 + t100 * t91 + (t162 + t73) * t305;
t5 = [-0.2e1 * ((t30 * t345 + 0.2e1 * t296 * t183 + (t137 * t268 * t297 + t136 * t199 - 0.2e1 * t299) * t306 + t269) * t130 + t14 * t201 + (t30 * t306 - t174 * t309 + (t340 + t341) * t299 + (-0.2e1 * t132 * t182 - t297 * t353) * t137) * t275 + t64 * t200 - t125 * t183 - t180 * t309 - t355 - t297 * t281 * t307 + 0.2e1 * t363 * t15 + (t269 + t183) * t129) * t322 + 0.4e1 * (t14 * t130 + t15 * t275 + t64 * t129 + t167 - 0.1e1 / 0.2e1 + t85) * t203 - 0.4e1 * (t41 * t201 - (-0.2e1 * qJD(3) - t282) * t243 * t304 + t137 * t182 + t217 * t308 + t363 * t49 + (((-qJD(5) - 0.2e1 * t289) * t306 - t174) * t130 - (-t289 * t296 - t180 + t353) * t304 + t129 * t219 - t287 / 0.2e1) * t136 + ((t282 + qJD(3)) * t340 + (qJD(3) * t296 - t279 * t282) * t130 + qJD(3) * t341 + t218) * t132) * t204 + 0.4e1 * t3 * t93 * t246 + ((t79 + (-t137 * t33 - 0.2e1 * t267) * t287 + (-t320 * t337 + t264) * t136 * t268) * t345 + t2 * t200 + (t225 * t360 + t158 * t288 + (0.4e1 * (-t288 * t50 + t59) * t137 - 0.2e1 * t32 * t287 - 0.4e1 * t192) * t132) * t138 - t7 * t290 + t57 + 0.2e1 * t136 * t193 + t102 * t197 + t259 * t359 + t325 + t357) * t130 + (t149 + 0.2e1 * t248 + t115 * t359 + (t144 - t145) * m(7) - t121 * t145 + t7 * t138 + t2 * t129 + Ifges(3,1) - Ifges(3,2)) * t201 + 0.2e1 * ((t50 * t233 + (-t267 + (t124 * t81 - t136 * t50) * t137) * qJD(5)) * t129 * t344 - 0.4e1 * t13 * t134 * t235 + (-t74 * t230 + t147 * t235 + (0.2e1 * t370 + (-0.2e1 * t264 + (t346 + 0.2e1) * t320) * t283) * t305) * t138 - t1 * t290 - (t100 * t237 + t81 * t227) * t303 + t227 * t95 + t239 * t342 + t88 * t100 * t288 + t354 * t18 + (((t91 + t315) * pkin(3) + t147 * t134) * t138 - t36 * t303 + t366) * t292) * t139 - (0.2e1 * t1 * t138 + t18 * t274 + t316 * t348 + t344 * t366) * t293 + (t370 + (-t136 * t264 + (t337 * t71 + t334) * t132) * qJD(5)) * t345 + (t141 + t149 - 0.2e1 * t263) * t200 + (t88 * t215 + t25 * t196 - t159 * t288 + (-(-t81 * t234 + t59 + (-pkin(5) * t288 + t260) * t106) * t137 + t25 * t287 + t192) * t276) * t138 - (t101 * t348 - t134 * t159 + t25 * t274) * t290 + 0.2e1 * (t101 * t134 - pkin(5)) * t230 + (mrSges(3,1) + (m(7) + t121) * pkin(3)) * t215 - 0.4e1 * (t172 + t189) * t3 - t325 + (-0.2e1 * t315 - (2 * t91)) * (pkin(2) * t239 + t101 * t288); -(t110 * t241 + ((t134 * t297 + t136 * t207) * t139 - (t136 * t208 - t138 * t297) * t135) * qJD(5) + t352 * (-t135 * t54 - t139 * t55)) * t255 + ((-t39 * t135 - t139 * t42 / 0.2e1) * qJD(2) + ((t139 * t198 + (-qJD(3) + t282 + 0.2e1 * t291) * t339) * t138 + ((-t291 + t218 + qJD(3) / 0.2e1) * t139 + t232 * t344) * t134) * t136 + ((-t128 * t281 + t134 * t219 + t288 * t338 + t217) * t139 + (-0.2e1 * t128 * t284 - t231 - t236 + t284) * t339) * t132) * t216 + ((t132 * t59 + t231 * t38) * t139 - t135 * (t134 * t287 * t38 + t11 * t132) + ((t207 * t81 + t17) * t139 - t135 * (t208 * t81 + t21)) * t283 + ((t254 + t132 * (-pkin(5) * t314 - t252 - t266)) * t139 + t135 * t251) * qJD(3) + (-t135 * t6 - t139 * t8) * qJD(2)) * t133 + (-t135 * t8 + t139 * t6) * t286 - (mrSges(4,3) * pkin(3) - Ifges(3,5) + t210) * t293 + t270 + t364 * (t39 * t139 + t339 * t42) + t365 * (-t135 * t55 + t139 * t54); (qJD(5) * t102 + t104 * t287 + t194) * t276 * t106 + t150 + 0.2e1 * t155 + 0.2e1 * (-t104 * t282 - t259) * t315; -(t27 * t309 + t161 * t110 + (t297 * t75 + 0.2e1 * t190) * qJD(5)) * t255 + ((-t280 * t28 + (0.2e1 * t241 - t298 / 0.2e1) * t137) * t136 + (qJD(5) * t165 * t280 + t219 * t75 - t27 * t338) * t132) * t216 + (-t56 * t318 + t28 * t319 - t38 * t241 + ((-qJD(3) * t21 + t59) * t139 - t135 * t11 + (-t135 * t20 - t139 * t21) * qJD(2)) * t132 + t4 * t283) * t133 + (t132 * t4 + t318 * t38) * t286 - t210 * t293 + t270 + t365 * (-t110 * t318 - t309 * t75) + t364 * (-t165 * t206 + t221 * t75); t155 + ((-t103 * t282 - t259) * t136 + (t194 + t103 * t287 + (t115 + t349) * qJD(5)) * t132) * t106 + t150; t150 + t79 + t357; (-t27 * t311 + t161 * t309 + (t209 * t75 + t297 * t318) * qJD(5)) * t322 + (t311 * t75 - t190) * t178 - t136 * t165 * t171 - (t26 * t292 + (t63 + (-qJD(2) * t87 - t78) * t135) * t106 + ((-t136 * t290 - t224) * t139 + (-t300 * t352 + t228) * t135) * t81) * t303 + t75 * t57 - (t177 + t367) * t27 + t354 * (t106 * t336 + t26 * t135 + t139 * t252) + (-t28 * t204 + (t246 * t324 - t172) * t165) * t132; ((t285 * t81 + t193) * t133 + t36 * t286) * t132 + t36 * t226 + t362; -t50 * t226 - t132 * (-t229 * t81 + t286 * t50) + t362; t57 + (qJD(5) * t209 * t94 - 0.2e1 * t107 * t311 * t329) * t131; t165 * t188 + t28 * t255 - t165 * t164 + ((t165 * t282 + t27) * t132 + t154 * t136) * t262 + (t28 * t90 + (-(-t135 * t78 + t63 + (-t135 * t87 + t139 * t82) * qJD(2)) * t136 + t34 * t285) * t106) * t133 + (t165 * t90 - t315 * t34) * t286 - pkin(7) * t146 + (t256 - t257) * (-t132 * t75 - t165 * t302); (-qJD(5) * t317 + t189) * t136 + (t102 * t132 * t282 + (-t136 * t281 + (t132 * t134 - t245) * qJD(3)) * pkin(3)) * t106 + (-(t90 - t248) * t133 + t153) * qJD(4) + t212; t94 * t187 + (-pkin(5) * t247 - t136 * t317) * qJD(5) + (-(t90 + t263) * t133 + t153) * qJD(4) + t212; pkin(7) * t250 + (-qJD(5) * t261 + t312 * t324) * t131; t178; t146; t58; t58; -t250; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11), t5(16); t5(2), t5(3), t5(5), t5(8), t5(12), t5(17); t5(4), t5(5), t5(6), t5(9), t5(13), t5(18); t5(7), t5(8), t5(9), t5(10), t5(14), t5(19); t5(11), t5(12), t5(13), t5(14), t5(15), t5(20); t5(16), t5(17), t5(18), t5(19), t5(20), t5(21);];
Mq = res;
