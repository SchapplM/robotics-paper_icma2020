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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350OL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350OL_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350OL_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:00:14
% EndTime: 2020-06-20 08:00:25
% DurationCPUTime: 11.55s
% Computational Cost: add. (28419->591), mult. (40855->837), div. (0->0), fcn. (51308->12), ass. (0->283)
t282 = qJ(2) + qJ(3);
t278 = sin(t282);
t290 = cos(qJ(4));
t292 = cos(qJ(1));
t332 = t292 * t290;
t285 = sin(qJ(4));
t287 = sin(qJ(1));
t336 = t287 * t285;
t245 = t278 * t336 + t332;
t333 = t292 * t285;
t335 = t287 * t290;
t247 = t278 * t333 - t335;
t279 = cos(t282);
t341 = t279 * t285;
t246 = t278 * t335 - t333;
t284 = sin(qJ(5));
t289 = cos(qJ(5));
t340 = t279 * t287;
t208 = -t246 * t284 + t289 * t340;
t209 = t246 * t289 + t284 * t340;
t151 = Icges(6,5) * t209 + Icges(6,6) * t208 + Icges(6,3) * t245;
t153 = Icges(6,4) * t209 + Icges(6,2) * t208 + Icges(6,6) * t245;
t155 = Icges(6,1) * t209 + Icges(6,4) * t208 + Icges(6,5) * t245;
t81 = t245 * t151 + t208 * t153 + t209 * t155;
t248 = t278 * t332 + t336;
t338 = t279 * t292;
t210 = -t248 * t284 + t289 * t338;
t211 = t248 * t289 + t284 * t338;
t152 = Icges(6,5) * t211 + Icges(6,6) * t210 + Icges(6,3) * t247;
t154 = Icges(6,4) * t211 + Icges(6,2) * t210 + Icges(6,6) * t247;
t156 = Icges(6,1) * t211 + Icges(6,4) * t210 + Icges(6,5) * t247;
t82 = t245 * t152 + t208 * t154 + t209 * t156;
t339 = t279 * t290;
t238 = -t278 * t289 - t284 * t339;
t239 = -t278 * t284 + t289 * t339;
t185 = Icges(6,5) * t239 + Icges(6,6) * t238 + Icges(6,3) * t341;
t186 = Icges(6,4) * t239 + Icges(6,2) * t238 + Icges(6,6) * t341;
t187 = Icges(6,1) * t239 + Icges(6,4) * t238 + Icges(6,5) * t341;
t96 = t245 * t185 + t208 * t186 + t209 * t187;
t31 = t81 * t245 + t82 * t247 + t96 * t341;
t283 = sin(qJ(6));
t288 = cos(qJ(6));
t176 = t209 * t283 + t245 * t288;
t177 = -t209 * t288 + t245 * t283;
t120 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t208;
t122 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t208;
t124 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t208;
t49 = t208 * t120 + t176 * t122 + t177 * t124;
t178 = t211 * t283 + t247 * t288;
t179 = -t211 * t288 + t247 * t283;
t121 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t210;
t123 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t210;
t125 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t210;
t50 = t208 * t121 + t176 * t123 + t177 * t125;
t200 = t239 * t283 + t288 * t341;
t201 = -t239 * t288 + t283 * t341;
t146 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t238;
t147 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t238;
t148 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t238;
t65 = t208 * t146 + t176 * t147 + t177 * t148;
t7 = t49 * t245 + t50 * t247 + t65 * t341;
t378 = t7 + t31;
t83 = t247 * t151 + t210 * t153 + t211 * t155;
t84 = t247 * t152 + t210 * t154 + t211 * t156;
t97 = t247 * t185 + t210 * t186 + t211 * t187;
t32 = t83 * t245 + t84 * t247 + t97 * t341;
t51 = t210 * t120 + t178 * t122 + t179 * t124;
t52 = t210 * t121 + t178 * t123 + t179 * t125;
t66 = t210 * t146 + t178 * t147 + t179 * t148;
t8 = t51 * t245 + t52 * t247 + t66 * t341;
t377 = t8 + t32;
t23 = -t50 * t287 + t49 * t292;
t376 = -t82 * t287 + t81 * t292 + t23;
t24 = -t52 * t287 + t51 * t292;
t375 = -t84 * t287 + t83 * t292 + t24;
t56 = t238 * t121 + t200 * t123 + t201 * t125;
t351 = t56 * t287;
t55 = t238 * t120 + t200 * t122 + t201 * t124;
t352 = t55 * t292;
t28 = -t351 + t352;
t90 = t152 * t341 + t238 * t154 + t239 * t156;
t349 = t90 * t287;
t89 = t151 * t341 + t238 * t153 + t239 * t155;
t350 = t89 * t292;
t374 = t28 - t349 + t350;
t189 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t340;
t191 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t340;
t193 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t340;
t106 = t189 * t340 - t245 * t191 + t246 * t193;
t190 = Icges(5,5) * t248 - Icges(5,6) * t247 + Icges(5,3) * t338;
t192 = Icges(5,4) * t248 - Icges(5,2) * t247 + Icges(5,6) * t338;
t194 = Icges(5,1) * t248 - Icges(5,4) * t247 + Icges(5,5) * t338;
t107 = t190 * t340 - t245 * t192 + t246 * t194;
t11 = -t65 * t278 + (t287 * t49 + t292 * t50) * t279;
t214 = -Icges(5,3) * t278 + (Icges(5,5) * t290 - Icges(5,6) * t285) * t279;
t215 = -Icges(5,6) * t278 + (Icges(5,4) * t290 - Icges(5,2) * t285) * t279;
t216 = -Icges(5,5) * t278 + (Icges(5,1) * t290 - Icges(5,4) * t285) * t279;
t137 = t214 * t340 - t245 * t215 + t246 * t216;
t35 = -t96 * t278 + (t287 * t81 + t292 * t82) * t279;
t373 = t11 + t35 - t137 * t278 + (t106 * t287 + t107 * t292) * t279;
t108 = t189 * t338 - t247 * t191 + t248 * t193;
t109 = t190 * t338 - t247 * t192 + t248 * t194;
t12 = -t66 * t278 + (t287 * t51 + t292 * t52) * t279;
t138 = t214 * t338 - t247 * t215 + t248 * t216;
t36 = -t97 * t278 + (t287 * t83 + t292 * t84) * t279;
t372 = t12 + t36 - t138 * t278 + (t108 * t287 + t109 * t292) * t279;
t371 = t106 * t292 - t107 * t287 + t376;
t370 = t108 * t292 - t109 * t287 + t375;
t103 = t185 * t341 + t238 * t186 + t239 * t187;
t76 = t238 * t146 + t200 * t147 + t201 * t148;
t369 = -t103 - t76;
t127 = t179 * rSges(7,1) + t178 * rSges(7,2) + t210 * rSges(7,3);
t368 = -t210 * pkin(6) - t127;
t305 = Icges(4,5) * t278 + Icges(4,6) * t279;
t220 = Icges(4,3) * t292 + t305 * t287;
t221 = -Icges(4,3) * t287 + t305 * t292;
t280 = t287 ^ 2;
t346 = Icges(4,4) * t278;
t307 = Icges(4,2) * t279 + t346;
t223 = -Icges(4,6) * t287 + t307 * t292;
t345 = Icges(4,4) * t279;
t309 = Icges(4,1) * t278 + t345;
t225 = -Icges(4,5) * t287 + t309 * t292;
t303 = -t223 * t279 - t225 * t278;
t222 = Icges(4,6) * t292 + t307 * t287;
t224 = Icges(4,5) * t292 + t309 * t287;
t304 = t222 * t279 + t224 * t278;
t367 = -t280 * t221 - (t304 * t292 + (-t220 + t303) * t287) * t292 - t370;
t281 = t292 ^ 2;
t366 = t208 / 0.2e1;
t365 = t210 / 0.2e1;
t364 = t238 / 0.2e1;
t363 = t245 / 0.2e1;
t362 = t247 / 0.2e1;
t361 = -t278 / 0.2e1;
t360 = -t287 / 0.2e1;
t359 = t292 / 0.2e1;
t358 = pkin(4) * t278;
t357 = pkin(6) * t238;
t355 = t292 * pkin(2);
t286 = sin(qJ(2));
t354 = rSges(3,1) * t286;
t353 = t292 * rSges(4,3);
t348 = Icges(3,4) * t286;
t291 = cos(qJ(2));
t347 = Icges(3,4) * t291;
t114 = -t278 * t189 + (-t191 * t285 + t193 * t290) * t279;
t344 = t114 * t292;
t115 = -t278 * t190 + (-t192 * t285 + t194 * t290) * t279;
t343 = t115 * t287;
t342 = t278 * t292;
t337 = t285 * t215;
t149 = t201 * rSges(7,1) + t200 * rSges(7,2) + t238 * rSges(7,3);
t144 = t287 * t149;
t188 = t239 * rSges(6,1) + t238 * rSges(6,2) + rSges(6,3) * t341;
t180 = t287 * t188;
t217 = -t278 * rSges(5,3) + (rSges(5,1) * t290 - rSges(5,2) * t285) * t279;
t212 = t287 * t217;
t334 = t291 * t292;
t145 = t292 * t149;
t181 = t292 * t188;
t213 = t292 * t217;
t315 = pkin(5) * t279 + t358;
t243 = t315 * t287;
t277 = t286 * pkin(3) + pkin(2);
t251 = (-pkin(2) + t277) * t287;
t331 = -t243 - t251;
t244 = pkin(4) * t342 + pkin(5) * t338;
t266 = t292 * t277;
t252 = t266 - t355;
t330 = -t244 - t252;
t258 = t279 * pkin(4) - t278 * pkin(5);
t249 = t287 * t258;
t274 = t287 * t291 * pkin(3);
t329 = t249 + t274;
t250 = t292 * t258;
t275 = pkin(3) * t334;
t328 = t250 + t275;
t327 = rSges(3,2) * t334 + t292 * t354;
t326 = t280 + t281;
t325 = t55 / 0.2e1 + t65 / 0.2e1;
t324 = t66 / 0.2e1 + t56 / 0.2e1;
t158 = t211 * rSges(6,1) + t210 * rSges(6,2) + t247 * rSges(6,3);
t196 = t248 * rSges(5,1) - t247 * rSges(5,2) + rSges(5,3) * t338;
t131 = t287 * t357 + t144 + t249;
t132 = t292 * t357 + t145 + t250;
t323 = t266 + t244;
t322 = t341 / 0.2e1;
t319 = (t281 * t220 + (t303 * t287 + (-t221 + t304) * t292) * t287 + t371) * t292;
t3 = t49 * t208 + t50 * t210 + t65 * t238;
t4 = t51 * t208 + t52 * t210 + t66 * t238;
t318 = t23 * t366 + t24 * t365 + t28 * t364 + t3 * t359 + t4 * t360;
t311 = -t177 * rSges(7,1) - t176 * rSges(7,2);
t126 = t208 * rSges(7,3) - t311;
t317 = -t208 * pkin(6) - t126 - t243;
t316 = -t244 + t368;
t231 = rSges(4,1) * t342 + rSges(4,2) * t338 - t287 * rSges(4,3);
t314 = -rSges(3,2) * t291 - t354;
t313 = rSges(4,1) * t278 + rSges(4,2) * t279;
t312 = -t246 * rSges(5,1) + t245 * rSges(5,2);
t310 = Icges(3,1) * t286 + t347;
t308 = Icges(3,2) * t291 + t348;
t306 = Icges(3,5) * t286 + Icges(3,6) * t291;
t255 = -Icges(4,2) * t278 + t345;
t256 = Icges(4,1) * t279 - t346;
t300 = t255 * t279 + t256 * t278;
t299 = t89 / 0.2e1 + t96 / 0.2e1 + t325;
t298 = t97 / 0.2e1 + t90 / 0.2e1 + t324;
t297 = (-t277 - t315) * t287;
t157 = t209 * rSges(6,1) + t208 * rSges(6,2) + t245 * rSges(6,3);
t296 = t374 * t322 + t378 * t359 + t377 * t360 + t375 * t362 + t376 * t363;
t295 = t367 * t287 + t319;
t294 = (-t343 + t344 + t374) * t361 + t372 * t360 + t373 * t359 + t371 * t340 / 0.2e1 + t370 * t338 / 0.2e1;
t254 = Icges(4,5) * t279 - Icges(4,6) * t278;
t293 = t352 / 0.2e1 - t351 / 0.2e1 + t350 / 0.2e1 - t349 / 0.2e1 - t343 / 0.2e1 + t344 / 0.2e1 + (-t278 * t223 + t279 * t225 - t287 * t254 + t300 * t292 + t138 + t66 + t97) * t360 + (-t278 * t222 + t279 * t224 + t292 * t254 + t300 * t287 + t137 + t65 + t96) * t359;
t265 = t292 * rSges(2,1) - t287 * rSges(2,2);
t264 = t291 * rSges(3,1) - t286 * rSges(3,2);
t263 = -t287 * rSges(2,1) - t292 * rSges(2,2);
t257 = t279 * rSges(4,1) - t278 * rSges(4,2);
t233 = -Icges(3,3) * t287 + t306 * t292;
t232 = Icges(3,3) * t292 + t306 * t287;
t230 = t313 * t287 + t353;
t227 = -t287 * rSges(3,3) + t327 + t355;
t226 = -t292 * rSges(3,3) + (-pkin(2) + t314) * t287;
t219 = t292 * t257 + t275;
t218 = t287 * t257 + t274;
t207 = t231 + t266;
t206 = -t353 + (-t277 - t313) * t287;
t205 = t216 * t339;
t199 = t250 + t213;
t198 = t249 + t212;
t197 = t314 * t280 - t292 * t327;
t195 = rSges(5,3) * t340 - t312;
t184 = t213 + t328;
t183 = t212 + t329;
t182 = -t287 * t230 - t292 * t231;
t166 = t250 + t181;
t165 = t249 + t180;
t164 = t323 + t196;
t163 = (-t358 - t277 + (-pkin(5) - rSges(5,3)) * t279) * t287 + t312;
t162 = t181 + t328;
t161 = t180 + t329;
t160 = -t278 * t196 - t279 * t213;
t159 = t278 * t195 + t279 * t212;
t150 = (-t231 - t252) * t292 + (-t230 - t251) * t287;
t143 = -t278 * t214 - t279 * t337 + t205;
t142 = (t195 * t292 - t196 * t287) * t279;
t140 = t323 + t158;
t139 = t297 - t157;
t130 = t275 + t132;
t129 = t274 + t131;
t128 = (-t196 - t244) * t292 + (-t195 - t243) * t287;
t119 = -t278 * t158 - t279 * t181;
t118 = t278 * t157 + t279 * t180;
t117 = t158 * t341 - t247 * t188;
t116 = -t157 * t341 + t245 * t188;
t111 = (-t196 + t330) * t292 + (-t195 + t331) * t287;
t105 = (t157 * t292 - t158 * t287) * t279;
t102 = t247 * t157 - t245 * t158;
t101 = t323 - t368;
t100 = (-pkin(6) - rSges(7,3)) * t208 + t297 + t311;
t99 = t103 * t341;
t98 = (-t158 - t244) * t292 + (-t157 - t243) * t287;
t93 = (-t158 + t330) * t292 + (-t157 + t331) * t287;
t92 = t238 * t127 - t210 * t149;
t91 = -t238 * t126 + t208 * t149;
t86 = -t279 * t145 - t278 * t127 + (-t210 * t278 - t238 * t338) * pkin(6);
t85 = t279 * t144 + t278 * t126 + (t208 * t278 + t238 * t340) * pkin(6);
t79 = t127 * t341 - t247 * t149 + (t210 * t341 - t238 * t247) * pkin(6);
t78 = -t126 * t341 + t245 * t149 + (-t208 * t341 + t238 * t245) * pkin(6);
t77 = t210 * t126 - t208 * t127;
t75 = t76 * t341;
t71 = t76 * t238;
t68 = (t126 * t292 - t127 * t287 + (t208 * t292 - t210 * t287) * pkin(6)) * t279;
t67 = t317 * t287 + t316 * t292;
t62 = t247 * t126 - t245 * t127 + (t208 * t247 - t210 * t245) * pkin(6);
t61 = (-t252 + t316) * t292 + (-t251 + t317) * t287;
t38 = -t103 * t278 + (t89 * t287 + t90 * t292) * t279;
t37 = t89 * t245 + t90 * t247 + t99;
t15 = -t76 * t278 + (t55 * t287 + t56 * t292) * t279;
t14 = t55 * t245 + t56 * t247 + t75;
t13 = t55 * t208 + t56 * t210 + t71;
t1 = [-t286 * (-Icges(3,2) * t286 + t347) + t291 * (Icges(3,1) * t291 - t348) + Icges(2,3) + t205 + (t256 - t337) * t279 + (-t214 - t255) * t278 + m(7) * (t100 ^ 2 + t101 ^ 2) + m(6) * (t139 ^ 2 + t140 ^ 2) + m(5) * (t163 ^ 2 + t164 ^ 2) + m(4) * (t206 ^ 2 + t207 ^ 2) + m(3) * (t226 ^ 2 + t227 ^ 2) + m(2) * (t263 ^ 2 + t265 ^ 2) - t369; (t280 / 0.2e1 + t281 / 0.2e1) * (Icges(3,5) * t291 - Icges(3,6) * t286) + m(7) * (t130 * t100 + t129 * t101) + m(6) * (t162 * t139 + t161 * t140) + m(5) * (t184 * t163 + t183 * t164) + m(4) * (t219 * t206 + t218 * t207) + m(3) * (t226 * t292 + t227 * t287) * t264 + (-t286 * (-Icges(3,6) * t287 + t308 * t292) + t291 * (-Icges(3,5) * t287 + t310 * t292)) * t360 + (-t286 * (Icges(3,6) * t292 + t308 * t287) + t291 * (Icges(3,5) * t292 + t310 * t287)) * t359 + t293; m(7) * (t129 ^ 2 + t130 ^ 2 + t61 ^ 2) + m(6) * (t161 ^ 2 + t162 ^ 2 + t93 ^ 2) + m(5) * (t111 ^ 2 + t183 ^ 2 + t184 ^ 2) + m(4) * (t150 ^ 2 + t218 ^ 2 + t219 ^ 2) + t292 * t281 * t232 + m(3) * (t326 * t264 ^ 2 + t197 ^ 2) + t319 + (-t280 * t233 + (t287 * t232 - t292 * t233) * t292 + t367) * t287; m(4) * (t206 * t292 + t207 * t287) * t257 + t293 + m(7) * (t132 * t100 + t131 * t101) + m(6) * (t166 * t139 + t165 * t140) + m(5) * (t199 * t163 + t198 * t164); m(7) * (t131 * t129 + t132 * t130 + t67 * t61) + m(6) * (t165 * t161 + t166 * t162 + t98 * t93) + m(5) * (t128 * t111 + t198 * t183 + t199 * t184) + m(4) * (t182 * t150 + (t218 * t287 + t219 * t292) * t257) + t295; m(7) * (t131 ^ 2 + t132 ^ 2 + t67 ^ 2) + m(6) * (t165 ^ 2 + t166 ^ 2 + t98 ^ 2) + m(5) * (t128 ^ 2 + t198 ^ 2 + t199 ^ 2) + m(4) * (t326 * t257 ^ 2 + t182 ^ 2) + t295; (-t143 + t369) * t278 + m(7) * (t100 * t85 + t101 * t86) + m(6) * (t118 * t139 + t119 * t140) + m(5) * (t159 * t163 + t160 * t164) + ((t115 / 0.2e1 + t138 / 0.2e1 + t298) * t292 + (t137 / 0.2e1 + t114 / 0.2e1 + t299) * t287) * t279; m(7) * (t86 * t129 + t85 * t130 + t68 * t61) + m(6) * (t105 * t93 + t118 * t162 + t119 * t161) + m(5) * (t142 * t111 + t159 * t184 + t160 * t183) + t294; m(7) * (t86 * t131 + t85 * t132 + t68 * t67) + m(6) * (t105 * t98 + t118 * t166 + t119 * t165) + m(5) * (t142 * t128 + t159 * t199 + t160 * t198) + t294; (t143 * t278 - t15 - t38) * t278 + m(7) * (t68 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(6) * (t105 ^ 2 + t118 ^ 2 + t119 ^ 2) + m(5) * (t142 ^ 2 + t159 ^ 2 + t160 ^ 2) + ((-t278 * t115 + t372) * t292 + (-t278 * t114 + t373) * t287) * t279; t75 + t99 + m(7) * (t100 * t78 + t101 * t79) + m(6) * (t116 * t139 + t117 * t140) + t298 * t247 + t299 * t245; m(7) * (t79 * t129 + t78 * t130 + t62 * t61) + m(6) * (t102 * t93 + t116 * t162 + t117 * t161) + t296; m(7) * (t79 * t131 + t78 * t132 + t62 * t67) + m(6) * (t102 * t98 + t116 * t166 + t117 * t165) + t296; (-t14 / 0.2e1 - t37 / 0.2e1) * t278 + (t12 / 0.2e1 + t36 / 0.2e1) * t247 + (t11 / 0.2e1 + t35 / 0.2e1) * t245 + m(7) * (t62 * t68 + t78 * t85 + t79 * t86) + m(6) * (t102 * t105 + t116 * t118 + t117 * t119) + ((t8 / 0.2e1 + t32 / 0.2e1) * t292 + (t7 / 0.2e1 + t31 / 0.2e1) * t287 + (t15 / 0.2e1 + t38 / 0.2e1) * t285) * t279; (t14 + t37) * t341 + t377 * t247 + t378 * t245 + m(7) * (t62 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(6) * (t102 ^ 2 + t116 ^ 2 + t117 ^ 2); m(7) * (t100 * t91 + t101 * t92) + t71 + t324 * t210 + t325 * t208; m(7) * (t92 * t129 + t91 * t130 + t77 * t61) + t318; m(7) * (t92 * t131 + t91 * t132 + t77 * t67) + t318; t12 * t365 + t15 * t364 + t11 * t366 + m(7) * (t68 * t77 + t85 * t91 + t86 * t92) + t13 * t361 + (t4 * t359 + t287 * t3 / 0.2e1) * t279; t13 * t322 + t3 * t363 + m(7) * (t62 * t77 + t78 * t91 + t79 * t92) + t14 * t364 + t7 * t366 + t8 * t365 + t4 * t362; t210 * t4 + t208 * t3 + t238 * t13 + m(7) * (t77 ^ 2 + t91 ^ 2 + t92 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
