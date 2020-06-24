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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 21:55:58
% EndTime: 2020-06-23 21:56:04
% DurationCPUTime: 6.14s
% Computational Cost: add. (11198->403), mult. (17242->624), div. (0->0), fcn. (20538->12), ass. (0->252)
t313 = pkin(6) + rSges(7,3);
t219 = qJ(2) + qJ(3);
t213 = sin(t219);
t227 = cos(qJ(4));
t229 = cos(qJ(1));
t268 = t227 * t229;
t222 = sin(qJ(4));
t224 = sin(qJ(1));
t274 = t222 * t224;
t173 = t213 * t274 + t268;
t270 = t224 * t227;
t273 = t222 * t229;
t176 = t213 * t273 - t270;
t214 = cos(t219);
t278 = t214 * t222;
t175 = t213 * t270 - t273;
t221 = sin(qJ(5));
t226 = cos(qJ(5));
t277 = t214 * t224;
t142 = t175 * t226 + t221 * t277;
t178 = t213 * t268 + t274;
t275 = t214 * t229;
t144 = t178 * t226 + t221 * t275;
t282 = t173 * t176;
t141 = -t175 * t221 + t226 * t277;
t143 = -t178 * t221 + t226 * t275;
t287 = t141 * t143;
t62 = Icges(6,1) * t142 * t144 + Icges(6,2) * t287 + Icges(6,3) * t282;
t261 = Icges(6,3) * t278;
t276 = t214 * t227;
t169 = -t213 * t226 - t221 * t276;
t286 = t141 * t169;
t170 = -t213 * t221 + t226 * t276;
t295 = Icges(6,1) * t170;
t82 = Icges(6,2) * t286 + t142 * t295 + t173 * t261;
t137 = t141 ^ 2;
t324 = t173 ^ 2;
t92 = Icges(6,1) * t142 ^ 2 + Icges(6,2) * t137 + Icges(6,3) * t324;
t27 = t173 * t92 + t176 * t62 + t82 * t278;
t220 = sin(qJ(6));
t225 = cos(qJ(6));
t117 = t142 * t220 + t173 * t225;
t118 = -t142 * t225 + t173 * t220;
t119 = t144 * t220 + t176 * t225;
t120 = -t144 * t225 + t176 * t220;
t35 = Icges(7,1) * t118 * t120 + Icges(7,2) * t117 * t119 + Icges(7,3) * t287;
t133 = t170 * t220 + t225 * t278;
t289 = Icges(7,2) * t133;
t134 = -t170 * t225 + t220 * t278;
t294 = Icges(7,1) * t134;
t51 = Icges(7,3) * t286 + t117 * t289 + t118 * t294;
t67 = Icges(7,1) * t118 ^ 2 + Icges(7,2) * t117 ^ 2 + Icges(7,3) * t137;
t6 = t173 * t67 + t176 * t35 + t51 * t278;
t335 = t6 + t27;
t285 = t143 * t169;
t83 = Icges(6,2) * t285 + t144 * t295 + t176 * t261;
t138 = t143 ^ 2;
t323 = t176 ^ 2;
t93 = Icges(6,1) * t144 ^ 2 + Icges(6,2) * t138 + Icges(6,3) * t323;
t28 = t173 * t62 + t176 * t93 + t83 * t278;
t52 = Icges(7,3) * t285 + t119 * t289 + t120 * t294;
t68 = Icges(7,1) * t120 ^ 2 + Icges(7,2) * t119 ^ 2 + Icges(7,3) * t138;
t7 = t173 * t35 + t176 * t68 + t52 * t278;
t334 = t7 + t28;
t50 = t51 * t229;
t17 = -t52 * t224 + t50;
t79 = t82 * t229;
t333 = -t83 * t224 + t17 + t79;
t303 = t224 * t35;
t23 = t229 * t67 - t303;
t302 = t224 * t62;
t332 = t229 * t92 + t23 - t302;
t299 = t229 * t35;
t24 = -t224 * t68 + t299;
t298 = t229 * t62;
t331 = -t224 * t93 + t24 + t298;
t10 = -t51 * t213 + (t224 * t67 + t299) * t214;
t280 = t213 * t224;
t290 = Icges(5,2) * t222;
t296 = Icges(5,1) * t227;
t103 = (-Icges(5,3) * t280 + t173 * t290 + t175 * t296) * t214;
t216 = t224 ^ 2;
t212 = t214 ^ 2;
t288 = Icges(5,3) * t212;
t123 = Icges(5,1) * t175 ^ 2 + Icges(5,2) * t324 + t216 * t288;
t269 = t224 * t229;
t94 = Icges(5,1) * t175 * t178 + Icges(5,2) * t282 + t269 * t288;
t297 = t229 * t94;
t31 = -t82 * t213 + (t224 * t92 + t298) * t214;
t330 = t10 + t31 - t103 * t213 + (t123 * t224 + t297) * t214;
t279 = t213 * t229;
t104 = (-Icges(5,3) * t279 + t176 * t290 + t178 * t296) * t214;
t11 = -t52 * t213 + (t229 * t68 + t303) * t214;
t218 = t229 ^ 2;
t124 = Icges(5,1) * t178 ^ 2 + Icges(5,2) * t323 + t218 * t288;
t301 = t224 * t94;
t32 = -t83 * t213 + (t229 * t93 + t302) * t214;
t329 = t11 + t32 - t104 * t213 + (t124 * t229 + t301) * t214;
t328 = t123 * t229 - t301 + t332;
t327 = -t124 * t224 + t297 + t331;
t326 = (t141 * t229 - t143 * t224) * t214;
t325 = t313 * t143;
t322 = m(3) * rSges(3,1);
t321 = t141 / 0.2e1;
t320 = t143 / 0.2e1;
t319 = t169 / 0.2e1;
t318 = t173 / 0.2e1;
t317 = t176 / 0.2e1;
t316 = -t213 / 0.2e1;
t315 = -t224 / 0.2e1;
t314 = t229 / 0.2e1;
t312 = pkin(2) * t229;
t311 = pkin(4) * t213;
t168 = t169 ^ 2;
t90 = Icges(7,1) * t134 ^ 2 + Icges(7,2) * t133 ^ 2 + Icges(7,3) * t168;
t13 = t51 * t173 + t52 * t176 + t90 * t278;
t281 = t212 * t222 ^ 2;
t122 = Icges(6,1) * t170 ^ 2 + Icges(6,2) * t168 + Icges(6,3) * t281;
t33 = t122 * t278 + t82 * t173 + t83 * t176;
t310 = t13 + t33;
t309 = t52 * t275 + t51 * t277;
t308 = t83 * t275 + t82 * t277;
t307 = t103 * t277 + t104 * t275;
t223 = sin(qJ(2));
t306 = rSges(3,1) * t223;
t305 = rSges(6,2) * t143;
t304 = t141 * rSges(6,2);
t300 = t229 * rSges(4,3);
t293 = Icges(4,4) * t213;
t292 = Icges(4,4) * t214;
t228 = cos(qJ(2));
t217 = t228 ^ 2;
t291 = Icges(3,2) * t217;
t284 = t169 * t224;
t283 = t169 * t229;
t272 = t223 * t224;
t271 = t223 * t229;
t267 = t228 * t229;
t252 = pkin(5) * t214 + t311;
t171 = t252 * t224;
t211 = pkin(3) * t223 + pkin(2);
t183 = (-pkin(2) + t211) * t224;
t266 = -t171 - t183;
t172 = pkin(4) * t279 + pkin(5) * t275;
t201 = t229 * t211;
t184 = t201 - t312;
t265 = -t172 - t184;
t190 = pkin(4) * t214 - pkin(5) * t213;
t181 = t224 * t190;
t127 = rSges(6,2) * t284 + t181;
t182 = t229 * t190;
t128 = rSges(6,2) * t283 + t182;
t264 = t216 + t218;
t263 = -t104 - t52 - t83;
t148 = t212 * t227 ^ 2 * Icges(5,1) + t213 ^ 2 * Icges(5,3) + Icges(5,2) * t281;
t262 = -t90 - t122 - t148;
t109 = t313 * t284 + t181;
t110 = t313 * t283 + t182;
t260 = t201 + t172;
t259 = t313 * t141;
t258 = t278 / 0.2e1;
t246 = Icges(4,5) * t213 + Icges(4,6) * t214;
t151 = Icges(4,3) * t229 + t246 * t224;
t152 = -Icges(4,3) * t224 + t246 * t229;
t247 = Icges(4,2) * t214 + t293;
t154 = -Icges(4,6) * t224 + t247 * t229;
t248 = Icges(4,1) * t213 + t292;
t156 = -Icges(4,5) * t224 + t248 * t229;
t242 = -t154 * t214 - t156 * t213;
t153 = Icges(4,6) * t229 + t247 * t224;
t155 = Icges(4,5) * t229 + t248 * t224;
t243 = t153 * t214 + t155 * t213;
t255 = (t218 * t151 + (t242 * t224 + (-t152 + t243) * t229) * t224 + t328) * t229;
t254 = -t216 * t152 - (t243 * t229 + (-t151 + t242) * t224) * t229 - t327;
t253 = t264 * t214 * rSges(5,3);
t146 = -rSges(5,3) * t280 + t181;
t147 = -rSges(5,3) * t279 + t182;
t164 = rSges(4,1) * t279 + rSges(4,2) * t275 - rSges(4,3) * t224;
t251 = -t171 - t259;
t250 = -t172 - t325;
t249 = rSges(4,1) * t213 + rSges(4,2) * t214;
t245 = t141 * t176 - t143 * t173;
t241 = -t224 * t171 - t229 * t172;
t187 = -Icges(4,2) * t213 + t292;
t188 = Icges(4,1) * t214 - t293;
t240 = t187 * t214 + t188 * t213;
t239 = -t141 * t278 + t169 * t173;
t238 = t141 * t213 + t169 * t277;
t237 = t143 * t278 - t169 * t176;
t236 = -t143 * t213 - t169 * t275;
t102 = t103 * t229;
t186 = Icges(4,5) * t214 - Icges(4,6) * t213;
t235 = t102 + t50 + t79 + (-t154 * t213 + t156 * t214 - t224 * t186 + t240 * t229) * t315 + (-t153 * t213 + t155 * t214 + t229 * t186 + t240 * t224) * t314;
t234 = (-t211 - t252) * t224;
t233 = t333 * t258 + t335 * t314 + t334 * t315 + t331 * t317 + t332 * t318;
t232 = t254 * t224 + t255;
t231 = (-t104 * t224 + t102 + t333) * t316 + t329 * t315 + t330 * t314 + t328 * t277 / 0.2e1 + t327 * t275 / 0.2e1;
t208 = pkin(3) * t267;
t207 = t224 * t228 * pkin(3);
t206 = rSges(3,1) * t271;
t197 = rSges(2,1) * t229 - rSges(2,2) * t224;
t196 = -rSges(2,1) * t224 - rSges(2,2) * t229;
t192 = Icges(3,5) * t271 - Icges(3,3) * t224;
t191 = Icges(3,5) * t272 + Icges(3,3) * t229;
t189 = rSges(4,1) * t214 - rSges(4,2) * t213;
t180 = -rSges(3,3) * t224 + t206 + t312;
t179 = -t229 * rSges(3,3) + (-pkin(2) - t306) * t224;
t163 = t249 * t224 + t300;
t150 = t189 * t229 + t208;
t149 = t189 * t224 + t207;
t145 = -t206 * t229 - t216 * t306;
t140 = t164 + t201;
t139 = -t300 + (-t211 - t249) * t224;
t136 = t147 + t208;
t135 = t146 + t207;
t130 = rSges(5,3) * t275 + t260;
t129 = (-t311 - t211 + (-pkin(5) - rSges(5,3)) * t214) * t224;
t126 = t208 + t128;
t125 = t207 + t127;
t121 = -t224 * t163 - t229 * t164;
t115 = t260 + t305;
t114 = t234 - t304;
t113 = -t253 + t241;
t106 = t236 * rSges(6,2);
t105 = t238 * rSges(6,2);
t101 = t208 + t110;
t100 = t207 + t109;
t97 = t237 * rSges(6,2);
t96 = t239 * rSges(6,2);
t95 = rSges(6,2) * t326;
t91 = (-t164 - t184) * t229 + (-t163 - t183) * t224;
t89 = t260 + t325;
t88 = -t259 + t234;
t87 = t266 * t224 + t265 * t229 - t253;
t85 = t245 * rSges(6,2);
t84 = (-t141 * t224 - t143 * t229) * rSges(6,2) + t241;
t81 = t313 * t236;
t80 = t313 * t238;
t66 = (t265 - t305) * t229 + (t266 - t304) * t224;
t65 = t313 * t237;
t64 = t313 * t239;
t63 = t313 * t326;
t60 = t251 * t224 + t250 * t229;
t58 = t313 * t245;
t55 = (-t184 + t250) * t229 + (-t183 + t251) * t224;
t34 = -t122 * t213 + t308;
t14 = -t90 * t213 + t309;
t12 = t141 * t51 + t143 * t52 + t169 * t90;
t3 = t141 * t35 + t143 * t68 + t169 * t52;
t2 = t141 * t67 + t143 * t35 + t169 * t51;
t1 = t17 * t319 + t2 * t314 + t23 * t321 + t24 * t320 + t3 * t315;
t4 = [t223 ^ 2 * Icges(3,2) + Icges(3,1) * t217 - t213 * t187 + t214 * t188 + Icges(2,3) + m(6) * (t114 ^ 2 + t115 ^ 2) + m(7) * (t88 ^ 2 + t89 ^ 2) + m(4) * (t139 ^ 2 + t140 ^ 2) + m(5) * (t129 ^ 2 + t130 ^ 2) + m(3) * (t179 ^ 2 + t180 ^ 2) + m(2) * (t196 ^ 2 + t197 ^ 2) - t262; m(6) * (t114 * t126 + t115 * t125) + m(7) * (t100 * t89 + t101 * t88) + m(4) * (t139 * t150 + t140 * t149) + m(5) * (t129 * t136 + t130 * t135) + ((-Icges(3,1) * t271 / 0.2e1 + Icges(3,5) * t224 + t180 * t322) * t228 + t263) * t224 + t235 + (Icges(3,1) * t272 / 0.2e1 + Icges(3,5) * t229 + t179 * t322) * t267; (t191 * t229 + t216 * t291) * t218 + m(7) * (t100 ^ 2 + t101 ^ 2 + t55 ^ 2) + m(6) * (t125 ^ 2 + t126 ^ 2 + t66 ^ 2) + m(5) * (t135 ^ 2 + t136 ^ 2 + t87 ^ 2) + m(4) * (t149 ^ 2 + t150 ^ 2 + t91 ^ 2) + m(3) * (rSges(3,1) ^ 2 * t217 * t264 + t145 ^ 2) + ((-t192 * t224 + t218 * t291) * t224 + (t191 * t224 - t192 * t229 - 0.2e1 * t269 * t291) * t229 + t254) * t224 + t255; t263 * t224 + m(6) * (t114 * t128 + t115 * t127) + m(7) * (t109 * t89 + t110 * t88) + m(5) * (t129 * t147 + t130 * t146) + m(4) * (t139 * t229 + t140 * t224) * t189 + t235; m(7) * (t100 * t109 + t101 * t110 + t55 * t60) + m(6) * (t125 * t127 + t126 * t128 + t66 * t84) + m(5) * (t113 * t87 + t135 * t146 + t136 * t147) + m(4) * (t121 * t91 + (t149 * t224 + t150 * t229) * t189) + t232; m(7) * (t109 ^ 2 + t110 ^ 2 + t60 ^ 2) + m(6) * (t127 ^ 2 + t128 ^ 2 + t84 ^ 2) + m(5) * (t113 ^ 2 + t146 ^ 2 + t147 ^ 2) + m(4) * (t189 ^ 2 * t264 + t121 ^ 2) + t232; m(6) * (t105 * t114 + t106 * t115) + m(7) * (t80 * t88 + t81 * t89) + t262 * t213 + t307 + t308 + t309; m(7) * (t100 * t81 + t101 * t80 + t55 * t63) + m(6) * (t105 * t126 + t106 * t125 + t66 * t95) + t231; m(7) * (t109 * t81 + t110 * t80 + t60 * t63) + m(6) * (t105 * t128 + t106 * t127 + t84 * t95) + t231; m(7) * (t63 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(6) * (t105 ^ 2 + t106 ^ 2 + t95 ^ 2) + (t148 * t213 - t14 - t307 - t34) * t213 + (t330 * t224 + t329 * t229) * t214; m(6) * (t114 * t96 + t115 * t97) + m(7) * (t64 * t88 + t65 * t89) + t310; m(7) * (t100 * t65 + t101 * t64 + t55 * t58) + m(6) * (t125 * t97 + t126 * t96 + t66 * t85) + t233; m(7) * (t109 * t65 + t110 * t64 + t58 * t60) + m(6) * (t127 * t97 + t128 * t96 + t84 * t85) + t233; (-t13 / 0.2e1 - t33 / 0.2e1) * t213 + (t11 / 0.2e1 + t32 / 0.2e1) * t176 + (t10 / 0.2e1 + t31 / 0.2e1) * t173 + m(7) * (t58 * t63 + t64 * t80 + t65 * t81) + m(6) * (t105 * t96 + t106 * t97 + t85 * t95) + ((t7 / 0.2e1 + t28 / 0.2e1) * t229 + (t6 / 0.2e1 + t27 / 0.2e1) * t224 + (t14 / 0.2e1 + t34 / 0.2e1) * t222) * t214; t310 * t278 + t334 * t176 + t335 * t173 + m(7) * (t58 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(6) * (t85 ^ 2 + t96 ^ 2 + t97 ^ 2); t12; t1; t1; t10 * t321 + t14 * t319 + t11 * t320 + t12 * t316 + (t3 * t314 + t224 * t2 / 0.2e1) * t214; t12 * t258 + t13 * t319 + t2 * t318 + t3 * t317 + t320 * t7 + t321 * t6; t12 * t169 + t141 * t2 + t143 * t3;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11), t4(16); t4(2), t4(3), t4(5), t4(8), t4(12), t4(17); t4(4), t4(5), t4(6), t4(9), t4(13), t4(18); t4(7), t4(8), t4(9), t4(10), t4(14), t4(19); t4(11), t4(12), t4(13), t4(14), t4(15), t4(20); t4(16), t4(17), t4(18), t4(19), t4(20), t4(21);];
Mq = res;
