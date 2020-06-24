% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = CloosQRC350DE_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:02:09
% EndTime: 2020-06-23 21:02:32
% DurationCPUTime: 17.70s
% Computational Cost: add. (10575->497), mult. (23710->728), div. (0->0), fcn. (17168->10), ass. (0->240)
t198 = sin(qJ(3));
t202 = cos(qJ(3));
t203 = cos(qJ(2));
t276 = qJD(1) * t203;
t199 = sin(qJ(2));
t277 = qJD(1) * t199;
t164 = t198 * t276 + t202 * t277;
t165 = -t198 * t277 + t202 * t276;
t139 = -pkin(4) * t165 + t164 * pkin(5);
t197 = sin(qJ(4));
t201 = cos(qJ(4));
t319 = pkin(3) * qJD(2);
t263 = t202 * t319;
t117 = -t139 * t197 + t201 * t263;
t196 = sin(qJ(5));
t200 = cos(qJ(5));
t264 = t198 * t319;
t270 = qJD(5) * t200;
t269 = qJD(5) * t201;
t273 = qJD(4) * t200;
t347 = t196 * t269 + t197 * t273;
t340 = pkin(4) * t270 + t347 * pkin(5) - t200 * t117 + t196 * t264;
t279 = t200 * t201;
t125 = t164 * t279 + t165 * t196;
t271 = qJD(5) * t196;
t248 = t197 * t271;
t354 = (t248 - t125) * pkin(6);
t116 = t139 * t201 + t197 * t263;
t194 = t196 * pkin(4);
t163 = t194 + (-pkin(5) * t200 - pkin(6)) * t201;
t192 = pkin(7) * qJD(5) - qJD(6);
t330 = pkin(6) * t200;
t252 = -pkin(5) - t330;
t272 = qJD(4) * t201;
t353 = -t163 * t192 + t252 * t272 - t116 + t354;
t169 = t252 * t197;
t274 = qJD(4) * t197;
t191 = pkin(6) * t274;
t287 = t164 * t197;
t352 = -pkin(6) * t287 - t169 * t192 - t191 - t340;
t127 = -pkin(3) * t276 + t139;
t188 = pkin(3) * t198 - pkin(5);
t189 = pkin(3) * t202 + pkin(4);
t156 = t188 * t279 + t196 * t189;
t148 = -pkin(6) * t201 + t156;
t242 = t188 - t330;
t318 = pkin(3) * qJD(3);
t260 = t202 * t318;
t351 = t127 * t201 + t148 * t192 - t197 * t260 - t242 * t272 - t354;
t226 = t201 * t260;
t261 = t198 * t318;
t108 = -t347 * t188 + t189 * t270 - t196 * t261 + t200 * t226;
t161 = t242 * t197;
t350 = -t161 * t192 - t108 - t191 + (-pkin(6) * t164 - t127 * t200) * t197;
t349 = t196 * t272 + t197 * t270;
t187 = -pkin(3) * t199 - pkin(2);
t278 = qJD(1) * t187;
t120 = -pkin(4) * t164 - pkin(5) * t165 + t278;
t195 = qJD(2) + qJD(3);
t171 = -pkin(5) * t195 + t264;
t256 = qJD(2) * t318;
t224 = t202 * t256;
t166 = t198 * t203 + t199 * t202;
t143 = t195 * t166;
t133 = t143 * qJD(1);
t167 = t198 * t199 - t202 * t203;
t142 = t195 * t167;
t134 = t142 * qJD(1);
t262 = t203 * t319;
t71 = pkin(4) * t134 + pkin(5) * t133 - qJD(1) * t262;
t45 = t197 * t224 - t120 * t274 + (qJD(4) * t171 + t71) * t201;
t97 = t120 * t201 + t171 * t197;
t215 = t197 * t45 + t97 * t272;
t98 = -t197 * t120 + t171 * t201;
t348 = t98 * t274 - t215;
t147 = -t165 * t201 - t197 * t195;
t160 = qJD(4) + t164;
t112 = -t147 * t196 + t200 * t160;
t109 = t112 - t192;
t311 = t109 * Ifges(7,3);
t336 = Ifges(6,2) * t112 + t311;
t113 = t147 * t200 + t160 * t196;
t314 = Ifges(6,1) * t113;
t346 = mrSges(6,2) * t97 + t314;
t146 = t165 * t197 - t195 * t201;
t87 = t146 * qJD(4) + t133 * t201;
t40 = -qJD(5) * t113 - t134 * t200 - t196 * t87;
t345 = (Ifges(6,2) + Ifges(7,3)) * t40;
t193 = pkin(7) * qJ(5) - qJ(6);
t185 = sin(t193);
t186 = cos(t193);
t344 = t352 * t185 + t186 * t353;
t343 = t185 * t353 - t352 * t186;
t342 = t185 * t351 + t186 * t350;
t341 = t185 * t350 - t186 * t351;
t225 = t198 * t256;
t44 = -qJD(4) * t97 - t197 * t71 + t201 * t224;
t172 = pkin(4) * t195 + t263;
t75 = t172 * t196 + t200 * t98;
t22 = -qJD(5) * t75 - t196 * t44 - t200 * t225;
t281 = t196 * t201;
t74 = t172 * t200 - t196 * t98;
t339 = t22 * (pkin(4) * t200 + pkin(5) * t281) + (-pkin(4) * t271 + (-t196 * t274 + t200 * t269) * pkin(5) + t117 * t196 + t200 * t264) * t74;
t282 = t196 * t197;
t257 = t74 * t282;
t338 = -t127 * t257 + t22 * (-t188 * t281 + t189 * t200) + t74 * ((-t188 * t269 - t261) * t200 + (-qJD(5) * t189 + t188 * t274 - t226) * t196);
t140 = -pkin(4) * t166 + pkin(5) * t167 + t187;
t85 = pkin(4) * t142 + pkin(5) * t143 - t262;
t337 = t140 * t271 - t200 * t85;
t335 = -m(6) / 0.2e1;
t334 = pkin(6) * m(7);
t332 = t164 / 0.2e1;
t331 = -t165 / 0.2e1;
t329 = m(4) * t187;
t267 = -qJD(5) + t146;
t66 = t113 * t186 - t185 * t267;
t327 = Ifges(7,1) * t66;
t216 = t113 * t185 + t186 * t267;
t326 = Ifges(7,2) * t216;
t39 = qJD(5) * t112 - t134 * t196 + t200 * t87;
t325 = t39 * Ifges(6,1);
t324 = t40 * Ifges(7,3);
t323 = t75 * mrSges(6,2);
t322 = t87 * Ifges(5,1);
t317 = mrSges(4,3) * t164;
t316 = mrSges(4,3) * t165;
t315 = Ifges(5,1) * t147;
t313 = Ifges(4,4) * t165;
t310 = t116 * t97;
t309 = t125 * t97;
t308 = (t166 * t200 - t167 * t281) * t40;
t306 = t267 * Ifges(6,3);
t305 = t185 * t216;
t304 = t186 * t66;
t302 = t197 * t97;
t300 = t200 * t97;
t299 = t201 * t44;
t88 = t133 * t197 - t165 * t272 - t195 * t274;
t298 = t201 * t88;
t297 = t201 * t97;
t296 = Ifges(3,5) * qJD(2);
t295 = t112 * t185;
t294 = t112 * t186;
t293 = t134 * t166;
t292 = t140 * t200;
t291 = t267 * t200;
t290 = t146 * t197;
t289 = t146 * t200;
t286 = t185 * t192;
t285 = t185 * t197;
t284 = t186 * t192;
t283 = t186 * t197;
t280 = t197 * t200;
t275 = qJD(4) * t140;
t266 = pkin(2) * t203 * mrSges(3,1);
t265 = t197 * t323;
t259 = Ifges(5,2) * t290;
t255 = t97 * t275;
t245 = t296 / 0.2e1;
t61 = pkin(6) * t267 + t75;
t238 = -pkin(6) * t39 + t192 * t61 - t45;
t21 = t172 * t270 + t200 * t44 + (-qJD(5) * t98 - t225) * t196;
t60 = pkin(6) * t113 + t97;
t239 = pkin(6) * t88 - t192 * t60 - t21;
t2 = t185 * t238 + t186 * t239;
t217 = t185 * t61 - t186 * t60;
t244 = -t192 * t217 + t2;
t25 = t185 * t60 + t186 * t61;
t3 = t185 * t239 - t186 * t238;
t243 = t192 * t25 + t3;
t233 = t192 * t267 - t39;
t234 = t113 * t192 - t88;
t10 = t185 * t234 + t186 * t233;
t241 = t192 * t216 - t10;
t11 = t185 * t233 - t186 * t234;
t240 = t192 * t66 + t11;
t220 = -pkin(6) * t167 - t292;
t136 = t166 * t196 + t167 * t279;
t77 = pkin(6) * t136 + t140 * t201;
t237 = -t192 * t77 - t220 * t272 - (-pkin(6) * t143 + t337) * t197;
t236 = t85 * t257 + (t22 * t282 + t349 * t74) * t140;
t105 = t220 * t197;
t212 = -qJD(5) * t166 - t143 * t201 + t167 * t274;
t222 = -t167 * t269 - t142;
t53 = t196 * t222 - t200 * t212;
t235 = -pkin(6) * t53 + t105 * t192 + t140 * t274 - t201 * t85;
t228 = qJD(4) * t188 - t127;
t227 = t192 * t200 - qJD(4);
t223 = -t167 * t192 * t197 - t53;
t221 = t140 * t45 + t85 * t97;
t219 = -Ifges(5,2) * t146 - t306;
t218 = t185 * t25 - t186 * t217;
t154 = -t185 * t201 + t186 * t280;
t78 = -t125 * t185 + t164 * t283;
t211 = -t248 + (-t192 + t273) * t201;
t81 = t186 * t211 - t227 * t285;
t214 = -t154 * t3 + t217 * t81 + t25 * t78;
t213 = -t136 * t192 + t143 * t197 + t167 * t272;
t210 = t185 * t243 + t186 * t244;
t209 = -t185 * t244 + t186 * t243;
t208 = t160 * t306 - t200 * t325 - t322;
t207 = t88 * Ifges(5,2) - t160 * t315 - t273 * t314;
t206 = -t299 + t348;
t128 = Ifges(4,2) * t164 + Ifges(4,6) * t195 - t313;
t158 = Ifges(4,4) * t164;
t129 = -Ifges(4,1) * t165 + Ifges(4,5) * t195 + t158;
t153 = t185 * t280 + t186 * t201;
t79 = -t125 * t186 - t164 * t285;
t82 = t185 * t211 + t227 * t283;
t205 = t346 * t248 + (-t164 * t297 + t98 * t287 + t348) * mrSges(5,3) - (-mrSges(4,1) * t165 + mrSges(4,2) * t164) * t278 - (Ifges(4,2) * t165 + t129 + t158) * t164 / 0.2e1 + (t2 * t153 - t217 * t79 - t25 * t82) * mrSges(7,3) + t164 * t265 + t164 * t259 + t263 * t317 + Ifges(6,3) * t298 - t125 * t314 + t165 * (Ifges(4,1) * t164 + t313) / 0.2e1 + t128 * t331 + (-t81 + t79) * t327 - t195 * (Ifges(4,5) * t164 + Ifges(4,6) * t165) / 0.2e1 + (t78 - t82) * t326 + t11 * Ifges(7,2) * t153 + t10 * Ifges(7,1) * t154 - t160 * Ifges(5,3) * t165 + (t265 + t259) * qJD(4) + t282 * t345 + Ifges(4,5) * t133 - Ifges(4,6) * t134 + (t164 * t281 - t165 * t200 + t349) * t336;
t204 = qJD(1) ^ 2;
t175 = -Ifges(3,1) * t276 + t296;
t168 = -pkin(5) * t279 + t194;
t150 = mrSges(4,1) * t195 + t316;
t149 = -mrSges(4,2) * t195 + t317;
t138 = -mrSges(4,1) * t164 - mrSges(4,2) * t165;
t122 = -t163 * t186 - t169 * t185;
t121 = -t163 * t185 + t169 * t186;
t107 = -t148 * t186 - t161 * t185;
t106 = -t148 * t185 + t161 * t186;
t93 = -t136 * t186 - t167 * t285;
t92 = -t136 * t185 + t167 * t283;
t84 = -t147 * t185 - t186 * t289;
t83 = t147 * t186 - t185 * t289;
t67 = pkin(6) * t289 + t98;
t64 = -pkin(6) * t147 - t300;
t54 = t196 * t212 + t200 * t222;
t50 = -t105 * t186 - t185 * t77;
t49 = -t105 * t185 + t186 * t77;
t48 = -pkin(7) * t66 - t295;
t47 = pkin(7) * t216 - t294;
t28 = -t185 * t67 - t186 * t64;
t27 = -t185 * t64 + t186 * t67;
t19 = t185 * t223 + t186 * t213;
t18 = -t185 * t213 + t186 * t223;
t7 = pkin(6) * t294 - pkin(7) * t25 - t185 * t74;
t6 = -pkin(6) * t295 + pkin(7) * t217 - t186 * t74;
t5 = t185 * t237 - t186 * t235;
t4 = t185 * t235 + t186 * t237;
t1 = [t195 * (Ifges(4,5) * t143 - Ifges(4,6) * t142) / 0.2e1 - t142 * t128 / 0.2e1 + t143 * t129 / 0.2e1 + (mrSges(4,1) * t134 + mrSges(4,2) * t133 + qJD(1) * (mrSges(4,1) * t142 + mrSges(4,2) * t143)) * t187 + (t109 * t54 + t308) * Ifges(7,3) + (-t142 * t160 - t293) * Ifges(5,3) + (t11 * t92 - t19 * t216) * Ifges(7,2) + (t112 * t54 + t308) * Ifges(6,2) + (-t142 * t332 - t293) * Ifges(4,2) + (t10 * t93 - t18 * t66) * Ifges(7,1) + (t113 * t53 + t136 * t39) * Ifges(6,1) + (t133 * t167 + t143 * t331) * Ifges(4,1) + (t136 * t45 + t53 * t97) * mrSges(6,2) + m(7) * (t2 * t50 - t217 * t5 - t25 * t4 + t3 * t49 + t236) + m(6) * t236 + (-t10 * t49 + t11 * t50 + t18 * t217 - t19 * t25 + t2 * t92 - t216 * t4 - t3 * t93 + t5 * t66) * mrSges(7,3) + (t133 * t166 - t134 * t167 - t142 * t331 + t143 * t332) * Ifges(4,4) + (0.2e1 * qJD(1) * t266 + (t175 / 0.2e1 + t245 + (-0.3e1 / 0.2e1 * Ifges(3,1) + (2 * Ifges(3,2))) * t276) * t199 + ((-t138 + (mrSges(4,1) * t166 - mrSges(4,2) * t167 - 0.2e1 * t329) * qJD(1)) * t203 + (-t198 * t142 - t202 * t143 + (t166 * t202 + t167 * t198) * qJD(3)) * mrSges(4,3)) * pkin(3)) * qJD(2) + (t143 * t315 + m(5) * (-t275 * t98 + t221) + m(6) * (-t140 * t273 * t75 + t221) + (t219 * qJD(4) + t322) * t167 + (t85 * t113 + t140 * t39 + (-t140 * t291 - t167 * t75) * qJD(4)) * mrSges(6,2) + (t140 * t87 + t97 * t143 + t85 * t147 + t45 * t167 + (-t140 * t146 - t167 * t98) * qJD(4)) * mrSges(5,3)) * t201 + (t219 * t143 + m(5) * (-t140 * t44 - t85 * t98 - t255) + m(6) * (-t21 * t292 + t337 * t75 - t255) + (-qJD(4) * t315 + (Ifges(5,2) + Ifges(6,3)) * t88) * t167 + (t140 * t88 - t98 * t143 - t85 * t146 - t44 * t167 + (-t140 * t147 - t167 * t97) * qJD(4)) * mrSges(5,3) + (-t85 * t291 - t75 * t143 - t21 * t167 + (-qJD(4) * t113 + t200 * t88 + t267 * t271) * t140) * mrSges(6,2)) * t197; ((qJD(1) * t138 + t204 * t329) * t203 + (-t202 * t133 + (-qJD(2) * t165 - t134) * t198) * mrSges(4,3) + ((-qJD(2) * mrSges(4,1) + m(5) * (-qJD(2) * t189 - t172) - t150) * t198 + (-qJD(2) * mrSges(4,2) + t149 + (t146 * t201 + t147 * t197) * mrSges(5,3) + m(5) * (t201 * t98 + t302)) * t202) * qJD(3)) * pkin(3) + ((-t146 * t228 + t188 * t87) * mrSges(5,3) + t208) * t197 + (t108 * t267 - t309 - t156 * t88 + (-t113 * t127 - t21 + (t113 * t188 - t300) * qJD(4)) * t201 + (t113 * t260 + t188 * t39 + (t127 * t267 - t45) * t200) * t197) * mrSges(6,2) + 0.2e1 * (-m(5) * (-t197 * t98 + t297) / 0.2e1 + (-t280 * t75 + t297) * t335) * t127 - m(5) * t206 * t188 + t205 + ((t147 * t228 - t188 * t88 - t44) * mrSges(5,3) + t207) * t201 + (-t10 * t106 + t107 * t11 - t216 * t342 + t341 * t66 + t214) * mrSges(7,3) - t204 * t266 + ((Ifges(3,1) / 0.2e1 - Ifges(3,2)) * t204 * t203 + (-t175 / 0.2e1 + t245) * qJD(1)) * t199 + (t106 * t3 + t107 * t2 - t217 * t341 - t25 * t342 + t338) * m(7) + (t108 * t75 + t156 * t21 + t188 * t215 + t302 * t260 + t338) * m(6); (-t116 * t147 - t117 * t146 - t299 + (-t87 * t197 + t298 + (-t147 * t201 + t290) * qJD(4)) * pkin(5)) * mrSges(5,3) + 0.2e1 * t215 * t335 * pkin(5) + t205 + t208 * t197 + ((-mrSges(4,2) * qJD(3) - t149) * t202 + (-mrSges(4,1) * qJD(3) + t150 - t316) * t198) * t319 + t207 * t201 + (-t10 * t121 + t11 * t122 + t216 * t343 + t344 * t66 + t214) * mrSges(7,3) + (-t113 * t116 - t309 - t168 * t88 + (-pkin(5) * t39 - t200 * t45) * t197 + t340 * t267 + (-t21 + (-pkin(5) * t113 - t300) * qJD(4)) * t201) * mrSges(6,2) + (t121 * t3 + t122 * t2 - t217 * t344 + t343 * t25 + t339) * m(7) + (t168 * t21 + t340 * t75 - t310 + t339) * m(6) + (-pkin(4) * t225 + pkin(5) * t206 - t117 * t98 + t172 * t264 - t310) * m(5); t84 * t327 + t83 * t326 - t134 * Ifges(5,3) - m(7) * (-t217 * t27 - t25 * t28) + (t323 + t306 + (-Ifges(5,1) + Ifges(5,2)) * t146) * t147 + (t45 * mrSges(6,2) + t325 + (-t11 * t185 + t216 * t284) * Ifges(7,2) + (-t10 * t186 - t286 * t66) * Ifges(7,1) + 0.2e1 * (t335 - m(7) / 0.2e1) * t97 * t74 + t336 * t267) * t196 + (-t146 * t314 + t345 + (m(6) * t75 + (t267 - t146) * mrSges(6,2)) * t97 + (Ifges(7,1) * t304 + Ifges(7,2) * t305 + t346) * qJD(5)) * t200 + (t209 * t196 + t218 * t270) * t334 + (-t217 * t84 + t25 * t83 - t27 * t66 + t28 * t216 + ((t304 + t305) * pkin(6) + t218) * t270 + ((-t185 * t240 + t186 * t241) * pkin(6) + t209) * t196) * mrSges(7,3) + (-m(6) * t97 - mrSges(6,2) * t113) * t98; -pkin(7) * t324 + t88 * Ifges(6,3) - m(7) * (-t217 * t7 - t25 * t6 - t75 * t74) + (t311 + (-Ifges(6,1) + Ifges(6,2)) * t112) * t113 + (t11 * t186 - (-t48 - t286) * t216) * Ifges(7,2) + (-t10 * t185 + (t47 + t284) * t66) * Ifges(7,1) + (-t112 * t97 - t267 * t74 - t21) * mrSges(6,2) + t210 * t334 + (-t217 * t47 + t25 * t48 + t6 * t216 - t7 * t66 + (t185 * t241 + t186 * t240) * pkin(6) + t210) * mrSges(7,3); t324 - (Ifges(7,1) - Ifges(7,2)) * t66 * t216;];
tauc = t1(:);
