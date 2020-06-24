% Calculate vector of inverse dynamics joint torques for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = CloosQRC350DE_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:02:09
% EndTime: 2020-06-23 21:02:37
% DurationCPUTime: 22.72s
% Computational Cost: add. (12180->553), mult. (25124->776), div. (0->0), fcn. (18210->10), ass. (0->282)
t226 = cos(qJ(3));
t395 = sin(qJ(2));
t279 = qJD(1) * t395;
t227 = cos(qJ(2));
t394 = sin(qJ(3));
t294 = t227 * t394;
t176 = qJD(1) * t294 + t226 * t279;
t255 = t395 * t394;
t319 = qJD(1) * t227;
t177 = -qJD(1) * t255 + t226 * t319;
t222 = sin(qJ(5));
t224 = cos(qJ(5));
t225 = cos(qJ(4));
t322 = t224 * t225;
t134 = t176 * t322 + t177 * t222;
t223 = sin(qJ(4));
t314 = qJD(5) * t222;
t286 = t223 * t314;
t426 = (t286 - t134) * pkin(6);
t306 = pkin(3) * t394;
t209 = t306 - pkin(5);
t393 = pkin(3) * t226;
t210 = pkin(4) + t393;
t163 = t209 * t322 + t222 * t210;
t153 = -pkin(6) * t225 + t163;
t217 = pkin(7) * qJD(5) - qJD(6);
t370 = pkin(3) * qJD(3);
t303 = t226 * t370;
t262 = t223 * t303;
t392 = pkin(6) * t224;
t276 = t209 - t392;
t315 = qJD(4) * t225;
t144 = -pkin(4) * t177 + t176 * pkin(5);
t305 = pkin(3) * t319;
t136 = t144 - t305;
t338 = t136 * t225;
t425 = t153 * t217 - t276 * t315 - t262 + t338 - t426;
t259 = qJD(3) * t306;
t261 = t225 * t303;
t313 = qJD(5) * t224;
t312 = qJD(5) * t225;
t316 = qJD(4) * t224;
t415 = t222 * t312 + t223 * t316;
t110 = -t415 * t209 + t210 * t313 - t222 * t259 + t224 * t261;
t172 = t276 * t223;
t317 = qJD(4) * t223;
t214 = pkin(6) * t317;
t424 = -t172 * t217 - t110 - t214 + (-pkin(6) * t176 - t136 * t224) * t223;
t371 = pkin(3) * qJD(2);
t304 = t226 * t371;
t126 = -t144 * t223 + t225 * t304;
t260 = qJD(2) * t306;
t408 = -pkin(4) * t313 - t415 * pkin(5) + t224 * t126 - t222 * t260;
t125 = t144 * t225 + t223 * t304;
t219 = t222 * pkin(4);
t174 = t219 + (-pkin(5) * t224 - pkin(6)) * t225;
t295 = -pkin(5) - t392;
t423 = -t174 * t217 + t295 * t315 - t125 + t426;
t187 = t295 * t223;
t332 = t176 * t223;
t422 = -pkin(6) * t332 - t187 * t217 - t214 + t408;
t208 = -pkin(3) * t395 - pkin(2);
t320 = qJD(1) * t208;
t129 = -pkin(4) * t176 - pkin(5) * t177 + t320;
t221 = qJD(2) + qJD(3);
t189 = -t221 * pkin(5) + t260;
t100 = -t223 * t129 + t189 * t225;
t181 = qJD(2) * t303 + qJDD(2) * t306;
t220 = qJDD(2) + qJDD(3);
t167 = -pkin(5) * t220 + t181;
t278 = qJD(2) * t395;
t185 = qJD(1) * t278 - t227 * qJDD(1);
t272 = t395 * qJDD(1);
t318 = qJD(2) * t227;
t277 = qJD(1) * t318;
t186 = t272 + t277;
t182 = t226 * t395 + t294;
t234 = t182 * qJD(3);
t121 = qJD(1) * t234 + t226 * t185 + t186 * t394;
t122 = t177 * qJD(3) - t185 * t394 + t226 * t186;
t175 = -pkin(3) * t277 + t208 * qJDD(1);
t61 = -pkin(4) * t122 + pkin(5) * t121 + t175;
t31 = t167 * t223 - t129 * t317 + (qJD(4) * t189 + t61) * t225;
t350 = t223 * t31;
t99 = t129 * t225 + t189 * t223;
t239 = t99 * t315 + t350;
t421 = t100 * t317 - t239;
t420 = t222 * t315 + t223 * t313;
t218 = pkin(7) * qJ(5) - qJ(6);
t205 = sin(t218);
t206 = cos(t218);
t419 = t425 * t205 + t424 * t206;
t418 = t424 * t205 - t425 * t206;
t180 = -qJD(2) * t259 + qJDD(2) * t393;
t166 = t220 * pkin(4) + t180;
t30 = -qJD(4) * t99 + t167 * t225 - t223 * t61;
t190 = pkin(4) * t221 + t304;
t80 = t100 * t224 + t190 * t222;
t18 = -qJD(5) * t80 + t166 * t224 - t222 * t30;
t325 = t222 * t223;
t79 = -t100 * t222 + t190 * t224;
t298 = t79 * t325;
t324 = t222 * t225;
t407 = -t209 * t317 + t261;
t417 = -t136 * t298 + t18 * (-t209 * t324 + t210 * t224) + t79 * ((-t209 * t312 - t259) * t224 + (-qJD(5) * t210 - t407) * t222);
t399 = pkin(6) * m(7);
t216 = mrSges(6,2) + mrSges(7,3) + t399;
t374 = m(7) + m(5) + m(6);
t165 = pkin(4) * t374 - t216 * t324 + mrSges(4,1);
t168 = pkin(5) * t374 + t216 * t224 + mrSges(4,2) + mrSges(5,3);
t416 = t165 * t226 - t168 * t394;
t171 = qJD(4) + t176;
t242 = t177 * t225 + t223 * t221;
t116 = t171 * t222 - t224 * t242;
t359 = t116 * Ifges(6,1);
t375 = t99 * mrSges(6,2);
t414 = t359 + t375;
t115 = t224 * t171 + t222 * t242;
t111 = t115 - t217;
t360 = t111 * Ifges(7,3);
t404 = Ifges(6,2) * t115 + t360;
t120 = qJDD(4) + t122;
t151 = t177 * t223 - t221 * t225;
t73 = qJD(4) * t151 + t121 * t225 - t220 * t223;
t37 = -qJD(5) * t116 + t224 * t120 - t222 * t73;
t33 = -pkin(7) * qJDD(5) + qJDD(6) + t37;
t384 = t33 * Ifges(7,3);
t413 = t37 * Ifges(6,2) + t384;
t412 = t422 * t205 + t423 * t206;
t411 = t423 * t205 - t422 * t206;
t410 = t18 * (pkin(4) * t224 + pkin(5) * t324) + (-pkin(4) * t314 + (-t222 * t317 + t224 * t312) * pkin(5) + t222 * t126 + t224 * t260) * t79;
t406 = t209 * t315 + t262;
t400 = pkin(3) * m(4);
t405 = -pkin(2) * mrSges(3,1) + t208 * t400;
t333 = t242 * t225;
t336 = t151 * t223;
t243 = t333 + t336;
t323 = t223 * t224;
t160 = -t205 * t225 + t206 * t323;
t62 = pkin(6) * t116 + t99;
t310 = -qJD(5) + t151;
t63 = pkin(6) * t310 + t80;
t247 = t205 * t63 - t206 * t62;
t25 = t205 * t62 + t206 * t63;
t36 = qJD(5) * t115 + t120 * t222 + t224 * t73;
t280 = -pkin(6) * t36 + t217 * t63 - t31;
t17 = -t100 * t314 + t222 * t166 + t190 * t313 + t224 * t30;
t74 = qJD(4) * t242 - t121 * t223 - t220 * t225;
t71 = qJDD(5) - t74;
t281 = pkin(6) * t71 - t217 * t62 - t17;
t3 = t205 * t281 - t206 * t280;
t327 = t206 * t223;
t83 = -t134 * t205 + t176 * t327;
t233 = -t286 + (-t217 + t316) * t225;
t263 = t217 * t224 - qJD(4);
t329 = t205 * t223;
t85 = t206 * t233 - t263 * t329;
t403 = -t3 * t160 + t247 * t85 + t25 * t83;
t268 = t217 * t310 - t36;
t269 = t116 * t217 - t71;
t10 = t205 * t269 + t206 * t268;
t11 = t205 * t268 - t206 * t269;
t353 = t177 * Ifges(4,4);
t137 = t176 * Ifges(4,2) + t221 * Ifges(4,6) - t353;
t169 = Ifges(4,4) * t176;
t138 = -t177 * Ifges(4,1) + t221 * Ifges(4,5) + t169;
t159 = t205 * t323 + t206 * t225;
t2 = t205 * t280 + t206 * t281;
t300 = Ifges(5,2) * t336;
t376 = t80 * mrSges(6,2);
t308 = t223 * t376;
t344 = t225 * t99;
t354 = t171 * Ifges(5,3);
t366 = mrSges(5,3) * t100;
t368 = mrSges(4,3) * t176;
t379 = t71 * Ifges(6,3);
t245 = t116 * t205 + t206 * t310;
t390 = Ifges(7,2) * t245;
t69 = t116 * t206 - t205 * t310;
t391 = Ifges(7,1) * t69;
t396 = -t177 / 0.2e1;
t84 = -t134 * t206 - t176 * t329;
t86 = t205 * t233 + t263 * t327;
t402 = t413 * t325 + t414 * t286 + (-t85 + t84) * t391 + (-t86 + t83) * t390 + (t308 + t300) * qJD(4) + (t176 * t324 - t177 * t224 + t420) * t404 + (-t176 * t344 + t421) * mrSges(5,3) - (Ifges(4,2) * t177 + t138 + t169) * t176 / 0.2e1 + (t2 * t159 - t247 * t84 - t25 * t86) * mrSges(7,3) - (-mrSges(4,1) * t177 + mrSges(4,2) * t176) * t320 + t332 * t366 + t304 * t368 + t225 * t379 + t137 * t396 + t176 * t308 + t176 * t300 + Ifges(4,3) * t220 - t221 * (Ifges(4,5) * t176 + Ifges(4,6) * t177) / 0.2e1 + t180 * mrSges(4,1) - t181 * mrSges(4,2) + Ifges(4,5) * t121 + Ifges(4,6) * t122 + t11 * Ifges(7,2) * t159 + t10 * Ifges(7,1) * t160 + t177 * (Ifges(4,1) * t176 + t353) / 0.2e1 - t177 * t354 - t134 * t359;
t401 = -m(6) / 0.2e1;
t389 = g(3) * t216;
t388 = t17 * mrSges(6,2);
t383 = t36 * Ifges(6,1);
t378 = t73 * Ifges(5,1);
t377 = t74 * Ifges(5,2);
t369 = mrSges(6,2) * t116;
t367 = mrSges(4,3) * t177;
t362 = Ifges(5,1) * t242;
t358 = t120 * Ifges(5,3);
t357 = t134 * t99;
t355 = t310 * Ifges(6,3);
t352 = t205 * t245;
t351 = t206 * t69;
t183 = t227 * t226 - t255;
t147 = t221 * t183;
t148 = qJD(2) * t182 + t234;
t89 = -pkin(3) * t318 - pkin(4) * t147 + pkin(5) * t148;
t349 = t224 * t89;
t348 = t224 * t99;
t347 = t225 * t30;
t346 = t225 * t74;
t345 = t225 * t89;
t343 = t73 * t223;
t342 = t99 * t125;
t341 = Ifges(3,5) * qJD(2);
t340 = t115 * t205;
t339 = t115 * t206;
t337 = t148 * t223;
t335 = t151 * t224;
t331 = t183 * t223;
t330 = t205 * t217;
t328 = t206 * t217;
t326 = t209 * t223;
t307 = t224 * t383;
t302 = mrSges(6,2) * t323;
t301 = t223 * t355;
t299 = Ifges(5,1) * t333;
t297 = Ifges(3,2) * t395;
t287 = t224 * t315;
t283 = -t217 * t247 + t2;
t282 = t217 * t25 + t3;
t275 = t217 * t245 - t10;
t274 = t217 * t69 + t11;
t145 = -t182 * pkin(4) - t183 * pkin(5) + t208;
t250 = pkin(6) * t183 - t145 * t224;
t141 = t182 * t222 - t183 * t322;
t82 = pkin(6) * t141 + t145 * t225;
t273 = -t217 * t82 - t250 * t315 - (-pkin(6) * t148 + t145 * t314 - t349) * t223;
t271 = t89 * t298 + (t18 * t325 + t420 * t79) * t145;
t107 = t250 * t223;
t235 = qJD(5) * t182 + t148 * t225 + t183 * t317;
t251 = t183 * t312 + t147;
t54 = t222 * t251 + t224 * t235;
t270 = -pkin(6) * t54 + t107 * t217 + t145 * t317 - t345;
t258 = t394 * t367;
t257 = t316 * t359;
t252 = t217 * t331 - t54;
t249 = -Ifges(5,2) * t151 - t355;
t248 = t205 * t25 - t206 * t247;
t246 = -t100 * t223 + t344;
t240 = -t323 * t80 + t344;
t238 = t225 * t31 - t317 * t99;
t237 = m(5) * t246;
t236 = t141 * t217 + t183 * t315 - t337;
t232 = t205 * t282 + t206 * t283;
t231 = -t205 * t283 + t206 * t282;
t193 = -Ifges(3,1) * t319 + t341;
t184 = -pkin(5) * t322 + t219;
t155 = mrSges(4,1) * t221 + t367;
t154 = -mrSges(4,2) * t221 + t368;
t143 = -mrSges(4,1) * t176 - mrSges(4,2) * t177;
t140 = t182 * t224 + t183 * t324;
t131 = -t174 * t206 - t187 * t205;
t130 = -t174 * t205 + t187 * t206;
t123 = (t165 * t394 + t168 * t226) * t227;
t109 = -t153 * t206 - t172 * t205;
t108 = -t153 * t205 + t172 * t206;
t95 = -t141 * t206 + t183 * t329;
t94 = -t141 * t205 - t183 * t327;
t88 = t205 * t242 - t206 * t335;
t87 = -t205 * t335 - t206 * t242;
t70 = pkin(6) * t335 + t100;
t67 = pkin(6) * t242 - t348;
t55 = -t222 * t235 + t224 * t251;
t51 = -t107 * t206 - t205 * t82;
t50 = -t107 * t205 + t206 * t82;
t49 = -pkin(7) * t69 - t340;
t48 = pkin(7) * t245 - t339;
t35 = -t205 * t70 - t206 * t67;
t34 = -t205 * t67 + t206 * t70;
t22 = t205 * t252 - t206 * t236;
t21 = t205 * t236 + t206 * t252;
t13 = pkin(6) * t339 - pkin(7) * t25 - t205 * t79;
t12 = -pkin(6) * t340 + pkin(7) * t247 - t206 * t79;
t5 = t205 * t273 - t206 * t270;
t4 = t205 * t270 + t206 * t273;
t1 = [(-Ifges(3,1) * t185 - Ifges(3,5) * qJDD(2) + (-pkin(3) * t143 + (t297 - Ifges(3,1) * t395 / 0.2e1 - t405) * qJD(1)) * qJD(2)) * t227 + (t116 * t345 + t31 * t141 + t99 * t54 + (-t148 * t80 - t310 * t349) * t223) * mrSges(6,2) + (m(5) * (-t100 * t315 - t223 * t30 + t238) + (-t223 * t74 + t225 * t73 + (-t151 * t225 + t223 * t242) * qJD(4)) * mrSges(5,3) + ((-t310 * t316 + t36) * t225 + (-qJD(4) * t116 + t224 * t71 + t310 * t314) * t223) * mrSges(6,2)) * t145 + m(7) * (t2 * t51 - t247 * t5 - t25 * t4 + t3 * t50 + t271) + (-t10 * t50 + t11 * t51 + t2 * t94 + t21 * t247 - t22 * t25 - t245 * t4 - t3 * t95 + t5 * t69) * mrSges(7,3) + (t11 * t94 - t22 * t245) * Ifges(7,2) + (-mrSges(4,2) * t175 + mrSges(4,3) * t180 - Ifges(4,1) * t121 - Ifges(4,4) * t122 - Ifges(4,5) * t220 + (-mrSges(5,3) * t31 - t378) * t225 + (t30 * mrSges(5,3) + t377 - t379 + t388) * t223 + ((mrSges(5,3) * t99 - t362) * t223 + (-t249 + t366 + t376) * t225) * qJD(4)) * t183 + (t271 + t240 * t89 + (-t17 * t323 + t238 + (t286 - t287) * t80) * t145) * m(6) + t237 * t89 + (t148 * t246 - t243 * t89) * mrSges(5,3) + (t10 * t95 - t21 * t69) * Ifges(7,1) + (t116 * t54 + t141 * t36) * Ifges(6,1) + (t115 * t55 + t140 * t37) * Ifges(6,2) + t147 * t354 + (m(3) * qJDD(1) * pkin(2) + (t272 + t186) * mrSges(3,1)) * pkin(2) + t249 * t337 + (Ifges(4,1) * t148 + Ifges(4,4) * t147) * t396 + (-mrSges(4,1) * t122 + mrSges(4,2) * t121 + qJD(1) * (-mrSges(4,1) * t147 + mrSges(4,2) * t148) + m(4) * t175) * t208 + (t111 * t55 + t140 * t33) * Ifges(7,3) - t148 * t299 + Ifges(2,3) * qJDD(1) + t221 * (Ifges(4,5) * t148 + Ifges(4,6) * t147) / 0.2e1 + t176 * (Ifges(4,4) * t148 + Ifges(4,2) * t147) / 0.2e1 + t147 * t137 / 0.2e1 + t148 * t138 / 0.2e1 + (t395 * t193 / 0.2e1 + Ifges(3,5) * t278 / 0.2e1 + (t147 * t394 - t148 * t226) * mrSges(4,3) * pkin(3)) * qJD(2) + t186 * t297 + (-mrSges(4,1) * t175 + mrSges(4,3) * t181 + Ifges(4,4) * t121 + Ifges(4,2) * t122 + Ifges(4,6) * t220 + t358) * t182; (-t108 * t10 + t109 * t11 - t419 * t245 + t418 * t69 + t403) * mrSges(7,3) + (t108 * t3 + t109 * t2 - t418 * t247 - t419 * t25 + t417) * m(7) + (t110 * t80 - t136 * t240 + t163 * t17 + t209 * t239 + t262 * t99 + t417) * m(6) + g(3) * ((mrSges(3,1) + (m(4) + t374) * pkin(3) + t416) * t395 + t123) + (-t338 + t406) * t369 + t402 + (t302 * t310 - t237) * t136 + (t110 * t310 - t163 * t71 + t36 * t326 - t357) * mrSges(6,2) + (t243 * t136 + t407 * t151 + t209 * t346 - t242 * t406 + t73 * t326 - t347) * mrSges(5,3) + (t299 + t301) * t176 - (t341 + t193) * t279 / 0.2e1 + (t405 + (Ifges(3,1) / 0.2e1 - Ifges(3,2)) * t395) * t227 * qJD(1) ^ 2 - t225 * t257 + (t180 * t226 + t181 * t394) * t400 + (mrSges(4,1) * t220 - mrSges(4,3) * t121) * t393 + t154 * t303 + t143 * t305 + Ifges(3,5) * t185 + Ifges(3,3) * qJDD(2) + m(5) * (t166 * t210 + (-t394 * t190 + (t100 * t225 + t223 * t99) * t226) * t370 + (qJD(4) * t246 + t347 + t350) * t209) + qJD(4) * t299 + qJD(4) * t301 - t31 * t302 + (-mrSges(4,2) * t220 + mrSges(4,3) * t122) * t306 - t223 * t307 - Ifges(5,1) * t343 - Ifges(5,2) * t346 - t155 * t259 - t258 * t371 - t287 * t375 - t225 * t388; (-t116 * t125 - t357 - t184 * t71 + (-pkin(5) * t36 - t224 * t31) * t223 - t408 * t310 + (-t17 + (-pkin(5) * t116 - t348) * qJD(4)) * t225) * mrSges(6,2) + (t125 * t242 - t126 * t151 - t347 + (qJD(4) * t243 - t343 - t346) * pkin(5)) * mrSges(5,3) + (-t130 * t10 + t131 * t11 + t245 * t411 + t412 * t69 + t403) * mrSges(7,3) + (t171 * t355 - t307 - t378) * t223 + (t171 * t362 - t257 - t377) * t225 + 0.2e1 * (m(5) * (-t347 + t421) / 0.2e1 + t239 * t401) * pkin(5) - m(5) * (t100 * t126 - t190 * t260 + t342) + m(5) * pkin(4) * t166 + (-t154 * t226 + t155 * t394 - t258) * t371 + g(3) * (t395 * t416 + t123) + (t130 * t3 + t131 * t2 - t247 * t412 + t25 * t411 + t410) * m(7) + (t17 * t184 - t408 * t80 - t342 + t410) * m(6) + t402; t88 * t391 + t87 * t390 + t358 - m(7) * (-t247 * t34 - t25 * t35) - (t376 + t355 + (-Ifges(5,1) + Ifges(5,2)) * t151) * t242 + (-t151 * t359 + (m(6) * t80 + (-t151 + t310) * mrSges(6,2)) * t99 + (Ifges(7,1) * t351 + Ifges(7,2) * t352 + t248 * t399 + t414) * qJD(5) + t413) * t224 + (-t331 * t389 + t31 * mrSges(6,2) + t383 + (-t11 * t205 + t245 * t328) * Ifges(7,2) + (-t10 * t206 - t330 * t69) * Ifges(7,1) + 0.2e1 * (-m(7) / 0.2e1 + t401) * t99 * t79 + t231 * t399 + t404 * t310) * t222 + (-t247 * t88 + t25 * t87 - t34 * t69 + t35 * t245 + ((t351 + t352) * pkin(6) + t248) * t313 + ((-t205 * t274 + t206 * t275) * pkin(6) + t231) * t222) * mrSges(7,3) + (-m(6) * t99 - t369) * t100; t379 - pkin(7) * t384 + ((-t222 * t394 + t226 * t322) * t227 - t395 * (t222 * t226 + t322 * t394)) * t389 - m(7) * (-t12 * t25 - t13 * t247 - t80 * t79) + (t360 + (-Ifges(6,1) + Ifges(6,2)) * t115) * t116 + (t11 * t206 - (-t49 - t330) * t245) * Ifges(7,2) + (-t10 * t205 + (t48 + t328) * t69) * Ifges(7,1) + (-t115 * t99 - t310 * t79 - t17) * mrSges(6,2) + t232 * t399 + (t12 * t245 - t13 * t69 - t247 * t48 + t25 * t49 + (t205 * t275 + t206 * t274) * pkin(6) + t232) * mrSges(7,3); t384 - (Ifges(7,1) - Ifges(7,2)) * t69 * t245;];
tau = t1;
