% Calculate vector of inverse dynamics joint torques for
% CloosQRC350OL
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = CloosQRC350OL_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:56:05
% EndTime: 2020-06-23 21:56:34
% DurationCPUTime: 20.15s
% Computational Cost: add. (11205->539), mult. (24318->767), div. (0->0), fcn. (17935->12), ass. (0->267)
t220 = sin(qJ(4));
t215 = qJD(2) + qJD(3);
t225 = cos(qJ(4));
t226 = cos(qJ(3));
t227 = cos(qJ(2));
t307 = qJD(1) * t227;
t221 = sin(qJ(3));
t222 = sin(qJ(2));
t312 = t221 * t222;
t386 = qJD(1) * t312 - t226 * t307;
t250 = t215 * t225 - t220 * t386;
t244 = qJD(5) + t250;
t239 = t220 * t244;
t396 = Ifges(6,3) * t239;
t182 = -t221 * t227 - t222 * t226;
t173 = t182 * qJD(1);
t140 = -pkin(4) * t386 + t173 * pkin(5);
t353 = pkin(3) * qJD(2);
t290 = t226 * t353;
t121 = -t140 * t220 + t225 * t290;
t219 = sin(qJ(5));
t224 = cos(qJ(5));
t291 = t221 * t353;
t302 = qJD(5) * t224;
t301 = qJD(5) * t225;
t305 = qJD(4) * t224;
t385 = t219 * t301 + t220 * t305;
t330 = -pkin(4) * t302 - pkin(5) * t385 + t224 * t121 - t219 * t291;
t309 = t224 * t225;
t128 = t173 * t309 + t219 * t386;
t303 = qJD(5) * t219;
t276 = t220 * t303;
t395 = (-t128 + t276) * pkin(6);
t213 = t222 * pkin(3);
t356 = pkin(2) + t213;
t193 = t356 * qJD(1);
t234 = -pkin(4) * t173 - pkin(5) * t386 + t193;
t122 = t220 * t234;
t352 = pkin(3) * qJD(3);
t284 = qJD(2) * t352;
t329 = pkin(3) * qJDD(2);
t180 = t221 * t329 + t226 * t284;
t214 = qJDD(2) + qJDD(3);
t164 = -pkin(5) * t214 + t180;
t247 = -pkin(5) * t215 + t291;
t296 = qJD(1) * qJD(2);
t185 = qJDD(1) * t227 - t222 * t296;
t295 = qJDD(1) * t222;
t186 = -t227 * t296 - t295;
t242 = t182 * qJD(3);
t117 = qJD(1) * t242 + t185 * t226 + t186 * t221;
t118 = qJD(3) * t386 - t185 * t221 + t226 * t186;
t217 = qJDD(1) * pkin(2);
t166 = -pkin(3) * t186 + t217;
t58 = -pkin(4) * t118 + pkin(5) * t117 + t166;
t28 = t220 * t164 - qJD(4) * t122 + (qJD(4) * t247 + t58) * t225;
t304 = qJD(4) * t225;
t91 = t220 * t247 + t225 * t234;
t246 = t220 * t28 + t91 * t304;
t306 = qJD(4) * t220;
t92 = t225 * t247 - t122;
t394 = t92 * t306 - t246;
t393 = t219 * t304 + t220 * t302;
t202 = pkin(3) * t221 - pkin(5);
t203 = pkin(3) * t226 + pkin(4);
t288 = t226 * t352;
t268 = t225 * t288;
t289 = t221 * t352;
t104 = -t202 * t385 + t203 * t302 - t219 * t289 + t224 * t268;
t130 = pkin(3) * t307 + t140;
t208 = pkin(6) * t306;
t392 = -(-pkin(6) * t173 - t130 * t224) * t220 + t208 + t104;
t269 = t220 * t288;
t368 = pkin(6) * t224;
t272 = t202 - t368;
t391 = -t130 * t225 + t272 * t304 + t269 + t395;
t320 = t173 * t220;
t390 = pkin(6) * t320 + t208 - t330;
t120 = t140 * t225 + t220 * t290;
t281 = -pkin(5) - t368;
t389 = t281 * t304 - t120 + t395;
t148 = -t220 * t215 - t225 * t386;
t298 = -qJD(4) - t173;
t109 = -t148 * t219 - t224 * t298;
t108 = qJD(6) + t109;
t348 = Ifges(7,3) * t108;
t378 = Ifges(6,2) * t109 + t348;
t110 = t148 * t224 - t219 * t298;
t114 = qJDD(4) + t118;
t69 = -qJD(4) * t250 + t117 * t225 - t214 * t220;
t33 = -qJD(5) * t110 + t224 * t114 - t219 * t69;
t31 = qJDD(6) + t33;
t362 = t31 * Ifges(7,3);
t388 = t33 * Ifges(6,2) + t362;
t387 = m(6) * t91;
t218 = sin(qJ(6));
t223 = cos(qJ(6));
t212 = t219 * pkin(4);
t172 = t212 + (-pkin(5) * t224 - pkin(6)) * t225;
t187 = t281 * t220;
t251 = t172 * t223 - t187 * t218;
t384 = qJD(6) * t251 + t218 * t390 + t223 * t389;
t136 = t172 * t218 + t187 * t223;
t383 = -qJD(6) * t136 - t218 * t389 + t223 * t390;
t156 = t202 * t309 + t219 * t203;
t149 = -pkin(6) * t225 + t156;
t170 = t272 * t220;
t252 = t149 * t223 - t170 * t218;
t382 = qJD(6) * t252 + t218 * t392 + t391 * t223;
t115 = t149 * t218 + t170 * t223;
t381 = -qJD(6) * t115 - t391 * t218 + t223 * t392;
t179 = -t221 * t284 + t226 * t329;
t163 = pkin(4) * t214 + t179;
t27 = -qJD(4) * t91 + t225 * t164 - t220 * t58;
t189 = pkin(4) * t215 + t290;
t78 = t219 * t189 + t224 * t92;
t16 = -qJD(5) * t78 + t163 * t224 - t219 * t27;
t315 = t219 * t225;
t77 = t189 * t224 - t219 * t92;
t380 = t16 * (pkin(4) * t224 + pkin(5) * t315) + (-pkin(4) * t303 + (-t219 * t306 + t224 * t301) * pkin(5) + t121 * t219 + t224 * t291) * t77;
t316 = t219 * t220;
t286 = t77 * t316;
t379 = -t130 * t286 + t16 * (-t202 * t315 + t203 * t224) + t77 * ((-t202 * t301 - t289) * t224 + (-qJD(5) * t203 + t202 * t306 - t268) * t219);
t70 = -qJD(4) * t148 - t117 * t220 - t214 * t225;
t321 = t148 * t225;
t324 = t250 * t220;
t253 = -t321 - t324;
t332 = t69 * t220;
t335 = t225 * t70;
t377 = qJD(4) * t253 - t332 - t335;
t376 = m(5) / 0.2e1;
t375 = -m(6) / 0.2e1;
t374 = m(6) / 0.2e1;
t373 = -m(5) - m(6);
t372 = pkin(2) * mrSges(3,1);
t370 = -t386 / 0.2e1;
t369 = pkin(6) * t109;
t367 = m(5) * t163;
t74 = t223 * t110 - t218 * t244;
t366 = Ifges(7,1) * t74;
t73 = t218 * t110 + t223 * t244;
t365 = Ifges(7,2) * t73;
t216 = qJ(2) + qJ(3);
t210 = sin(t216);
t211 = cos(t216);
t364 = g(3) * (t210 * t219 - t211 * t309);
t363 = g(3) * t211;
t67 = qJDD(5) - t70;
t360 = t67 * Ifges(6,3);
t359 = t78 * mrSges(6,2);
t358 = t91 * mrSges(6,2);
t357 = t92 * mrSges(5,3);
t351 = mrSges(4,3) * t173;
t350 = mrSges(4,3) * t386;
t32 = qJD(5) * t109 + t114 * t219 + t224 * t69;
t327 = qJD(6) * t73;
t10 = t218 * t67 - t223 * t32 + t327;
t347 = t10 * t223;
t326 = qJD(6) * t74;
t11 = t218 * t32 + t223 * t67 + t326;
t346 = t11 * t218;
t345 = t110 * Ifges(6,1);
t344 = t114 * Ifges(5,3);
t343 = t120 * t91;
t341 = t298 * Ifges(5,3);
t340 = t386 * Ifges(4,4);
t339 = t218 * t73;
t337 = t223 * t74;
t336 = t224 * t91;
t334 = t225 * t91;
t333 = t27 * t225;
t331 = t91 * t128;
t144 = qJD(2) * t182 + t242;
t325 = t144 * t220;
t323 = t250 * t224;
t317 = t218 * t220;
t314 = t220 * t223;
t313 = t220 * t224;
t311 = t223 * t224;
t310 = t223 * t225;
t308 = t227 * qJD(1) ^ 2;
t300 = qJD(6) * t218;
t299 = qJD(6) * t223;
t294 = Ifges(3,2) - Ifges(3,1) / 0.2e1;
t293 = t220 * t359;
t292 = t220 * t363;
t287 = Ifges(5,2) * t324;
t285 = t220 * t226 * t91;
t277 = t224 * t304;
t249 = -t226 * t227 + t312;
t141 = -pkin(4) * t182 - pkin(5) * t249 + t356;
t145 = t215 * t249;
t82 = -pkin(4) * t145 + pkin(5) * t144 + t227 * t353;
t271 = t82 * t286 + (t16 * t316 + t393 * t77) * t141;
t267 = -qJD(6) * t224 - qJD(4);
t266 = qJD(6) + t305;
t265 = -pkin(4) * t210 - pkin(5) * t211;
t236 = -qJD(5) * t182 - t144 * t225 - t249 * t306;
t261 = t249 * t301 + t145;
t52 = t219 * t261 - t224 * t236;
t264 = -qJD(6) * t220 * t249 - t52;
t263 = Ifges(3,1) * t307 / 0.2e1 + Ifges(3,5) * qJD(2);
t262 = -t333 + t363;
t260 = pkin(6) * t249 - t141 * t224;
t230 = t110 * pkin(6) + t91;
t59 = -pkin(6) * t244 + t78;
t22 = t218 * t59 + t223 * t230;
t23 = -t218 * t230 + t223 * t59;
t259 = -t218 * t23 + t22 * t223;
t257 = -t220 * t92 + t334;
t103 = t260 * t220;
t135 = t182 * t219 - t249 * t309;
t80 = pkin(6) * t135 + t141 * t225;
t256 = t103 * t223 - t218 * t80;
t45 = t103 * t218 + t223 * t80;
t254 = t148 * t220 - t225 * t250;
t248 = -t213 + t265;
t245 = t225 * t28 - t306 * t91;
t15 = t219 * t163 + t189 * t302 + t224 * t27 - t303 * t92;
t243 = -Ifges(7,1) * t337 - Ifges(7,2) * t339;
t161 = t210 * t315 - t211 * t224;
t241 = mrSges(4,1) * t210 + mrSges(4,2) * t211 + (-mrSges(6,2) - mrSges(7,3)) * t161;
t124 = t266 * t310 + (t218 * t267 - t223 * t303) * t220;
t176 = t218 * t225 + t220 * t311;
t6 = -pkin(6) * t67 + t15;
t7 = t32 * pkin(6) + t28;
t3 = qJD(6) * t23 + t218 * t6 + t223 * t7;
t86 = t128 * t218 + t173 * t314;
t240 = -t124 * t22 - t176 * t3 + t23 * t86;
t238 = t224 * t244;
t2 = qJD(6) * t22 + t218 * t7 - t223 * t6;
t237 = t2 * t218 + t223 * t3 - t292;
t235 = qJD(6) * t135 - t249 * t304 + t325;
t233 = -Ifges(6,3) * t244 + t359;
t232 = -t333 + t394;
t231 = t257 * t376 + (-t313 * t78 + t334) * t374;
t125 = t267 * t314 + (-t225 * t266 + t276) * t218;
t131 = t173 * Ifges(4,2) + t215 * Ifges(4,6) - t340;
t167 = Ifges(4,4) * t173;
t132 = -Ifges(4,1) * t386 + Ifges(4,5) * t215 + t167;
t175 = -t218 * t313 + t310;
t87 = -t128 * t223 + t173 * t317;
t229 = (t293 - t287) * qJD(4) + (-t125 * t23 + t175 * t2 + t22 * t87) * mrSges(7,3) + (t173 * t315 - t224 * t386 + t393) * t378 + (-t173 * t334 + t394) * mrSges(5,3) + t173 * t293 - t173 * t287 - (Ifges(4,2) * t386 + t132 + t167) * t173 / 0.2e1 - t193 * (-mrSges(4,1) * t386 + mrSges(4,2) * t173) + t388 * t316 + t131 * t370 + t290 * t351 + t320 * t357 + t225 * t360 + t298 * t396 + Ifges(4,3) * t214 + t179 * mrSges(4,1) - t180 * mrSges(4,2) + (-t124 + t87) * t366 - Ifges(5,2) * t335 + Ifges(4,5) * t117 + Ifges(4,6) * t118 + (t298 * t321 - t332) * Ifges(5,1) + (t125 - t86) * t365 + t11 * Ifges(7,2) * t175 + t10 * Ifges(7,1) * t176 - t215 * (Ifges(4,5) * t173 + Ifges(4,6) * t386) / 0.2e1 + t386 * (Ifges(4,1) * t173 + t340) / 0.2e1 + t386 * t341 + (t345 + t358) * t276 + (-t32 * t313 + (-t128 - t277) * t110) * Ifges(6,1);
t184 = -pkin(5) * t309 + t212;
t159 = t161 * pkin(6);
t151 = mrSges(4,1) * t215 + t350;
t150 = -mrSges(4,2) * t215 + t351;
t139 = -mrSges(4,1) * t173 - mrSges(4,2) * t386;
t134 = t182 * t224 + t249 * t315;
t100 = -t135 * t223 - t249 * t317;
t99 = t135 * t218 - t249 * t314;
t89 = t148 * t218 + t250 * t311;
t88 = t148 * t223 - t218 * t323;
t66 = -pkin(6) * t323 + t92;
t64 = -pkin(6) * t148 - t336;
t53 = t219 * t236 + t224 * t261;
t48 = t218 * t369 - t223 * t77;
t47 = t218 * t77 + t223 * t369;
t35 = t218 * t66 - t223 * t64;
t34 = t218 * t64 + t223 * t66;
t24 = t260 * t304 + (-pkin(6) * t144 + t141 * t303 - t224 * t82) * t220;
t21 = pkin(6) * t52 - t141 * t306 + t225 * t82;
t20 = -t218 * t264 + t223 * t235;
t19 = t218 * t235 + t223 * t264;
t5 = qJD(6) * t256 + t21 * t223 + t218 * t24;
t4 = qJD(6) * t45 + t21 * t218 - t223 * t24;
t1 = [(m(5) * (-t220 * t27 - t304 * t92 + t245) + m(6) * (-t15 * t313 + t245 + (t276 - t277) * t78) + (-qJD(4) * t254 - t220 * t70 + t225 * t69) * mrSges(5,3) + (-t239 * t303 + t225 * t32 + t67 * t313 + (-t220 * t110 + t225 * t238) * qJD(4)) * mrSges(6,2)) * t141 + (t109 * t53 + t134 * t33) * Ifges(6,2) + (-t10 * t45 - t100 * t3 - t11 * t256 - t19 * t22 + t2 * t99 - t20 * t23 + t4 * t73 + t5 * t74) * mrSges(7,3) + m(7) * (-t2 * t256 + t22 * t5 - t23 * t4 + t3 * t45 + t271) + (t10 * t100 - t19 * t74) * Ifges(7,1) - (t166 * mrSges(4,2) - t179 * mrSges(4,3) + Ifges(4,1) * t117 + Ifges(4,4) * t118 + Ifges(4,5) * t214 + (mrSges(5,3) * t28 + t69 * Ifges(5,1)) * t225 + (-t15 * mrSges(6,2) - mrSges(5,3) * t27 - Ifges(5,2) * t70 + t360) * t220 + ((-mrSges(5,3) * t91 - Ifges(5,1) * t148) * t220 + (Ifges(5,2) * t250 - t233 - t357) * t225) * qJD(4)) * t249 + (t144 * t257 - t253 * t82) * mrSges(5,3) + Ifges(5,1) * t144 * t321 + (t110 * t52 + t135 * t32) * Ifges(6,1) + (t11 * t99 + t20 * t73) * Ifges(7,2) + (t108 * t53 + t134 * t31) * Ifges(7,3) + 0.2e1 * t231 * t82 + (Ifges(4,1) * t144 + Ifges(4,4) * t145) * t370 - t145 * t341 + t144 * t396 + Ifges(2,3) * qJDD(1) + t215 * (Ifges(4,5) * t144 + Ifges(4,6) * t145) / 0.2e1 + (Ifges(3,1) * t185 + Ifges(3,5) * qJDD(2) + ((m(4) * t193 + t139) * pkin(3) + (t294 * t222 + t372) * qJD(1)) * qJD(2)) * t227 + t193 * (-mrSges(4,1) * t145 + mrSges(4,2) * t144) + t173 * (Ifges(4,4) * t144 + Ifges(4,2) * t145) / 0.2e1 + (m(4) * t166 - mrSges(4,1) * t118 + mrSges(4,2) * t117) * t356 + (-mrSges(4,1) * t166 + mrSges(4,3) * t180 + Ifges(4,4) * t117 + Ifges(4,2) * t118 + Ifges(4,6) * t214 + t344) * t182 + t144 * t132 / 0.2e1 + t145 * t131 / 0.2e1 + (-t78 * t325 + t28 * t135 + t91 * t52 + (t225 * t110 + t220 * t238) * t82) * mrSges(6,2) + (m(3) * t217 + (-t186 + t295) * mrSges(3,1)) * pkin(2) + (-t263 * t222 + (-t144 * t226 + t145 * t221) * mrSges(4,3) * pkin(3)) * qJD(2) - t186 * Ifges(3,2) * t222 + m(6) * t271 + t144 * t287; (t253 * t130 - t202 * t377 + t262) * mrSges(5,3) + (-t331 - t156 * t67 - t15 * t225 - t104 * t244 + t220 * t202 * t32 + (-t130 * t239 - t246) * t224 + (t269 + (qJD(4) * t202 - t130) * t225) * t110) * mrSges(6,2) + (t222 * mrSges(3,1) + t373 * t248 + t241) * g(3) + t203 * t367 + (-t115 * t10 - t11 * t252 - t381 * t73 + t382 * t74 + t240) * mrSges(7,3) + 0.2e1 * (t246 * t374 - m(5) * t232 / 0.2e1) * t202 - 0.2e1 * t231 * t130 + (qJD(1) * t263 - t294 * t308) * t222 + t229 + Ifges(3,5) * t185 - t308 * t372 + Ifges(3,3) * qJDD(2) + (t115 * t3 - t252 * t2 + (-t159 - t248) * g(3) + t381 * t23 + t382 * t22 + t379) * m(7) + (t104 * t78 + t15 * t156 + t379) * m(6) + (0.2e1 * (t285 * t374 + (t225 * t226 * t92 - t189 * t221 + t285) * t376) * qJD(3) + t226 * (mrSges(4,1) * t214 - mrSges(4,3) * t117) - t139 * t307 + (-mrSges(4,2) * t214 + (-qJD(2) * t386 + t118) * mrSges(4,3)) * t221 + (-t221 * t151 + (t254 * mrSges(5,3) + t150) * t226) * qJD(3) + (g(3) * t222 + t179 * t226 + t180 * t221 - t193 * t307) * m(4)) * pkin(3); (-t10 * t136 - t11 * t251 - t383 * t73 + t384 * t74 + t240) * mrSges(7,3) + 0.2e1 * (t232 * t376 + t246 * t375) * pkin(5) + (t373 * t265 + t241) * g(3) + pkin(4) * t367 + (-t120 * t110 - t331 - t184 * t67 + t330 * qJD(5) + (-pkin(5) * t32 - t28 * t224 - t330 * t386) * t220 + (-t15 + t330 * t215 + (-pkin(5) * t110 - t336) * qJD(4)) * t225) * mrSges(6,2) + t229 + (-t150 * t226 + (t151 - t350) * t221) * t353 + (pkin(5) * t377 - t120 * t148 + t121 * t250 + t262) * mrSges(5,3) - m(5) * (t121 * t92 - t189 * t291 + t343) + ((-t159 - t265) * g(3) + t136 * t3 - t251 * t2 + t383 * t23 + t384 * t22 + t380) * m(7) + (t15 * t184 - t330 * t78 - t343 + t380) * m(6); t89 * t366 - t88 * t365 + t344 - m(7) * (t22 * t34 - t23 * t35) + (-(-Ifges(5,1) + Ifges(5,2)) * t250 + t233) * t148 + (t22 * t89 + t23 * t88 - t34 * t74 - t35 * t73) * mrSges(7,3) + (t250 * t345 + t78 * t387 + (t345 + t259 * mrSges(7,3) + (m(7) * t259 + (t337 + t339) * mrSges(7,3)) * pkin(6) - t243) * qJD(5) + t388) * t224 + (t32 * Ifges(6,1) + (t299 * t73 + t346) * Ifges(7,2) + (-t300 * t74 - t347) * Ifges(7,1) + (t28 - t292) * mrSges(6,2) + 0.2e1 * (t375 - m(7) / 0.2e1) * t91 * t77 + ((-t218 * t22 - t223 * t23) * qJD(6) + t237) * mrSges(7,3) + ((-t22 * t300 - t23 * t299 + t237) * m(7) + (-t347 + t346 + (-t218 * t74 + t223 * t73) * qJD(6)) * mrSges(7,3)) * pkin(6) - t378 * t244) * t219 + (-mrSges(6,2) * t110 - t387) * t92; t110 * t348 + t360 + (t11 * t223 - t300 * t73) * Ifges(7,2) + (t10 * t218 - t299 * t74) * Ifges(7,1) + (t244 * t77 - t15 - t364) * mrSges(6,2) + (-t358 + (-Ifges(6,1) + Ifges(6,2)) * t110 + t243) * t109 + (-t364 - t47 * t74 - t48 * t73 + (t2 - t108 * t22 + (t11 - t326) * pkin(6)) * t223 + (-t3 + t108 * t23 + (t10 - t327) * pkin(6)) * t218) * mrSges(7,3) + (-t22 * t47 + t23 * t48 + t78 * t77 + (t2 * t223 - t218 * t3 - t22 * t299 + t23 * t300 - t364) * pkin(6)) * m(7); t362 + (Ifges(7,1) - Ifges(7,2)) * t74 * t73;];
tau = t1;
