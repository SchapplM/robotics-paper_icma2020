% Calculate vector of centrifugal and Coriolis load on the joints for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = CloosQRC350OL_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:00:27
% EndTime: 2020-06-20 08:02:03
% DurationCPUTime: 56.52s
% Computational Cost: add. (17343->794), mult. (40829->1152), div. (0->0), fcn. (31393->10), ass. (0->352)
t300 = sin(qJ(3));
t301 = sin(qJ(2));
t305 = cos(qJ(3));
t306 = cos(qJ(2));
t272 = -t300 * t306 - t305 * t301;
t395 = qJD(2) + qJD(3);
t242 = t395 * t272;
t230 = t242 * qJD(1);
t406 = qJD(1) * t306;
t415 = t300 * t301;
t269 = -qJD(1) * t415 + t305 * t406;
t299 = sin(qJ(4));
t366 = t299 * t395;
t304 = cos(qJ(4));
t400 = qJD(4) * t304;
t173 = -qJD(4) * t366 + t230 * t299 + t269 * t400;
t268 = t272 * qJD(1);
t292 = t301 * pkin(3) + pkin(2);
t282 = t292 * qJD(1);
t313 = -pkin(4) * t268 + pkin(5) * t269 + t282;
t211 = t299 * t313;
t469 = pkin(3) * qJD(2);
t391 = t300 * t469;
t317 = -pkin(5) * t395 + t391;
t185 = t304 * t317 - t211;
t390 = t305 * t469;
t278 = pkin(4) * t395 + t390;
t298 = sin(qJ(5));
t303 = cos(qJ(5));
t468 = pkin(3) * qJD(3);
t389 = qJD(2) * t468;
t362 = t300 * t389;
t397 = qJD(5) * t303;
t330 = -t305 * t306 + t415;
t243 = t395 * t330;
t231 = t243 * qJD(1);
t405 = qJD(2) * t306;
t394 = pkin(3) * t405;
t152 = -t231 * pkin(4) + t230 * pkin(5) + qJD(1) * t394;
t361 = t305 * t389;
t184 = t299 * t317 + t304 * t313;
t403 = qJD(4) * t184;
t91 = -t299 * t152 + t304 * t361 - t403;
t46 = t278 * t397 + t303 * t91 + (-qJD(5) * t185 - t362) * t298;
t29 = -pkin(6) * t173 + t46;
t297 = sin(qJ(6));
t247 = -t299 * t269 - t304 * t395;
t172 = qJD(4) * t247 + t304 * t230;
t248 = t269 * t304 - t366;
t262 = qJD(4) + t268;
t200 = -t248 * t298 + t303 * t262;
t85 = qJD(5) * t200 + t172 * t303 + t231 * t298;
t92 = t299 * t361 - qJD(4) * t211 + (qJD(4) * t317 + t152) * t304;
t30 = t85 * pkin(6) + t92;
t302 = cos(qJ(6));
t158 = t185 * t303 + t278 * t298;
t246 = qJD(5) - t247;
t126 = -pkin(6) * t246 + t158;
t201 = t248 * t303 + t262 * t298;
t311 = t201 * pkin(6) + t184;
t56 = t297 * t126 + t302 * t311;
t3 = qJD(6) * t56 - t302 * t29 + t297 * t30;
t57 = t302 * t126 - t297 * t311;
t4 = qJD(6) * t57 + t297 * t29 + t302 * t30;
t358 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t86 = -qJD(5) * t201 - t172 * t298 + t231 * t303;
t507 = t86 / 0.2e1;
t495 = t173 / 0.2e1;
t508 = t85 / 0.2e1;
t563 = Ifges(6,4) * t508 + Ifges(6,6) * t495;
t151 = t201 * t302 - t297 * t246;
t28 = qJD(6) * t151 + t173 * t302 + t297 * t85;
t25 = Ifges(7,6) * t28;
t150 = t201 * t297 + t246 * t302;
t27 = qJD(6) * t150 + t173 * t297 - t302 * t85;
t26 = Ifges(7,5) * t27;
t8 = Ifges(7,3) * t86 + t25 + t26;
t559 = Ifges(6,2) * t507 + t8 / 0.2e1 + t563;
t567 = -t92 * mrSges(6,1) + t46 * mrSges(6,3) + t358 + t559 + t563;
t412 = t303 * t304;
t216 = t268 * t412 - t269 * t298;
t398 = qJD(5) * t299;
t381 = t298 * t398;
t566 = -t216 + t381;
t565 = t566 * pkin(6);
t238 = t269 * pkin(4) + t268 * pkin(5);
t220 = pkin(3) * t406 + t238;
t290 = pkin(3) * t300 - pkin(5);
t475 = pkin(6) * t303;
t372 = t290 - t475;
t392 = t305 * t468;
t562 = t220 * t304 - t299 * t392 - t372 * t400 - t565;
t291 = pkin(3) * t305 + pkin(4);
t365 = t304 * t392;
t404 = qJD(3) * t300;
t393 = pkin(3) * t404;
t396 = qJD(5) * t304;
t401 = qJD(4) * t303;
t542 = t298 * t396 + t299 * t401;
t194 = -t290 * t542 + t291 * t397 - t298 * t393 + t303 * t365;
t402 = qJD(4) * t299;
t294 = pkin(6) * t402;
t561 = -(-pkin(6) * t268 - t220 * t303) * t299 + t294 + t194;
t210 = -t238 * t299 + t304 * t390;
t526 = pkin(4) * t397 + pkin(5) * t542 - t303 * t210 + t298 * t391;
t209 = t238 * t304 + t299 * t390;
t386 = -pkin(5) - t475;
t558 = t386 * t400 - t209 + t565;
t424 = t268 * t299;
t557 = -pkin(6) * t424 - t294 - t526;
t399 = qJD(5) * t298;
t556 = t210 * t298 + t303 * t391 - pkin(4) * t399 + (-t298 * t402 + t303 * t396) * pkin(5);
t439 = t299 * t92;
t327 = t184 * t400 + t439;
t555 = -t185 * t402 + t304 * t91 + t327;
t448 = t173 * Ifges(6,5);
t472 = t86 * Ifges(6,4);
t474 = t85 * Ifges(6,1);
t21 = t448 + t472 + t474;
t553 = t92 * mrSges(6,2) + t21 / 0.2e1;
t128 = -mrSges(5,2) * t231 - mrSges(5,3) * t173;
t206 = -mrSges(5,2) * t262 + mrSges(5,3) * t247;
t436 = mrSges(5,1) * t231 + mrSges(6,1) * t86 - mrSges(6,2) * t85 - mrSges(5,3) * t172;
t409 = -mrSges(5,1) * t262 - mrSges(6,1) * t200 + mrSges(6,2) * t201 + mrSges(5,3) * t248;
t527 = t409 * t304;
t552 = m(5) * t555 + qJD(4) * (-t299 * t206 + t527) + t304 * t128 - t436 * t299;
t451 = Ifges(6,6) * t246;
t452 = Ifges(6,2) * t200;
t460 = Ifges(6,4) * t201;
t116 = t451 + t452 + t460;
t199 = qJD(6) + t200;
t449 = Ifges(7,3) * t199;
t450 = Ifges(7,6) * t150;
t453 = Ifges(7,5) * t151;
t58 = t449 + t450 - t453;
t551 = t158 * mrSges(6,3) + t56 * mrSges(7,1) + t57 * mrSges(7,2) + t116 / 0.2e1 + t58 / 0.2e1 - t184 * mrSges(6,1);
t258 = t290 * t412 + t298 * t291;
t249 = -pkin(6) * t304 + t258;
t263 = t372 * t299;
t332 = t249 * t302 - t263 * t297;
t549 = qJD(6) * t332 + t561 * t297 - t302 * t562;
t203 = t249 * t297 + t263 * t302;
t548 = -qJD(6) * t203 + t297 * t562 + t561 * t302;
t464 = Ifges(4,4) * t269;
t547 = t248 * Ifges(5,5) + Ifges(4,2) * t268 + t395 * Ifges(4,6) + t247 * Ifges(5,6) + t262 * Ifges(5,3) + t464;
t417 = t299 * t302;
t177 = t216 * t297 + t268 * t417;
t363 = qJD(6) + t401;
t364 = -qJD(6) * t303 - qJD(4);
t213 = t364 * t417 + (-t304 * t363 + t381) * t297;
t546 = t177 - t213;
t421 = t297 * t299;
t178 = -t216 * t302 + t268 * t421;
t413 = t302 * t304;
t212 = t363 * t413 + (t297 * t364 - t302 * t399) * t299;
t545 = t178 - t212;
t419 = t298 * t304;
t215 = -t268 * t419 - t269 * t303;
t318 = t298 * t400 + t299 * t397;
t544 = t215 - t318;
t543 = t303 * t400 - t566;
t541 = t402 + t424;
t103 = -mrSges(7,2) * t199 + mrSges(7,3) * t150;
t104 = mrSges(7,1) * t199 + mrSges(7,3) * t151;
t337 = t297 * t57 - t302 * t56;
t326 = t337 * mrSges(7,3);
t537 = -(m(7) * t337 - t297 * t103 - t302 * t104) * pkin(6) - t326;
t153 = -mrSges(6,2) * t246 + mrSges(6,3) * t200;
t434 = mrSges(6,1) * t246 - mrSges(7,1) * t150 - mrSges(7,2) * t151 - mrSges(6,3) * t201;
t316 = -t303 * t153 + t298 * t434 - t206;
t428 = t184 * t304;
t334 = -t185 * t299 + t428;
t416 = t299 * t303;
t536 = -m(5) * t334 - m(6) * (-t158 * t416 + t428) - t527 - t299 * t316;
t535 = (m(4) * t282 - mrSges(4,1) * t268 + mrSges(4,2) * t269) * pkin(3);
t534 = t46 * mrSges(6,2);
t532 = -mrSges(5,1) * t304 + mrSges(5,2) * t299 - mrSges(4,1);
t295 = t298 * pkin(4);
t265 = t295 + (-pkin(5) * t303 - pkin(6)) * t304;
t276 = t386 * t299;
t331 = t265 * t302 - t276 * t297;
t531 = qJD(6) * t331 - t297 * t557 + t302 * t558;
t234 = t265 * t297 + t276 * t302;
t530 = qJD(6) * t234 + t297 * t558 + t302 * t557;
t352 = mrSges(5,1) * t299 + mrSges(5,2) * t304;
t529 = t278 * t352;
t157 = -t185 * t298 + t278 * t303;
t274 = pkin(4) * t303 + pkin(5) * t419;
t47 = -qJD(5) * t158 - t298 * t91 - t303 * t362;
t525 = t157 * t556 + t47 * t274;
t454 = Ifges(6,5) * t246;
t198 = Ifges(6,4) * t200;
t467 = Ifges(6,1) * t201;
t117 = t198 + t454 + t467;
t505 = -t117 / 0.2e1;
t322 = t505 - t454 / 0.2e1 - t184 * mrSges(6,2);
t350 = mrSges(7,1) * t297 + mrSges(7,2) * t302;
t329 = mrSges(6,3) + t350;
t339 = Ifges(7,5) * t302 - Ifges(7,6) * t297;
t455 = Ifges(7,4) * t302;
t343 = -Ifges(7,2) * t297 + t455;
t456 = Ifges(7,4) * t297;
t347 = Ifges(7,1) * t302 - t456;
t479 = t302 / 0.2e1;
t492 = t199 / 0.2e1;
t501 = -t151 / 0.2e1;
t502 = t150 / 0.2e1;
t457 = Ifges(7,4) * t151;
t59 = Ifges(7,2) * t150 + Ifges(7,6) * t199 - t457;
t513 = -t59 / 0.2e1;
t145 = Ifges(7,4) * t150;
t60 = -Ifges(7,1) * t151 + Ifges(7,5) * t199 + t145;
t521 = t297 * t513 + t339 * t492 + t343 * t502 + t347 * t501 + t479 * t60;
t523 = t329 * t157 + t322 + t521;
t520 = -m(6) / 0.2e1;
t9 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t86;
t519 = t9 / 0.2e1;
t518 = Ifges(6,2) / 0.2e1;
t10 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t86;
t517 = t10 / 0.2e1;
t515 = t27 / 0.2e1;
t514 = t28 / 0.2e1;
t512 = t59 / 0.2e1;
t511 = -t60 / 0.2e1;
t510 = t60 / 0.2e1;
t509 = -Ifges(5,4) * t172 / 0.2e1 + Ifges(5,2) * t495 - Ifges(5,6) * t231 / 0.2e1;
t444 = t246 * Ifges(6,3);
t445 = t201 * Ifges(6,5);
t446 = t200 * Ifges(6,6);
t115 = t444 + t445 + t446;
t506 = -t115 / 0.2e1;
t503 = -t150 / 0.2e1;
t500 = t151 / 0.2e1;
t442 = t262 * Ifges(5,6);
t463 = Ifges(5,4) * t248;
t164 = t247 * Ifges(5,2) + t442 + t463;
t499 = t164 / 0.2e1;
t245 = Ifges(5,4) * t247;
t443 = t262 * Ifges(5,5);
t165 = Ifges(5,1) * t248 + t245 + t443;
t498 = -t165 / 0.2e1;
t497 = t172 / 0.2e1;
t496 = -t173 / 0.2e1;
t494 = -t198 / 0.2e1;
t493 = -t199 / 0.2e1;
t491 = t200 / 0.2e1;
t490 = t201 / 0.2e1;
t488 = t231 / 0.2e1;
t485 = t246 / 0.2e1;
t483 = t268 / 0.2e1;
t482 = t269 / 0.2e1;
t481 = t297 / 0.2e1;
t480 = t301 / 0.2e1;
t478 = t304 / 0.2e1;
t477 = -t306 / 0.2e1;
t476 = pkin(6) * t200;
t470 = mrSges(6,1) * t173 - mrSges(7,1) * t28 + mrSges(7,2) * t27 - mrSges(6,3) * t85;
t466 = Ifges(3,4) * t301;
t465 = Ifges(3,4) * t306;
t462 = Ifges(5,4) * t299;
t461 = Ifges(5,4) * t304;
t459 = Ifges(6,4) * t298;
t458 = Ifges(6,4) * t303;
t441 = t268 * mrSges(4,3);
t440 = t269 * mrSges(4,3);
t195 = (-t290 * t396 - t393) * t303 + (-qJD(5) * t291 + t290 * t402 - t365) * t298;
t257 = -t290 * t419 + t291 * t303;
t437 = t157 * t195 + t47 * t257;
t433 = Ifges(3,5) * qJD(2);
t432 = Ifges(3,6) * qJD(2);
t430 = t184 * t209;
t429 = t184 * t299;
t427 = t242 * t299;
t426 = t247 * t303;
t423 = t268 * t304;
t422 = t330 * t299;
t420 = t298 * t299;
t414 = t302 * t303;
t410 = t304 * t242;
t407 = mrSges(4,1) * t395 - mrSges(5,1) * t247 + mrSges(5,2) * t248 - t440;
t19 = Ifges(6,5) * t85 + Ifges(6,6) * t86 + Ifges(6,3) * t173;
t388 = t157 * t420;
t387 = Ifges(5,5) * t172 - Ifges(5,6) * t173 + Ifges(5,3) * t231;
t378 = t299 * t499;
t376 = -t402 / 0.2e1;
t375 = -t400 / 0.2e1;
t374 = t400 / 0.2e1;
t373 = t398 / 0.2e1;
t171 = -t243 * pkin(4) + t242 * pkin(5) + t394;
t239 = -pkin(4) * t272 - pkin(5) * t330 + t292;
t367 = t171 * t388 + (t157 * t318 + t420 * t47) * t239;
t357 = -t297 * t4 + t3 * t302;
t356 = t297 * t3 + t302 * t4;
t355 = t330 * t396 + t243;
t324 = t330 * t402 + t410;
t315 = -qJD(5) * t272 - t324;
t118 = t298 * t355 - t303 * t315;
t354 = -qJD(6) * t422 - t118;
t353 = pkin(6) * t330 - t239 * t303;
t351 = -mrSges(7,1) * t302 + mrSges(7,2) * t297;
t349 = Ifges(5,1) * t304 - t462;
t348 = -Ifges(6,1) * t303 + t459;
t346 = Ifges(7,1) * t297 + t455;
t345 = -Ifges(5,2) * t299 + t461;
t344 = Ifges(6,2) * t298 - t458;
t342 = Ifges(7,2) * t302 + t456;
t341 = Ifges(5,5) * t304 - Ifges(5,6) * t299;
t340 = -Ifges(6,5) * t303 + Ifges(6,6) * t298;
t338 = Ifges(7,5) * t297 + Ifges(7,6) * t302;
t336 = -t297 * t56 - t302 * t57;
t233 = t272 * t298 - t330 * t412;
t166 = pkin(6) * t233 + t239 * t304;
t192 = t353 * t299;
t108 = t166 * t302 + t192 * t297;
t109 = t166 * t297 - t192 * t302;
t333 = t185 * t304 + t429;
t325 = -t330 * t400 + t427;
t314 = qJD(6) * t233 + t325;
t310 = t451 / 0.2e1 + t450 / 0.2e1 - t453 / 0.2e1 + t449 / 0.2e1 + t460 / 0.2e1 + t551;
t309 = t452 / 0.2e1 + t310;
t260 = Ifges(4,4) * t268;
t222 = Ifges(4,1) * t269 + t395 * Ifges(4,5) + t260;
t270 = -t297 * t416 + t413;
t271 = t297 * t304 + t299 * t414;
t72 = Ifges(5,1) * t172 - Ifges(5,4) * t173 + Ifges(5,5) * t231;
t307 = t420 * t559 + t3 * (-mrSges(7,2) * t420 + mrSges(7,3) * t270) + t46 * (-mrSges(6,2) * t304 + mrSges(6,3) * t420) + t4 * (mrSges(7,1) * t420 - mrSges(7,3) * t271) - t185 * mrSges(5,2) * t269 + (t116 + t58) * (t298 * t374 + t303 * t373 - t215 / 0.2e1) + (-Ifges(5,5) * t299 - Ifges(5,6) * t304) * t488 + (Ifges(7,5) * t271 + Ifges(6,6) * t304 + Ifges(7,6) * t270 + Ifges(7,3) * t420 + t299 * t344) * t507 - t246 * (Ifges(6,5) * t216 + Ifges(6,6) * t215 + Ifges(6,3) * t424) / 0.2e1 - t200 * (Ifges(6,4) * t216 + Ifges(6,2) * t215 + Ifges(6,6) * t424) / 0.2e1 - t201 * (Ifges(6,1) * t216 + Ifges(6,4) * t215 + Ifges(6,5) * t424) / 0.2e1 + (mrSges(6,1) * t304 - mrSges(7,1) * t270 + mrSges(7,2) * t271 + mrSges(6,3) * t416) * t47 + (-mrSges(7,2) * t544 + mrSges(7,3) * t546) * t57 + (-mrSges(6,1) * t541 + mrSges(7,1) * t546 - mrSges(7,2) * t545 + mrSges(6,3) * t543) * t157 + t547 * t482 + (-mrSges(6,1) * t298 - mrSges(6,2) * t303) * t439 + t390 * t441 - (t247 * t345 + t248 * t349 + t262 * t341) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t269 + t222 + t260) * t268 / 0.2e1 - t282 * (mrSges(4,1) * t269 + mrSges(4,2) * t268) + ((Ifges(6,5) * t298 + Ifges(6,6) * t303) * t485 + (Ifges(6,1) * t298 + t458) * t490 + (Ifges(6,2) * t303 + t459) * t491) * t398 + t19 * t478 + (-mrSges(5,1) * t269 + mrSges(6,1) * t544 - mrSges(6,2) * t543) * t184 + (mrSges(6,2) * t541 - mrSges(6,3) * t544) * t158 + (-mrSges(7,1) * t544 + mrSges(7,3) * t545) * t56 - t269 * (Ifges(4,1) * t268 - t464) / 0.2e1 - t395 * (Ifges(4,5) * t268 - Ifges(4,6) * t269) / 0.2e1 - t21 * t416 / 0.2e1 + (Ifges(7,5) * t212 + Ifges(7,6) * t213 + Ifges(7,3) * t318) * t492 + (Ifges(7,5) * t178 + Ifges(7,6) * t177 + Ifges(7,3) * t215) * t493 + (Ifges(6,3) * t304 + t299 * t340) * t495 + (-Ifges(5,2) * t304 - t462) * t496 + (-Ifges(5,1) * t299 - t461) * t497 + t423 * t498 + (Ifges(7,1) * t178 + Ifges(7,4) * t177 + Ifges(7,5) * t215) * t500 + t165 * t375 + t115 * t376 + t268 * t378 + (-t184 * t423 + t185 * t424 - t555) * mrSges(5,3) - t268 * t529 + (Ifges(7,1) * t212 + Ifges(7,4) * t213 + Ifges(7,5) * t318) * t501 + (Ifges(7,4) * t212 + Ifges(7,2) * t213 + Ifges(7,6) * t318) * t502 + (Ifges(7,4) * t178 + Ifges(7,2) * t177 + Ifges(7,6) * t215) * t503 + t216 * t505 + t424 * t506 + (Ifges(6,5) * t304 + t299 * t348) * t508 + t304 * t509 + t212 * t510 + t178 * t511 + t213 * t512 + t177 * t513 + (Ifges(7,4) * t271 + Ifges(7,2) * t270 + Ifges(7,6) * t420) * t514 + (Ifges(7,1) * t271 + Ifges(7,4) * t270 + Ifges(7,5) * t420) * t515 + t271 * t517 + t270 * t519 + ((-Ifges(6,3) * t299 + t304 * t340) * t485 + (-Ifges(6,5) * t299 + t304 * t348) * t490 + (-Ifges(6,6) * t299 + t304 * t344) * t491 + t378 - t529) * qJD(4) - t262 * (-Ifges(5,3) * t269 + t268 * t341) / 0.2e1 - t247 * (-Ifges(5,6) * t269 + t268 * t345) / 0.2e1 + Ifges(4,5) * t230 + Ifges(4,6) * t231 - t248 * (-Ifges(5,5) * t269 + t268 * t349) / 0.2e1 - t299 * t72 / 0.2e1 + (t298 * t373 + t303 * t375) * t117;
t275 = -pkin(5) * t412 + t295;
t267 = t433 + (t306 * Ifges(3,1) - t466) * qJD(1);
t266 = t432 + (-Ifges(3,2) * t301 + t465) * qJD(1);
t251 = -mrSges(4,2) * t395 + t441;
t189 = -t233 * t302 - t330 * t421;
t188 = t233 * t297 - t330 * t417;
t182 = -t247 * t414 + t248 * t297;
t181 = t248 * t302 + t297 * t426;
t142 = pkin(6) * t426 + t185;
t138 = -pkin(6) * t248 - t184 * t303;
t120 = t220 * t388;
t111 = -t157 * t302 + t297 * t476;
t110 = t157 * t297 + t302 * t476;
t107 = mrSges(5,1) * t173 + mrSges(5,2) * t172;
t66 = -t138 * t302 + t142 * t297;
t65 = t138 * t297 + t142 * t302;
t61 = t353 * t400 + (-pkin(6) * t242 - t171 * t303 + t239 * t399) * t299;
t54 = -mrSges(6,2) * t173 + mrSges(6,3) * t86;
t50 = pkin(6) * t118 + t171 * t304 - t239 * t402;
t49 = -t297 * t354 + t302 * t314;
t48 = t297 * t314 + t302 * t354;
t15 = -mrSges(7,2) * t86 + mrSges(7,3) * t28;
t14 = mrSges(7,1) * t86 - mrSges(7,3) * t27;
t13 = -qJD(6) * t109 + t297 * t61 + t302 * t50;
t12 = qJD(6) * t108 + t297 * t50 - t302 * t61;
t1 = [(t188 * t3 - t189 * t4 - t48 * t56 - t49 * t57) * mrSges(7,3) + (t118 * t184 - t158 * t325) * mrSges(6,2) + ((Ifges(6,2) + Ifges(7,3)) * t507 + Ifges(7,6) * t514 + Ifges(7,5) * t515 + t567) * (t272 * t303 + t330 * t419) + (Ifges(7,5) * t48 + Ifges(7,6) * t49) * t492 + (Ifges(6,4) * t490 + Ifges(7,5) * t501 + Ifges(6,2) * t491 + Ifges(6,6) * t485 + Ifges(7,6) * t502 + Ifges(7,3) * t492 + t551) * (t298 * t315 + t303 * t355) + (Ifges(7,4) * t48 + Ifges(7,2) * t49) * t502 + (Ifges(7,4) * t189 + Ifges(7,2) * t188) * t514 + (Ifges(6,5) * t118 + Ifges(6,3) * t325) * t485 + (Ifges(6,4) * t118 + Ifges(6,6) * t325) * t491 + (Ifges(7,1) * t48 + Ifges(7,4) * t49) * t501 + (Ifges(7,1) * t189 + Ifges(7,4) * t188) * t515 + (Ifges(6,1) * t118 + Ifges(6,5) * t325) * t490 + ((t316 * qJD(4) + m(5) * (-qJD(4) * t185 + t92) + m(6) * (-t158 * t401 + t92) - t436) * t304 + (-t303 * t54 - t128 + t470 * t298 - t409 * qJD(4) + (t298 * t153 + t303 * t434) * qJD(5) + m(5) * (-t91 - t403) + m(6) * (t158 * t399 - t303 * t46 - t403)) * t299) * t239 + (t272 * t231 + t243 * t483) * Ifges(4,2) + t547 * t243 / 0.2e1 + (Ifges(5,3) * t272 - t330 * t341) * t488 + (t115 - t164) * (t427 / 0.2e1 - t330 * t374) + (t352 * t330 * t404 + (-t305 * t242 + t300 * t243 + (t272 * t305 - t300 * t330) * qJD(3)) * mrSges(4,3)) * t469 + (-t330 * t376 + t410 / 0.2e1) * t165 + (Ifges(5,6) * t272 - t330 * t345) * t496 + (Ifges(5,5) * t272 - t330 * t349) * t497 + (t334 * t242 - (-qJD(4) * t333 - t91 * t299 + t92 * t304) * t330) * mrSges(5,3) + (-t230 * t330 + t242 * t482) * Ifges(4,1) + (t272 * t230 - t231 * t330 + t242 * t483 + t243 * t482) * Ifges(4,4) + (-t266 / 0.2e1 - t432 / 0.2e1 + t535 + (0.2e1 * pkin(2) * mrSges(3,1) - 0.3e1 / 0.2e1 * t465 + (-0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(3,2)) * t301 + (m(4) * t292 - mrSges(4,1) * t272 - mrSges(4,2) * t330) * pkin(3)) * qJD(1)) * t405 - t330 * t72 * t478 + (Ifges(7,5) * t189 + Ifges(7,6) * t188) * t507 + t395 * (Ifges(4,5) * t242 + Ifges(4,6) * t243) / 0.2e1 - t536 * t171 + t272 * t387 / 0.2e1 + m(7) * (t108 * t4 + t109 * t3 - t12 * t57 + t13 * t56 + t367) + m(6) * t367 + (-t267 / 0.2e1 - t433 / 0.2e1 + (-0.2e1 * pkin(2) * mrSges(3,2) + 0.3e1 / 0.2e1 * t466) * qJD(1)) * qJD(2) * t301 + t12 * t103 + t13 * t104 + t108 * t14 + t109 * t15 + t118 * t117 / 0.2e1 + (-t47 * mrSges(6,3) + Ifges(6,1) * t508 + Ifges(6,4) * t507 + Ifges(6,5) * t495 + t553) * t233 + (-t185 * t243 - t91 * t272 + t278 * t324) * mrSges(5,2) + (-t184 * t243 - t92 * t272 + t278 * t325) * mrSges(5,1) + (mrSges(6,1) * t325 - mrSges(7,1) * t49 + mrSges(7,2) * t48 - mrSges(6,3) * t118) * t157 + t262 * (Ifges(5,5) * t324 - Ifges(5,6) * t325 + Ifges(5,3) * t243) / 0.2e1 + t248 * (Ifges(5,1) * t324 - Ifges(5,4) * t325 + Ifges(5,5) * t243) / 0.2e1 + t247 * (Ifges(5,4) * t324 - Ifges(5,2) * t325 + Ifges(5,6) * t243) / 0.2e1 + t48 * t510 + t49 * t512 + t189 * t517 + t188 * t519 + t47 * (-mrSges(7,1) * t188 + mrSges(7,2) * t189) - (t19 / 0.2e1 - t534 + mrSges(6,1) * t47 + Ifges(6,3) * t495 + Ifges(6,6) * t507 + Ifges(6,5) * t508 + t509) * t422 + t242 * t222 / 0.2e1 + t282 * (-mrSges(4,1) * t243 + mrSges(4,2) * t242) + t292 * (-mrSges(4,1) * t231 + mrSges(4,2) * t230); (t267 * t480 + t306 * t266 / 0.2e1 + ((-Ifges(3,1) * t301 - t465) * t477 + (-Ifges(3,2) * t306 - t466) * t480 - pkin(2) * (mrSges(3,1) * t306 - mrSges(3,2) * t301)) * qJD(1) - t306 * t535 + (-Ifges(3,5) * t301 / 0.2e1 + Ifges(3,6) * t477) * qJD(2)) * qJD(1) + t549 * t104 + m(6) * (t158 * t194 + t258 * t46 + t437) + t470 * t257 - t548 * t103 + ((-t305 * t230 + (qJD(2) * t269 + t231) * t300) * mrSges(4,3) + ((-m(5) * t278 - t407 + (-m(5) * t291 + t532) * qJD(2)) * t300 + (m(5) * t333 + m(6) * t429 - qJD(2) * mrSges(4,2) + t304 * t206 + t409 * t299 + t251) * t305) * qJD(3)) * pkin(3) + (m(6) * t327 + t552) * t290 - m(6) * t120 + t307 + t434 * t195 + t194 * t153 + t203 * t14 - t332 * t15 + t258 * t54 + t291 * t107 + (t203 * t4 - t332 * t3 + t548 * t57 + t549 * t56 - t120 + t437) * m(7) + t536 * t220; t526 * t153 + t531 * t104 + t530 * t103 + t470 * t274 + t307 + ((-mrSges(4,2) * qJD(3) - t251) * t305 + (qJD(3) * t532 + t407 + t440) * t300) * t469 - m(5) * (t185 * t210 - t278 * t391 + t430) - t210 * t206 + t234 * t14 - t331 * t15 + t275 * t54 - t409 * t209 + t434 * t556 + (t234 * t4 - t3 * t331 - t530 * t57 + t531 * t56 + t525) * m(7) + (t158 * t526 + t275 * t46 - t430 + t525) * m(6) + (-m(5) * t362 + t107) * pkin(4) + (0.2e1 * t327 * t520 - t552) * pkin(5); (t181 * t57 + t182 * t56) * mrSges(7,3) + t387 + (t185 * mrSges(5,3) + t506 + t499 + t442 / 0.2e1 - t278 * mrSges(5,1) + t158 * mrSges(6,2) - t157 * mrSges(6,1) - t446 / 0.2e1 - t445 / 0.2e1 - t444 / 0.2e1 + t463 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t247) * t248 + (-t443 / 0.2e1 - t245 / 0.2e1 - t278 * mrSges(5,2) - t184 * mrSges(5,3) + t498) * t247 + (Ifges(7,1) * t182 + Ifges(7,4) * t181) * t500 - t91 * mrSges(5,2) - t92 * mrSges(5,1) - m(7) * (t56 * t65 - t57 * t66) + (-m(6) * t184 - t409) * t185 + (Ifges(7,5) * t182 + Ifges(7,6) * t181) * t493 - t66 * t103 - t65 * t104 + t181 * t513 + t182 * t511 - t157 * (-mrSges(7,1) * t181 + mrSges(7,2) * t182) + (Ifges(7,4) * t182 + Ifges(7,2) * t181) * t503 + t184 * t206 + (-t309 * t298 + (t198 / 0.2e1 + t467 / 0.2e1 - t523 + t537) * t303) * qJD(5) + (t474 / 0.2e1 + t472 / 0.2e1 + t448 / 0.2e1 - t86 * t339 / 0.2e1 - t27 * t347 / 0.2e1 - t28 * t343 / 0.2e1 + t9 * t481 - t302 * t10 / 0.2e1 - t329 * t47 + t356 * mrSges(7,3) + (0.2e1 * (t520 - m(7) / 0.2e1) * t157 - t434) * t184 + (m(7) * t356 + t302 * t14 + t297 * t15) * pkin(6) + t309 * t247 + (t60 * t481 + t59 * t479 + t157 * t351 + t342 * t502 + t346 * t501 + t338 * t492 + t336 * mrSges(7,3) + (m(7) * t336 + t302 * t103 - t297 * t104) * pkin(6)) * qJD(6) + t553) * t298 + (t26 / 0.2e1 + t25 / 0.2e1 + (t518 + Ifges(7,3) / 0.2e1) * t86 + (m(6) * t158 + t153) * t184 + (t157 * mrSges(6,3) + t494 - t467 / 0.2e1 + t322) * t247 + t567) * t303; t310 * t201 + t19 - m(7) * (t110 * t56 - t111 * t57 - t158 * t157) - t534 + t434 * t158 - t110 * t104 - t111 * t103 - t157 * t153 + (t494 + (t518 - Ifges(6,1) / 0.2e1) * t201 + t326 + t523) * t200 + t338 * t507 + t342 * t514 + t346 * t515 + (t157 * t350 + t521 - t537) * qJD(6) + (mrSges(6,1) + t351) * t47 + t357 * mrSges(7,3) + (m(7) * t357 - t297 * t14 + t302 * t15) * pkin(6) + t10 * t481 + t9 * t479; -t157 * (-mrSges(7,1) * t151 + mrSges(7,2) * t150) + (Ifges(7,1) * t150 + t457) * t500 + t59 * t501 + (Ifges(7,5) * t150 + Ifges(7,6) * t151) * t493 - t56 * t103 - t57 * t104 + (t150 * t56 + t151 * t57) * mrSges(7,3) + t358 + t8 + (Ifges(7,2) * t151 + t145 + t60) * t503;];
tauc = t1(:);
