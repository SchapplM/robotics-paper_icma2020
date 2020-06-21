% Calculate time derivative of joint inertia matrix for
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350OL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:00:29
% EndTime: 2020-06-20 08:01:14
% DurationCPUTime: 23.85s
% Computational Cost: add. (9646->728), mult. (23489->1081), div. (0->0), fcn. (22717->10), ass. (0->341)
t314 = sin(qJ(5));
t315 = sin(qJ(4));
t319 = cos(qJ(5));
t418 = qJD(5) * t319;
t320 = cos(qJ(4));
t421 = qJD(4) * t320;
t334 = t314 * t421 + t315 * t418;
t313 = sin(qJ(6));
t318 = cos(qJ(6));
t431 = t318 * t320;
t433 = t315 * t319;
t247 = -t313 * t433 + t431;
t432 = t318 * t319;
t248 = t313 * t320 + t315 * t432;
t470 = Ifges(6,4) * t319;
t349 = -Ifges(6,2) * t314 + t470;
t441 = t314 * t315;
t530 = Ifges(7,5) * t248 + Ifges(6,6) * t320 + Ifges(7,6) * t247 + Ifges(7,3) * t441 - t349 * t315;
t426 = mrSges(6,1) * t320 - mrSges(7,1) * t247 + mrSges(7,2) * t248 + mrSges(6,3) * t433;
t464 = Ifges(7,6) * t313;
t230 = Ifges(7,3) * t319 + (-Ifges(7,5) * t318 + t464) * t314;
t471 = Ifges(6,4) * t314;
t282 = Ifges(6,2) * t319 + t471;
t537 = t282 + t230;
t381 = t418 / 0.2e1;
t384 = t421 / 0.2e1;
t527 = t314 * t384 + t315 * t381;
t536 = 2 * pkin(6);
t316 = sin(qJ(3));
t317 = sin(qJ(2));
t321 = cos(qJ(3));
t322 = cos(qJ(2));
t340 = t316 * t317 - t321 * t322;
t253 = -t316 * t322 - t317 * t321;
t519 = qJD(2) + qJD(3);
t198 = t519 * t253;
t438 = t315 * t198;
t338 = -t340 * t421 + t438;
t423 = qJD(4) * t315;
t451 = t198 * t320;
t337 = -t340 * t423 - t451;
t331 = -qJD(5) * t253 + t337;
t199 = t519 * t340;
t417 = qJD(5) * t320;
t359 = t340 * t417 + t199;
t60 = t359 * t314 - t331 * t319;
t61 = t331 * t314 + t359 * t319;
t18 = Ifges(6,4) * t60 + Ifges(6,2) * t61 + t338 * Ifges(6,6);
t430 = t319 * t320;
t179 = t253 * t314 - t340 * t430;
t330 = qJD(6) * t179 + t338;
t448 = t340 * t315;
t360 = -qJD(6) * t448 - t60;
t25 = t330 * t313 + t360 * t318;
t26 = -t360 * t313 + t330 * t318;
t4 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t61;
t535 = t18 + t4;
t434 = t315 * t318;
t120 = t179 * t313 - t340 * t434;
t121 = -t179 * t318 - t313 * t448;
t439 = t314 * t320;
t178 = t253 * t319 + t340 * t439;
t42 = Ifges(7,5) * t121 + Ifges(7,6) * t120 + Ifges(7,3) * t178;
t89 = Ifges(6,4) * t179 + Ifges(6,2) * t178 - Ifges(6,6) * t448;
t534 = t89 + t42;
t419 = qJD(5) * t315;
t422 = qJD(4) * t319;
t371 = qJD(6) + t422;
t372 = -qJD(6) * t319 - qJD(4);
t420 = qJD(5) * t314;
t157 = t371 * t431 + (t372 * t313 - t318 * t420) * t315;
t397 = t314 * t419;
t158 = t372 * t434 + (-t371 * t320 + t397) * t313;
t69 = Ifges(7,5) * t157 + Ifges(7,6) * t158 + t334 * Ifges(7,3);
t533 = t282 * t419 + (-Ifges(6,6) * t315 - t349 * t320) * qJD(4) + t69;
t335 = t319 * t421 - t397;
t532 = -mrSges(6,1) * t423 - mrSges(7,1) * t158 + mrSges(7,2) * t157 + t335 * mrSges(6,3);
t414 = qJD(6) * t318;
t332 = t313 * t418 + t314 * t414;
t415 = qJD(6) * t314;
t393 = t313 * t415;
t170 = Ifges(7,5) * t393 + (-Ifges(7,5) * t432 - Ifges(7,3) * t314) * qJD(5) + t332 * Ifges(7,6);
t261 = t349 * qJD(5);
t531 = t170 + t261;
t171 = Ifges(6,5) * t397 + (-Ifges(6,5) * t430 - Ifges(6,3) * t315) * qJD(4) + t334 * Ifges(6,6);
t472 = Ifges(5,4) * t320;
t350 = -Ifges(5,2) * t315 + t472;
t529 = t350 * qJD(4) + t171;
t385 = -t421 / 0.2e1;
t388 = t441 / 0.2e1;
t528 = qJD(5) * t388 + t319 * t385;
t526 = t314 * t417 + t315 * t422;
t301 = t317 * pkin(3) + pkin(2);
t189 = -pkin(4) * t253 - pkin(5) * t340 + t301;
t108 = pkin(3) * qJD(2) * t322 - pkin(4) * t199 + pkin(5) * t198;
t456 = t108 * t319;
t520 = t189 * t420 - t456;
t525 = t108 * t441 + t334 * t189;
t524 = 0.2e1 * t314;
t522 = -mrSges(5,1) * t199 - mrSges(6,1) * t61 + mrSges(6,2) * t60 - t337 * mrSges(5,3);
t476 = mrSges(5,2) * t315;
t521 = mrSges(5,1) * t320 + mrSges(4,1) - t476;
t446 = t340 * t320;
t429 = -mrSges(5,1) * t253 - mrSges(6,1) * t178 + mrSges(6,2) * t179 - mrSges(5,3) * t446;
t462 = mrSges(6,1) * t319 - mrSges(6,2) * t314 + mrSges(5,1);
t518 = t462 * t320 - t476;
t221 = mrSges(6,2) * t423 + t334 * mrSges(6,3);
t474 = mrSges(6,2) * t319;
t355 = mrSges(6,1) * t314 + t474;
t250 = t355 * t315;
t267 = -mrSges(6,2) * t320 + mrSges(6,3) * t441;
t453 = t189 * t319;
t517 = qJD(4) * t189 * t250 - t221 * t453 - t267 * t456;
t193 = -mrSges(5,2) * t253 + mrSges(5,3) * t448;
t516 = -qJD(4) * t193 + t522;
t515 = -0.2e1 * pkin(5);
t514 = 2 * m(5);
t513 = 2 * m(6);
t512 = 2 * m(7);
t113 = t334 * mrSges(7,1) - mrSges(7,3) * t157;
t511 = 0.2e1 * t113;
t114 = -t334 * mrSges(7,2) + mrSges(7,3) * t158;
t510 = 0.2e1 * t114;
t216 = -mrSges(7,2) * t441 + mrSges(7,3) * t247;
t509 = 0.2e1 * t216;
t217 = mrSges(7,1) * t441 - mrSges(7,3) * t248;
t508 = 0.2e1 * t217;
t507 = 0.2e1 * t221;
t475 = mrSges(5,2) * t320;
t257 = (-mrSges(5,1) * t315 - t475) * qJD(4);
t506 = 0.2e1 * t257;
t505 = 0.2e1 * t267;
t504 = 0.2e1 * t301;
t502 = -0.2e1 * t319;
t501 = t60 / 0.2e1;
t500 = t61 / 0.2e1;
t469 = Ifges(7,4) * t313;
t281 = Ifges(7,2) * t318 + t469;
t468 = Ifges(7,4) * t318;
t348 = -Ifges(7,2) * t313 + t468;
t172 = t281 * t415 + (-Ifges(7,6) * t314 - t348 * t319) * qJD(5);
t499 = t172 / 0.2e1;
t284 = Ifges(7,1) * t313 + t468;
t351 = Ifges(7,1) * t318 - t469;
t174 = t284 * t415 + (-Ifges(7,5) * t314 - t351 * t319) * qJD(5);
t498 = t174 / 0.2e1;
t497 = t178 / 0.2e1;
t496 = t179 / 0.2e1;
t232 = Ifges(7,6) * t319 - t348 * t314;
t495 = t232 / 0.2e1;
t234 = Ifges(7,5) * t319 - t351 * t314;
t494 = t234 / 0.2e1;
t493 = t247 / 0.2e1;
t492 = t248 / 0.2e1;
t260 = t348 * qJD(6);
t491 = t260 / 0.2e1;
t263 = t351 * qJD(6);
t490 = t263 / 0.2e1;
t489 = t281 / 0.2e1;
t488 = t284 / 0.2e1;
t487 = t313 / 0.2e1;
t486 = t314 / 0.2e1;
t485 = t318 / 0.2e1;
t484 = t319 / 0.2e1;
t483 = t320 / 0.2e1;
t482 = pkin(5) * t320;
t481 = pkin(6) * t114;
t479 = pkin(6) * t319;
t478 = t338 * mrSges(6,1) - mrSges(7,1) * t26 + mrSges(7,2) * t25 - mrSges(6,3) * t60;
t477 = pkin(3) * qJD(3);
t473 = Ifges(5,4) * t315;
t467 = Ifges(5,5) * t199;
t466 = Ifges(5,5) * t253;
t465 = Ifges(5,6) * t253;
t463 = t199 * Ifges(5,6);
t461 = -mrSges(6,1) * t448 - mrSges(7,1) * t120 + mrSges(7,2) * t121 - mrSges(6,3) * t179;
t183 = t189 ^ 2;
t309 = t314 ^ 2;
t400 = t315 * t421;
t310 = t315 ^ 2;
t458 = t108 * t189;
t407 = t310 * t458;
t460 = (t183 * t400 + t407) * t309;
t276 = -mrSges(7,1) * t318 + mrSges(7,2) * t313;
t459 = t276 + mrSges(6,1);
t457 = t108 * t315;
t455 = t108 * t320;
t452 = t189 * t320;
t222 = pkin(4) * t418 + t526 * pkin(5);
t450 = t222 * t319;
t299 = pkin(3) * t316 - pkin(5);
t449 = t250 * t299;
t447 = t340 * t316;
t445 = t299 * t310;
t312 = t320 ^ 2;
t444 = t299 * t312;
t443 = t299 * t320;
t442 = t313 * t314;
t440 = t314 * t318;
t231 = Ifges(6,3) * t320 + (-Ifges(6,5) * t319 + Ifges(6,6) * t314) * t315;
t437 = t315 * t231;
t436 = t315 * t250;
t256 = t355 * qJD(5);
t435 = t315 * t256;
t300 = pkin(3) * t321 + pkin(4);
t410 = t321 * t477;
t374 = t320 * t410;
t411 = t316 * t477;
t129 = (-t299 * t417 - t411) * t319 + (-qJD(5) * t300 + t299 * t423 - t374) * t314;
t223 = -pkin(4) * t420 + (-t314 * t423 + t319 * t417) * pkin(5);
t225 = -t299 * t439 + t300 * t319;
t270 = pkin(4) * t319 + pkin(5) * t439;
t428 = t270 * t129 + t223 * t225;
t427 = Ifges(5,5) * t451 + Ifges(5,3) * t199;
t226 = t299 * t430 + t314 * t300;
t424 = -t310 - t312;
t416 = qJD(6) * t313;
t280 = Ifges(6,5) * t314 + Ifges(6,6) * t319;
t413 = t280 / 0.2e1 - Ifges(5,6);
t412 = 0.2e1 * t317;
t409 = m(7) * t418;
t408 = Ifges(5,5) * t421;
t405 = t189 * t441;
t404 = t267 * t430;
t403 = -pkin(5) - t479;
t394 = t217 * t416;
t390 = t318 * t418;
t389 = t442 / 0.2e1;
t387 = -t433 / 0.2e1;
t386 = -t423 / 0.2e1;
t382 = -t418 / 0.2e1;
t380 = t414 / 0.2e1;
t379 = t299 - t479;
t378 = t424 * mrSges(5,3);
t285 = Ifges(6,1) * t314 + t470;
t352 = Ifges(6,1) * t319 - t471;
t175 = t285 * t419 + (-Ifges(6,5) * t315 - t352 * t320) * qJD(4);
t353 = Ifges(5,1) * t320 - t473;
t265 = t353 * qJD(4);
t377 = -t319 * t175 + t265;
t235 = Ifges(6,5) * t320 - t352 * t315;
t286 = -Ifges(5,1) * t315 - t472;
t376 = -t235 * t319 - t286;
t375 = t129 * t405 + t525 * t225;
t17 = Ifges(6,5) * t60 + Ifges(6,6) * t61 + t338 * Ifges(6,3);
t373 = t223 * t405 + t525 * t270;
t370 = 0.2e1 * t461;
t369 = 0.2e1 * t532;
t368 = 0.2e1 * t426;
t146 = -t340 * t350 + t465;
t88 = Ifges(6,5) * t179 + Ifges(6,6) * t178 - Ifges(6,3) * t448;
t361 = -t146 + t88 - t465;
t27 = pkin(6) * t60 - t189 * t423 + t455;
t356 = pkin(6) * t340 - t453;
t28 = t356 * t421 + (-pkin(6) * t198 + t520) * t315;
t103 = pkin(6) * t179 + t452;
t122 = t356 * t315;
t51 = t103 * t313 - t122 * t318;
t10 = -qJD(6) * t51 + t27 * t318 + t28 * t313;
t50 = t103 * t318 + t122 * t313;
t9 = qJD(6) * t50 + t27 * t313 - t28 * t318;
t358 = -t10 * t313 + t9 * t318;
t311 = t319 ^ 2;
t357 = mrSges(5,2) + (-t309 - t311) * mrSges(6,3);
t354 = mrSges(7,1) * t313 + mrSges(7,2) * t318;
t128 = -t526 * t299 + t300 * t418 - t314 * t411 + t319 * t374;
t306 = pkin(6) * t423;
t119 = t306 + t128;
t213 = -pkin(6) * t320 + t226;
t244 = t379 * t315;
t144 = t213 * t313 + t244 * t318;
t297 = pkin(6) * t397;
t177 = t315 * t410 + t379 * t421 + t297;
t40 = qJD(6) * t144 - t119 * t318 + t177 * t313;
t342 = t213 * t318 - t244 * t313;
t41 = t342 * qJD(6) + t119 * t313 + t177 * t318;
t347 = -t313 * t41 + t318 * t40;
t346 = t313 * t51 + t318 * t50;
t307 = t314 * pkin(4);
t245 = t307 + (-pkin(5) * t319 - pkin(6)) * t320;
t272 = t403 * t315;
t181 = t245 * t313 + t272 * t318;
t204 = t306 + t222;
t224 = t403 * t421 + t297;
t85 = qJD(6) * t181 - t204 * t318 + t224 * t313;
t341 = t245 * t318 - t272 * t313;
t86 = t341 * qJD(6) + t204 * t313 + t224 * t318;
t345 = -t313 * t86 + t318 * t85;
t344 = t144 * t318 - t313 * t342;
t343 = t181 * t318 - t313 * t341;
t43 = Ifges(7,4) * t121 + Ifges(7,2) * t120 + Ifges(7,6) * t178;
t44 = Ifges(7,1) * t121 + Ifges(7,4) * t120 + Ifges(7,5) * t178;
t339 = t44 * t485 - t313 * t43 / 0.2e1;
t333 = t390 - t393;
t132 = mrSges(6,2) * t448 + mrSges(6,3) * t178;
t329 = t132 * t502 + t314 * t370 - 0.2e1 * t193;
t328 = -t313 * t113 - t216 * t416 - t217 * t414;
t167 = Ifges(7,4) * t248 + Ifges(7,2) * t247 + Ifges(7,6) * t441;
t168 = Ifges(7,1) * t248 + Ifges(7,4) * t247 + Ifges(7,5) * t441;
t283 = -Ifges(5,2) * t320 - t473;
t70 = Ifges(7,4) * t157 + Ifges(7,2) * t158 + t334 * Ifges(7,6);
t71 = Ifges(7,1) * t157 + Ifges(7,4) * t158 + t334 * Ifges(7,5);
t327 = t157 * t168 + t158 * t167 + t235 * t397 + t247 * t70 + t248 * t71 + t283 * t423 + t529 * t320 + t530 * t334 + t533 * t441;
t303 = Ifges(7,5) * t414;
t258 = -Ifges(7,6) * t416 + t303;
t279 = Ifges(7,5) * t313 + Ifges(7,6) * t318;
t326 = t157 * t488 + t158 * t489 - t167 * t416 / 0.2e1 + t318 * t481 + t70 * t485 + t71 * t487 + t168 * t380 + t171 + t247 * t491 + t248 * t490 + t258 * t388 + t527 * t279;
t304 = Ifges(6,5) * t418;
t259 = -Ifges(6,6) * t420 + t304;
t264 = t352 * qJD(5);
t302 = Ifges(5,6) * t423;
t325 = t175 * t486 + t174 * t492 + t172 * t493 + t157 * t494 + t158 * t495 + t442 * t481 + t259 * t483 + t302 - t71 * t440 / 0.2e1 + t70 * t389 + t280 * t386 + t264 * t387 + t235 * t381 + t533 * t484 + t531 * t388 - t530 * t420 / 0.2e1 + t528 * t285 + (qJD(6) * t389 + t318 * t382) * t168 + (t313 * t381 + t314 * t380) * t167 + t537 * t527 + (t113 * t440 + t332 * t216 + t217 * t390) * pkin(6);
t147 = -t340 * t353 + t466;
t188 = -t334 * mrSges(6,1) - t335 * mrSges(6,2);
t19 = Ifges(6,1) * t60 + Ifges(6,4) * t61 + t338 * Ifges(6,5);
t5 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t61;
t52 = -t337 * Ifges(5,4) - t338 * Ifges(5,2) + t463;
t53 = -t337 * Ifges(5,1) - t338 * Ifges(5,4) + t467;
t6 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t61;
t90 = Ifges(6,1) * t179 + Ifges(6,4) * t178 - Ifges(6,5) * t448;
t324 = t527 * t534 + t426 * t525 + (t437 / 0.2e1 + Ifges(4,5)) * t198 + t528 * t90 - t529 * t448 / 0.2e1 + t530 * t500 + t235 * t501 + t6 * t492 + t5 * t493 + t175 * t496 + t17 * t483 + t286 * t451 / 0.2e1 + t265 * t446 / 0.2e1 - t283 * t438 / 0.2e1 + t146 * t423 / 0.2e1 + t253 * (t302 - t408) / 0.2e1 + t188 * t452 - t250 * t455 + t147 * t385 + t88 * t386 + t19 * t387 + t532 * t405 + t533 * t497 + t199 * (-Ifges(5,5) * t315 - Ifges(5,6) * t320) / 0.2e1 - t320 * t52 / 0.2e1 - t315 * t53 / 0.2e1 + t9 * t216 + t10 * t217 + Ifges(4,6) * t199 + t26 * t167 / 0.2e1 + t25 * t168 / 0.2e1 + t157 * t44 / 0.2e1 + t158 * t43 / 0.2e1 + t121 * t71 / 0.2e1 + t50 * t113 + t51 * t114 + t120 * t70 / 0.2e1 + t535 * t388 + (-t231 * t384 - t283 * t385 - t286 * t386) * t340 + t189 * t267 * t397;
t271 = -pkin(5) * t430 + t307;
t268 = mrSges(7,1) * t319 + mrSges(7,3) * t440;
t266 = -mrSges(7,2) * t319 + mrSges(7,3) * t442;
t255 = t354 * qJD(6);
t252 = t410 * t445;
t249 = t354 * t314;
t219 = mrSges(7,2) * t420 + t332 * mrSges(7,3);
t218 = -mrSges(7,1) * t420 + t333 * mrSges(7,3);
t187 = -t332 * mrSges(7,1) - t333 * mrSges(7,2);
t180 = t270 * t223;
t98 = t225 * t129;
t97 = -mrSges(5,2) * t199 - t338 * mrSges(5,3);
t82 = mrSges(7,1) * t178 - mrSges(7,3) * t121;
t81 = -mrSges(7,2) * t178 + mrSges(7,3) * t120;
t75 = t338 * mrSges(5,1) - t337 * mrSges(5,2);
t74 = t312 * t458;
t34 = -t338 * mrSges(6,2) + mrSges(6,3) * t61;
t13 = -mrSges(7,2) * t61 + mrSges(7,3) * t26;
t12 = mrSges(7,1) * t61 - mrSges(7,3) * t25;
t1 = [(-mrSges(4,1) * t199 + mrSges(4,2) * t198) * t504 - 0.2e1 * t198 * t340 * Ifges(4,1) + t179 * t19 + t121 * t6 + t120 * t5 + t60 * t90 + 0.2e1 * t9 * t81 + 0.2e1 * t10 * t82 + 0.2e1 * t50 * t12 + 0.2e1 * t51 * t13 + t26 * t43 + t25 * t44 + t534 * t61 + t535 * t178 + (((2 * Ifges(4,2)) + Ifges(5,3)) * t199 + t427) * t253 + (t183 * t310 * t314 * t418 + t10 * t50 + t51 * t9 + t460) * t512 + (t74 + t407) * t514 + (t311 * t407 + t460 + t74) * t513 + 0.2e1 * (t198 * t253 - t199 * t340) * Ifges(4,4) + ((-pkin(2) * mrSges(3,2) + Ifges(3,4) * t317) * t412 + (m(4) * pkin(3) * t504 + 0.2e1 * pkin(3) * (-mrSges(4,1) * t253 - mrSges(4,2) * t340) + 0.2e1 * pkin(2) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t322 + (-Ifges(3,1) + Ifges(3,2)) * t412) * t322) * qJD(2) + (t198 * t147 - (t53 + t467) * t340 + (t329 * t189 - t340 * t361) * qJD(4) + 0.2e1 * t522 * t189 + 0.2e1 * t429 * t108) * t320 + ((t311 - 0.1e1) * t183 * t421 * t513 + t361 * t198 - (-t463 + t17 - t52 + (-t147 - t466) * qJD(4)) * t340 + t329 * t108 + (t34 * t502 - 0.2e1 * t97 + (t132 * t524 + t319 * t370) * qJD(5) - 0.2e1 * t429 * qJD(4) + t478 * t524) * t189) * t315; t324 + t478 * t225 + (m(6) * (t108 * t443 - t128 * t453 + t189 * t374 + t226 * t520) + mrSges(5,1) * t340 * t411 + t429 * t410 + t516 * t299 + t517) * t315 + (t429 * t443 + (-t404 + m(6) * (-t226 * t430 + t444 - t445)) * t189) * qJD(4) + t97 * t443 + (-Ifges(3,5) * t317 - Ifges(3,6) * t322) * qJD(2) + t300 * t75 + t226 * t34 + t128 * t132 + t144 * t12 - t342 * t13 + t40 * t81 + t41 * t82 + ((mrSges(5,2) * t447 + t193 * t321) * t320 * qJD(3) + (-t198 * t321 + t199 * t316 + (t253 * t321 - t447) * qJD(3)) * mrSges(4,3)) * pkin(3) + m(6) * t375 + m(7) * (t10 * t144 - t342 * t9 + t40 * t51 + t41 * t50 + t375) + t461 * t129; t327 + t300 * t506 + t128 * t505 + t226 * t507 + t40 * t509 + t41 * t508 + t144 * t511 - t342 * t510 + (-t437 + (t376 - 0.2e1 * t449) * t320) * qJD(4) + (0.2e1 * t188 * t299 + t377) * t315 + t252 * t514 + t225 * t369 + t129 * t368 + ((t444 * t514 - 0.2e1 * mrSges(4,2) + 0.2e1 * t378 - 0.2e1 * t436) * t321 + (-t300 * t514 - 0.2e1 * t521) * t316) * t477 + (t299 ^ 2 * t400 + t128 * t226 + t252 + t98) * t513 + (t144 * t41 - t342 * t40 + t98) * t512; t324 + m(6) * t373 + m(7) * (t10 * t181 - t341 * t9 + t50 * t86 + t51 * t85 + t373) + (-t429 * t482 + (-t404 + m(6) * (-t271 * t430 + (t310 - t312) * pkin(5))) * t189) * qJD(4) + (m(6) * (-t189 * t450 + t520 * t271) + (-m(6) * t455 - t516) * pkin(5) + t517) * t315 + t478 * t270 + t461 * t223 - t97 * t482 + t271 * t34 + t222 * t132 + t181 * t12 - t341 * t13 + t85 * t81 + t86 * t82 + pkin(4) * t75; t327 + (-t437 + (-t449 + (-0.2e1 * m(6) * t299 * t315 + t250) * pkin(5) + t376) * t320) * qJD(4) + (t222 + t128) * t267 + (t300 + pkin(4)) * t257 + (t271 + t226) * t221 + (t86 + t41) * t217 + (t85 + t40) * t216 + (-t341 - t342) * t114 + (t181 + t144) * t113 + m(7) * (t144 * t86 + t181 * t41 - t341 * t40 - t342 * t85 + t428) + ((-pkin(5) + t299) * t188 + t377) * t315 + (-t521 * t316 + (-m(6) * t310 * pkin(5) - mrSges(4,2) + t378 - t436) * t321 + (t424 * pkin(5) * t321 - pkin(4) * t316) * m(5)) * t477 + m(6) * (t128 * t271 + t222 * t226 + t428) + (t223 + t129) * t426 + t532 * (t270 + t225); t327 + (t188 * t515 + t377) * t315 + (pkin(5) ^ 2 * t400 + t222 * t271 + t180) * t513 + (t181 * t86 - t341 * t85 + t180) * t512 + (-t437 + (-t250 * t515 + t376) * t320) * qJD(4) + t222 * t505 + t271 * t507 + pkin(4) * t506 + t85 * t509 + t86 * t508 + t181 * t511 - t341 * t510 + t270 * t369 + t223 * t368; t285 * t501 + t264 * t496 + t9 * t266 + t10 * t268 + t26 * t495 + t25 * t494 + t51 * t219 + t50 * t218 + t121 * t498 + t120 * t499 + (t282 / 0.2e1 + t230 / 0.2e1) * t61 + (t261 / 0.2e1 + t170 / 0.2e1) * t178 + (t189 * t256 - t462 * t108 + (t357 * t189 - t340 * t413) * qJD(4)) * t320 + (t18 / 0.2e1 + t4 / 0.2e1 + (t90 / 0.2e1 + (m(7) * t346 + t313 * t81 + t318 * t82) * pkin(6) - t339) * qJD(5)) * t319 + (-t340 * t259 / 0.2e1 - t189 * t249 * t418 + t413 * t198 + t357 * t108 + (Ifges(5,5) * t340 + t462 * t189) * qJD(4)) * t315 + (-t318 * t6 / 0.2e1 + t5 * t487 - t249 * t457 + t19 / 0.2e1 + (t187 * t315 - t249 * t421) * t189 + (t43 * t485 + t44 * t487) * qJD(6) + (-t89 / 0.2e1 - t42 / 0.2e1) * qJD(5) + (m(7) * (t10 * t318 + t313 * t9 + t51 * t414 - t50 * t416) + t81 * t414 - t82 * t416 + t318 * t12 + t313 * t13) * pkin(6)) * t314 + t427; t325 + (t344 * t409 + (m(7) * (-t144 * t416 + t313 * t40 + t318 * t41 - t342 * t414) - t394) * t314) * pkin(6) + (-t462 * t315 - t475) * t410 + t299 * t435 + (t128 * t319 - t129 * t314 + (-t225 * t319 - t226 * t314) * qJD(5)) * mrSges(6,3) + (-Ifges(5,5) * t320 - t299 * t518) * qJD(4) + t40 * t266 + t41 * t268 - t129 * t249 - t342 * t219 + t225 * t187 + t144 * t218; t325 + (t343 * t409 + (m(7) * (-t181 * t416 + t313 * t85 + t318 * t86 - t341 * t414) - t394) * t314) * pkin(6) - t408 + t85 * t266 + t86 * t268 + t270 * t187 - t223 * t249 - t341 * t219 + t181 * t218 + (qJD(4) * t518 - t435) * pkin(5) + (t450 - t223 * t314 + (-t270 * t319 - t271 * t314) * qJD(5)) * mrSges(6,3); ((t313 * t232 - t318 * t234 + t285 + (t266 * t313 + t268 * t318) * t536) * qJD(5) + t531) * t319 + (t313 * t172 - t318 * t174 + t264 + (t232 * t318 + t234 * t313) * qJD(6) + ((t313 ^ 2 + t318 ^ 2) * (pkin(6) ^ 2) * t319 * t512 - t537) * qJD(5) + (t318 * t218 + t313 * t219 + (t266 * t318 - t268 * t313) * qJD(6)) * t536) * t314; t5 * t485 + t6 * t487 + t279 * t500 + t26 * t489 + t25 * t488 + t120 * t491 + t121 * t490 + t258 * t497 + t339 * qJD(6) + (t459 * t314 + t474) * t457 + (-t346 * qJD(6) + t358) * mrSges(7,3) + ((mrSges(6,2) * t421 + t459 * t419) * t319 + ((-mrSges(6,2) * qJD(5) + t255) * t315 + t459 * t421) * t314) * t189 + (m(7) * (-t50 * t414 - t51 * t416 + t358) + t318 * t13 - t313 * t12 - t82 * t414 - t81 * t416) * pkin(6) + t17; (m(7) * (-t144 * t414 + t342 * t416 + t347) + t328) * pkin(6) + (-t344 * qJD(6) + t347) * mrSges(7,3) + t459 * t129 + t326 + t225 * t255 - t128 * mrSges(6,2); (m(7) * (-t181 * t414 + t341 * t416 + t345) + t328) * pkin(6) + (-t343 * qJD(6) + t345) * mrSges(7,3) + t459 * t223 + t326 + t270 * t255 - t222 * mrSges(6,2); t258 * t484 + t304 + (-t279 / 0.2e1 - Ifges(6,6)) * t420 + (pkin(6) * t219 + t499 + t284 * t382 - t314 * t263 / 0.2e1 + (-pkin(6) * t268 + t281 * t486 + t494) * qJD(6)) * t318 + (-pkin(6) * t218 + t498 + t281 * t381 + t260 * t486 + (-pkin(6) * t266 - t232 / 0.2e1 + t284 * t486) * qJD(6)) * t313; t260 * t318 + t263 * t313 + (-t281 * t313 + t284 * t318) * qJD(6); mrSges(7,1) * t10 - mrSges(7,2) * t9 + t4; mrSges(7,1) * t41 - mrSges(7,2) * t40 + t69; mrSges(7,1) * t86 - mrSges(7,2) * t85 + t69; (t333 * mrSges(7,1) - t332 * mrSges(7,2)) * pkin(6) + t170; t303 + (t276 * pkin(6) - t464) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
