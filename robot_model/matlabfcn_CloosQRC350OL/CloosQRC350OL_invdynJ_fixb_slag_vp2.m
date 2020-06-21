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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 08:00:28
% EndTime: 2020-06-20 08:02:04
% DurationCPUTime: 57.25s
% Computational Cost: add. (20316->920), mult. (43264->1311), div. (0->0), fcn. (33426->12), ass. (0->392)
t338 = sin(qJ(3));
t339 = sin(qJ(2));
t343 = cos(qJ(3));
t344 = cos(qJ(2));
t297 = -t338 * t344 - t343 * t339;
t288 = t297 * qJD(1);
t336 = sin(qJ(5));
t341 = cos(qJ(5));
t342 = cos(qJ(4));
t451 = t341 * t342;
t447 = qJD(1) * t344;
t455 = t338 * t339;
t619 = qJD(1) * t455 - t343 * t447;
t225 = t288 * t451 + t336 * t619;
t337 = sin(qJ(4));
t441 = qJD(5) * t337;
t420 = t336 * t441;
t634 = t225 - t420;
t633 = t634 * pkin(6);
t244 = -pkin(4) * t619 + t288 * pkin(5);
t526 = pkin(3) * qJD(2);
t432 = t343 * t526;
t219 = -t244 * t337 + t342 * t432;
t433 = t338 * t526;
t440 = qJD(5) * t341;
t439 = qJD(5) * t342;
t444 = qJD(4) * t341;
t618 = t336 * t439 + t337 * t444;
t604 = pkin(4) * t440 + pkin(5) * t618 - t341 * t219 + t336 * t433;
t436 = qJD(2) + qJD(3);
t255 = -t337 * t436 - t342 * t619;
t282 = qJD(4) + t288;
t205 = -t255 * t336 + t341 * t282;
t254 = t337 * t619 - t342 * t436;
t253 = qJD(5) - t254;
t504 = Ifges(6,6) * t253;
t206 = t255 * t341 + t282 * t336;
t516 = Ifges(6,4) * t206;
t119 = Ifges(6,2) * t205 + t504 + t516;
t203 = qJD(6) + t205;
t502 = Ifges(7,3) * t203;
t335 = sin(qJ(6));
t340 = cos(qJ(6));
t160 = t206 * t335 + t253 * t340;
t503 = Ifges(7,6) * t160;
t161 = t206 * t340 - t335 * t253;
t507 = Ifges(7,5) * t161;
t63 = t502 + t503 - t507;
t625 = t119 + t63;
t218 = t244 * t342 + t337 * t432;
t537 = pkin(6) * t341;
t428 = -pkin(5) - t537;
t443 = qJD(4) * t342;
t632 = -t428 * t443 + t218 + t633;
t229 = pkin(3) * t447 + t244;
t320 = pkin(3) * t338 - pkin(5);
t405 = t320 - t537;
t525 = pkin(3) * qJD(3);
t434 = t343 * t525;
t631 = t229 * t342 - t337 * t434 - t405 * t443 + t633;
t445 = qJD(4) * t337;
t326 = pkin(6) * t445;
t464 = t288 * t337;
t630 = pkin(6) * t464 + t326 + t604;
t321 = pkin(3) * t343 + pkin(4);
t401 = t342 * t434;
t435 = t338 * t525;
t197 = -t320 * t618 + t321 * t440 - t336 * t435 + t341 * t401;
t629 = -(-pkin(6) * t288 - t229 * t341) * t337 + t326 + t197;
t442 = qJD(5) * t336;
t628 = t219 * t336 + t341 * t433 - pkin(4) * t442 + (-t336 * t445 + t341 * t439) * pkin(5);
t331 = t339 * pkin(3);
t529 = pkin(2) + t331;
t309 = t529 * qJD(1);
t350 = -pkin(4) * t288 - pkin(5) * t619 + t309;
t220 = t337 * t350;
t354 = -pkin(5) * t436 + t433;
t187 = t342 * t354 - t220;
t186 = t337 * t354 + t342 * t350;
t437 = qJD(1) * qJD(2);
t301 = qJDD(1) * t344 - t339 * t437;
t302 = -t339 * qJDD(1) - t344 * t437;
t361 = t297 * qJD(3);
t213 = qJD(1) * t361 + t343 * t301 + t338 * t302;
t214 = qJD(3) * t619 - t338 * t301 + t343 * t302;
t279 = qJDD(1) * pkin(2) - pkin(3) * t302;
t127 = -pkin(4) * t214 + pkin(5) * t213 + t279;
t431 = qJD(2) * t525;
t474 = pkin(3) * qJDD(2);
t295 = t338 * t474 + t343 * t431;
t332 = qJDD(2) + qJDD(3);
t277 = -pkin(5) * t332 + t295;
t73 = t337 * t277 - qJD(4) * t220 + (qJD(4) * t354 + t127) * t342;
t481 = t337 * t73;
t365 = t186 * t443 + t481;
t627 = -t187 * t445 + t365;
t154 = -qJD(4) * t255 - t337 * t213 - t342 * t332;
t151 = qJDD(5) - t154;
t153 = qJD(4) * t254 + t342 * t213 - t337 * t332;
t210 = qJDD(4) + t214;
t80 = qJD(5) * t205 + t153 * t341 + t210 * t336;
t29 = qJD(6) * t160 + t151 * t335 - t340 * t80;
t30 = qJD(6) * t161 + t151 * t340 + t335 * t80;
t81 = -qJD(5) * t206 - t153 * t336 + t341 * t210;
t77 = qJDD(6) + t81;
t8 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t77;
t613 = t80 * Ifges(6,4) + t81 * Ifges(6,2) + t151 * Ifges(6,6) + t8;
t626 = t613 / 0.2e1;
t457 = t337 * t340;
t179 = t225 * t335 + t288 * t457;
t399 = qJD(6) + t444;
t400 = -qJD(6) * t341 - qJD(4);
t222 = t400 * t457 + (-t342 * t399 + t420) * t335;
t624 = t179 - t222;
t461 = t335 * t337;
t180 = -t225 * t340 + t288 * t461;
t453 = t340 * t342;
t221 = t399 * t453 + (t335 * t400 - t340 * t442) * t337;
t623 = t180 - t221;
t458 = t336 * t342;
t224 = -t288 * t458 + t341 * t619;
t424 = t336 * t443;
t357 = t337 * t440 + t424;
t622 = t224 - t357;
t621 = t341 * t443 + t634;
t520 = Ifges(4,4) * t619;
t620 = t255 * Ifges(5,5) + Ifges(4,2) * t288 + Ifges(4,6) * t436 + t254 * Ifges(5,6) + t282 * Ifges(5,3) - t520;
t617 = t445 + t464;
t305 = pkin(4) * t436 + t432;
t167 = t187 * t341 + t305 * t336;
t162 = -mrSges(6,2) * t253 + mrSges(6,3) * t205;
t216 = -mrSges(5,2) * t282 + mrSges(5,3) * t254;
t403 = -t341 * t162 - t216;
t476 = mrSges(6,1) * t253 - mrSges(7,1) * t160 - mrSges(7,2) * t161 - mrSges(6,3) * t206;
t353 = t336 * t476 + t403;
t456 = t337 * t341;
t470 = t186 * t342;
t450 = -mrSges(5,1) * t282 - mrSges(6,1) * t205 + mrSges(6,2) * t206 + mrSges(5,3) * t255;
t605 = t342 * t450;
t616 = -m(6) * (-t167 * t456 + t470) - m(5) * (-t187 * t337 + t470) - t605 - t337 * t353;
t497 = t167 * mrSges(6,2);
t615 = t305 * mrSges(5,1) - t187 * mrSges(5,3) - t497;
t333 = qJ(2) + qJ(3);
t328 = sin(t333);
t541 = pkin(4) * t328;
t266 = t320 * t451 + t336 * t321;
t256 = -pkin(6) * t342 + t266;
t283 = t405 * t337;
t211 = t256 * t335 + t283 * t340;
t612 = -qJD(6) * t211 + t631 * t335 + t629 * t340;
t370 = t256 * t340 - t283 * t335;
t611 = qJD(6) * t370 + t629 * t335 - t631 * t340;
t330 = t336 * pkin(4);
t285 = t330 + (-pkin(5) * t341 - pkin(6)) * t342;
t303 = t428 * t337;
t369 = t285 * t340 - t303 * t335;
t610 = qJD(6) * t369 + t630 * t335 - t632 * t340;
t294 = -t338 * t431 + t343 * t474;
t276 = pkin(4) * t332 + t294;
t609 = m(5) * t276 - mrSges(5,1) * t154 + mrSges(5,2) * t153;
t240 = t285 * t335 + t303 * t340;
t608 = -qJD(6) * t240 + t632 * t335 + t630 * t340;
t391 = mrSges(5,1) * t337 + mrSges(5,2) * t342;
t607 = t305 * t391;
t329 = cos(t333);
t396 = -pkin(5) * t329 - t541;
t603 = -t331 + t396;
t438 = qJD(6) * t336;
t355 = t335 * t440 + t340 * t438;
t446 = qJD(4) * t186;
t72 = -t337 * t127 + t342 * t277 - t446;
t37 = -t187 * t442 + t336 * t276 + t305 * t440 + t341 * t72;
t601 = t167 * t442 - t341 * t37;
t166 = -t187 * t336 + t305 * t341;
t299 = pkin(4) * t341 + pkin(5) * t458;
t38 = -qJD(5) * t167 + t276 * t341 - t336 * t72;
t600 = t628 * t166 + t38 * t299;
t198 = (-t320 * t439 - t435) * t341 + (-qJD(5) * t321 + t320 * t445 - t401) * t336;
t265 = -t320 * t458 + t321 * t341;
t460 = t336 * t337;
t430 = t166 * t460;
t599 = t166 * t198 - t229 * t430 + t38 * t265;
t597 = -t305 * mrSges(5,2) - t186 * mrSges(5,3);
t272 = t328 * t458 - t329 * t341;
t273 = -t328 * t451 - t329 * t336;
t523 = mrSges(5,2) * t337;
t596 = (mrSges(4,2) + mrSges(5,3)) * t329 - (t273 * t335 - t328 * t457) * mrSges(7,2) - (-t273 * t340 - t328 * t461) * mrSges(7,1) - t273 * mrSges(6,1) - t328 * t523 + (-mrSges(7,3) - mrSges(6,2)) * t272;
t521 = Ifges(3,4) * t344;
t522 = Ifges(3,4) * t339;
t542 = t344 / 0.2e1;
t595 = -((-Ifges(3,1) * t339 - t521) * t542 - t339 * (-Ifges(3,2) * t344 - t522) / 0.2e1 + pkin(2) * (mrSges(3,1) * t344 - mrSges(3,2) * t339)) * qJD(1) + t339 * (Ifges(3,5) * qJD(2) + (t344 * Ifges(3,1) - t522) * qJD(1)) / 0.2e1 + (Ifges(3,6) * qJD(2) + (-Ifges(3,2) * t339 + t521) * qJD(1)) * t542;
t106 = -mrSges(5,2) * t210 + mrSges(5,3) * t154;
t594 = t450 * qJD(4) + t106;
t593 = -t73 * mrSges(5,1) - t72 * mrSges(5,2) + Ifges(5,5) * t153 + Ifges(5,6) * t154 + Ifges(5,3) * t210;
t581 = -t210 * Ifges(5,6) / 0.2e1 - t154 * Ifges(5,2) / 0.2e1 - t153 * Ifges(5,4) / 0.2e1;
t131 = -pkin(6) * t253 + t167;
t346 = t206 * pkin(6) + t186;
t61 = t335 * t131 + t340 * t346;
t62 = t340 * t131 - t335 * t346;
t372 = t335 * t62 - t340 * t61;
t374 = Ifges(7,5) * t340 - Ifges(7,6) * t335;
t511 = Ifges(7,4) * t340;
t379 = -Ifges(7,2) * t335 + t511;
t512 = Ifges(7,4) * t335;
t384 = Ifges(7,1) * t340 - t512;
t513 = Ifges(7,4) * t161;
t64 = Ifges(7,2) * t160 + Ifges(7,6) * t203 - t513;
t483 = t335 * t64;
t544 = t340 / 0.2e1;
t561 = t203 / 0.2e1;
t566 = -t161 / 0.2e1;
t567 = t160 / 0.2e1;
t155 = Ifges(7,4) * t160;
t65 = -Ifges(7,1) * t161 + Ifges(7,5) * t203 + t155;
t591 = t374 * t561 + t379 * t567 + t384 * t566 + t372 * mrSges(7,3) - t483 / 0.2e1 + t65 * t544;
t590 = m(5) / 0.2e1;
t589 = m(6) / 0.2e1;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t30 + Ifges(7,6) * t77;
t588 = t9 / 0.2e1;
t587 = pkin(6) * m(7);
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t30 + Ifges(7,5) * t77;
t586 = t10 / 0.2e1;
t148 = Ifges(6,3) * t151;
t75 = Ifges(6,6) * t81;
t76 = Ifges(6,5) * t80;
t19 = t76 + t75 + t148;
t585 = t19 / 0.2e1;
t21 = t80 * Ifges(6,1) + t81 * Ifges(6,4) + t151 * Ifges(6,5);
t584 = t21 / 0.2e1;
t583 = t29 / 0.2e1;
t582 = t30 / 0.2e1;
t580 = -t64 / 0.2e1;
t579 = t64 / 0.2e1;
t578 = -t65 / 0.2e1;
t577 = t65 / 0.2e1;
t576 = t77 / 0.2e1;
t575 = t80 / 0.2e1;
t574 = t81 / 0.2e1;
t573 = -m(6) - m(5);
t490 = t253 * Ifges(6,3);
t493 = t206 * Ifges(6,5);
t494 = t205 * Ifges(6,6);
t118 = t490 + t493 + t494;
t572 = t118 / 0.2e1;
t202 = Ifges(6,4) * t205;
t508 = Ifges(6,5) * t253;
t120 = Ifges(6,1) * t206 + t202 + t508;
t571 = -t120 / 0.2e1;
t569 = t151 / 0.2e1;
t568 = -t160 / 0.2e1;
t565 = t161 / 0.2e1;
t486 = t282 * Ifges(5,6);
t489 = t254 * Ifges(5,2);
t519 = Ifges(5,4) * t255;
t171 = t486 + t489 + t519;
t564 = -t171 / 0.2e1;
t252 = Ifges(5,4) * t254;
t487 = t282 * Ifges(5,5);
t488 = t255 * Ifges(5,1);
t172 = t252 + t487 + t488;
t563 = -t172 / 0.2e1;
t562 = -t203 / 0.2e1;
t560 = -t205 / 0.2e1;
t559 = t205 / 0.2e1;
t558 = -t206 / 0.2e1;
t557 = t206 / 0.2e1;
t553 = -t253 / 0.2e1;
t552 = t253 / 0.2e1;
t551 = -t254 / 0.2e1;
t550 = -t255 / 0.2e1;
t549 = t255 / 0.2e1;
t548 = -t282 / 0.2e1;
t546 = -t619 / 0.2e1;
t545 = -t337 / 0.2e1;
t108 = mrSges(7,1) * t203 + mrSges(7,3) * t161;
t539 = pkin(6) * t108;
t538 = pkin(6) * t205;
t535 = g(3) * t328;
t534 = g(3) * t329;
t533 = t37 * mrSges(6,2);
t532 = t72 * mrSges(5,3);
t527 = mrSges(6,1) * t151 - mrSges(7,1) * t30 + mrSges(7,2) * t29 - mrSges(6,3) * t80;
t524 = mrSges(5,1) * t342;
t518 = Ifges(5,4) * t337;
t517 = Ifges(5,4) * t342;
t515 = Ifges(6,4) * t336;
t514 = Ifges(6,4) * t341;
t510 = Ifges(3,5) * t339;
t248 = qJD(2) * t297 + t361;
t509 = Ifges(4,5) * t248;
t506 = Ifges(3,6) * t344;
t368 = -t343 * t344 + t455;
t249 = t436 * t368;
t505 = Ifges(4,6) * t249;
t501 = t153 * Ifges(5,1);
t499 = t154 * Ifges(5,4);
t492 = t210 * Ifges(5,5);
t485 = t288 * mrSges(4,3);
t484 = t619 * mrSges(4,3);
t482 = t336 * t38;
t479 = t342 * t72;
t475 = -mrSges(5,1) * t210 - mrSges(6,1) * t81 + mrSges(6,2) * t80 + mrSges(5,3) * t153;
t472 = t186 * t218;
t471 = t186 * t336;
t469 = t248 * t337;
t468 = t248 * t342;
t467 = t254 * t336;
t466 = t254 * t341;
t463 = t288 * t342;
t462 = t335 * t336;
t459 = t336 * t340;
t454 = t340 * t341;
t449 = -mrSges(4,1) * t436 + mrSges(5,1) * t254 - mrSges(5,2) * t255 - t484;
t429 = t186 * t337 * t343;
t417 = t340 * t440;
t416 = t335 * t438;
t413 = t118 * t545;
t412 = t171 * t337 / 0.2e1;
t407 = t440 / 0.2e1;
t404 = m(4) * t309 - mrSges(4,1) * t288 - mrSges(4,2) * t619;
t175 = -t249 * pkin(4) + t248 * pkin(5) + t344 * t526;
t245 = -pkin(4) * t297 - pkin(5) * t368 + t529;
t402 = t175 * t430 + (t166 * t357 + t38 * t460) * t245;
t394 = t368 * t439 + t249;
t352 = -qJD(5) * t297 - t368 * t445 - t468;
t121 = t336 * t394 - t341 * t352;
t393 = -qJD(6) * t337 * t368 - t121;
t392 = pkin(6) * t368 - t245 * t341;
t390 = -mrSges(6,1) * t341 + mrSges(6,2) * t336;
t389 = mrSges(6,1) * t336 + mrSges(6,2) * t341;
t388 = mrSges(7,1) * t335 + mrSges(7,2) * t340;
t387 = Ifges(5,1) * t342 - t518;
t386 = Ifges(6,1) * t341 - t515;
t385 = Ifges(6,1) * t336 + t514;
t383 = Ifges(7,1) * t335 + t511;
t382 = -Ifges(5,2) * t337 + t517;
t381 = -Ifges(6,2) * t336 + t514;
t380 = Ifges(6,2) * t341 + t515;
t378 = Ifges(7,2) * t340 + t512;
t377 = Ifges(5,5) * t342 - Ifges(5,6) * t337;
t376 = Ifges(6,5) * t341 - Ifges(6,6) * t336;
t375 = Ifges(6,5) * t336 + Ifges(6,6) * t341;
t373 = Ifges(7,5) * t335 + Ifges(7,6) * t340;
t239 = t297 * t336 - t368 * t451;
t173 = pkin(6) * t239 + t245 * t342;
t195 = t392 * t337;
t111 = t173 * t340 + t195 * t335;
t112 = t173 * t335 - t195 * t340;
t367 = mrSges(7,1) * t340 - mrSges(7,2) * t335 - mrSges(6,1);
t291 = t335 * t342 + t337 * t454;
t290 = -t335 * t456 + t453;
t362 = t186 * t389;
t275 = -t328 * t336 + t329 * t451;
t24 = -pkin(6) * t151 + t37;
t25 = t80 * pkin(6) + t73;
t3 = t61 * qJD(6) - t340 * t24 + t335 * t25;
t4 = t62 * qJD(6) + t335 * t24 + t340 * t25;
t360 = g(3) * t275 + t3 * t340 - t335 * t4;
t356 = -t416 + t417;
t351 = qJD(6) * t239 - t368 * t443 + t469;
t349 = t479 + t627;
t280 = Ifges(4,4) * t288;
t231 = -Ifges(4,1) * t619 + Ifges(4,5) * t436 + t280;
t60 = t492 + t499 + t501;
t345 = (t375 * t552 + t380 * t559 + t385 * t557) * t441 + (-t186 * t463 + t187 * t464 - t627) * mrSges(5,3) - t288 * t607 + (-t607 + (-Ifges(6,3) * t337 - t342 * t376) * t552 + (-Ifges(6,5) * t337 - t342 * t386) * t557 + (-Ifges(6,6) * t337 - t342 * t381) * t559 + t412 + t413) * qJD(4) + t620 * t546 + (mrSges(6,2) * t617 - mrSges(6,3) * t622) * t167 + (-mrSges(7,1) * t622 + mrSges(7,3) * t623) * t61 + (-mrSges(7,2) * t622 + mrSges(7,3) * t624) * t62 + (-mrSges(6,1) * t617 + mrSges(7,1) * t624 - mrSges(7,2) * t623 + mrSges(6,3) * t621) * t166 + t625 * (t337 * t407 + t424 / 0.2e1 - t224 / 0.2e1) + (mrSges(5,1) * t619 + mrSges(6,1) * t622 - mrSges(6,2) * t621) * t186 - t309 * (-mrSges(4,1) * t619 + mrSges(4,2) * t288) - (Ifges(4,2) * t619 + t231 + t280) * t288 / 0.2e1 + (Ifges(5,3) * t619 + t288 * t377) * t548 + (Ifges(5,5) * t619 + t288 * t387) * t550 + (Ifges(5,6) * t619 + t288 * t382) * t551 - t436 * (Ifges(4,5) * t288 + Ifges(4,6) * t619) / 0.2e1 + t619 * (Ifges(4,1) * t288 + t520) / 0.2e1 + (mrSges(6,1) * t342 - mrSges(7,1) * t290 + mrSges(7,2) * t291 + mrSges(6,3) * t456) * t38 + t187 * mrSges(5,2) * t619 - t21 * t456 / 0.2e1 - t389 * t481 - (t120 * t341 + t172) * t443 / 0.2e1 - (t254 * t382 + t255 * t387 + t282 * t377) * qJD(4) / 0.2e1 + t3 * (-mrSges(7,2) * t460 + mrSges(7,3) * t290) + t37 * (-mrSges(6,2) * t342 + mrSges(6,3) * t460) + t4 * (mrSges(7,1) * t460 - mrSges(7,3) * t291) + t291 * t586 + t290 * t588 + t342 * t585 + (Ifges(6,6) * t342 - t337 * t381) * t574 + (Ifges(6,5) * t342 - t337 * t386) * t575 + (Ifges(7,5) * t291 + Ifges(7,6) * t290 + Ifges(7,3) * t460) * t576 + t221 * t577 + t180 * t578 + t222 * t579 + t179 * t580 + t342 * t581 + (Ifges(7,4) * t291 + Ifges(7,2) * t290 + Ifges(7,6) * t460) * t582 + (Ifges(7,1) * t291 + Ifges(7,4) * t290 + Ifges(7,5) * t460) * t583 + t463 * t563 + (Ifges(7,1) * t180 + Ifges(7,4) * t179 + Ifges(7,5) * t224) * t565 + (Ifges(7,1) * t221 + Ifges(7,4) * t222 + Ifges(7,5) * t357) * t566 + (Ifges(7,4) * t221 + Ifges(7,2) * t222 + Ifges(7,6) * t357) * t567 + (Ifges(7,4) * t180 + Ifges(7,2) * t179 + Ifges(7,6) * t224) * t568 + (Ifges(6,3) * t342 - t337 * t376) * t569 + t225 * t571 + (Ifges(6,5) * t225 + Ifges(6,6) * t224 + Ifges(6,3) * t464) * t553 + (Ifges(6,1) * t225 + Ifges(6,4) * t224 + Ifges(6,5) * t464) * t558 + (Ifges(6,4) * t225 + Ifges(6,2) * t224 + Ifges(6,6) * t464) * t560 + (Ifges(7,5) * t221 + Ifges(7,6) * t222 + Ifges(7,3) * t357) * t561 + (Ifges(7,5) * t180 + Ifges(7,6) * t179 + Ifges(7,3) * t224) * t562 + t60 * t545 + t153 * (-Ifges(5,1) * t337 - t517) / 0.2e1 + t154 * (-Ifges(5,2) * t342 - t518) / 0.2e1 + t294 * mrSges(4,1) - t295 * mrSges(4,2) + t288 * t412 + t288 * t413 + t276 * (-t523 + t524) + Ifges(4,6) * t214 + Ifges(4,5) * t213 + t432 * t485 + Ifges(4,3) * t332 + (qJD(5) * t120 + t613) * t460 / 0.2e1 + t210 * (-Ifges(5,5) * t337 - Ifges(5,6) * t342) / 0.2e1;
t308 = -t339 * mrSges(3,1) - mrSges(3,2) * t344;
t300 = -pkin(5) * t451 + t330;
t270 = t272 * pkin(6);
t258 = -mrSges(4,2) * t436 + t485;
t238 = t297 * t341 + t368 * t458;
t192 = -t239 * t340 - t368 * t461;
t191 = t239 * t335 - t368 * t457;
t184 = -t254 * t454 + t255 * t335;
t183 = t255 * t340 + t335 * t466;
t147 = pkin(6) * t466 + t187;
t143 = -pkin(6) * t255 - t186 * t341;
t122 = t336 * t352 + t341 * t394;
t114 = -t166 * t340 + t335 * t538;
t113 = t166 * t335 + t340 * t538;
t107 = -mrSges(7,2) * t203 + mrSges(7,3) * t160;
t83 = -t143 * t340 + t147 * t335;
t82 = t143 * t335 + t147 * t340;
t69 = t392 * t443 + (-pkin(6) * t248 - t175 * t341 + t245 * t442) * t337;
t53 = pkin(6) * t121 + t175 * t342 - t245 * t445;
t51 = -t335 * t393 + t340 * t351;
t50 = t335 * t351 + t340 * t393;
t44 = -mrSges(6,2) * t151 + mrSges(6,3) * t81;
t15 = -mrSges(7,2) * t77 + mrSges(7,3) * t30;
t14 = mrSges(7,1) * t77 - mrSges(7,3) * t29;
t13 = -qJD(6) * t112 + t335 * t69 + t340 * t53;
t12 = qJD(6) * t111 + t335 * t53 - t340 * t69;
t1 = [((t353 * qJD(4) + m(5) * (-qJD(4) * t187 + t73) + m(6) * (-t167 * t444 + t73) + t475) * t342 + (-t341 * t44 + t527 * t336 + (t336 * t162 + t341 * t476) * qJD(5) + m(5) * (-t72 - t446) + m(6) * (-t446 + t601) - t594) * t337) * t245 + (t509 / 0.2e1 + t505 / 0.2e1 + (-t510 / 0.2e1 - t506 / 0.2e1) * qJD(2) + (t404 * t344 + (-t248 * t343 + t249 * t338) * mrSges(4,3)) * pkin(3) - t595) * qJD(2) - (t279 * mrSges(4,2) - t294 * mrSges(4,3) + Ifges(4,1) * t213 + Ifges(4,4) * t214 + Ifges(4,5) * t332 + (t492 / 0.2e1 + t501 / 0.2e1 + t499 / 0.2e1 + t276 * mrSges(5,2) + t60 / 0.2e1 + t73 * mrSges(5,3)) * t342 + (t276 * mrSges(5,1) + 0.2e1 * t581 + t585 - t532 + t148 / 0.2e1 + t76 / 0.2e1 + t75 / 0.2e1 + t38 * mrSges(6,1) - t533) * t337 + ((-t487 / 0.2e1 - t252 / 0.2e1 - t488 / 0.2e1 + t563 + t597) * t337 + (t493 / 0.2e1 + t166 * mrSges(6,1) + t490 / 0.2e1 + t494 / 0.2e1 - t486 / 0.2e1 - t489 / 0.2e1 - t519 / 0.2e1 + t564 + t572 + t615) * t342) * qJD(4)) * t368 - t616 * t175 + t625 * t122 / 0.2e1 + t620 * t249 / 0.2e1 + (-t279 * mrSges(4,1) + t295 * mrSges(4,3) + Ifges(4,4) * t213 + Ifges(4,2) * t214 + Ifges(4,6) * t332 + t593) * t297 + t238 * t626 + t254 * (Ifges(5,4) * t468 - Ifges(5,2) * t469 + Ifges(5,6) * t249) / 0.2e1 + t187 * (-mrSges(5,2) * t249 - mrSges(5,3) * t469) + (mrSges(6,1) * t469 - mrSges(7,1) * t51 + mrSges(7,2) * t50) * t166 + t305 * (mrSges(5,1) * t469 + mrSges(5,2) * t468) + t172 * t468 / 0.2e1 + (-mrSges(5,1) * t249 - mrSges(6,1) * t122 + mrSges(6,2) * t121 + mrSges(5,3) * t468) * t186 + t282 * (Ifges(5,5) * t468 - Ifges(5,6) * t469 + Ifges(5,3) * t249) / 0.2e1 + (-mrSges(3,1) * t302 + mrSges(3,2) * t301 + (pkin(2) * m(3) - t308) * qJDD(1)) * pkin(2) + (-t121 * t166 + t122 * t167 + t238 * t37 - t239 * t38) * mrSges(6,3) + (m(4) * t279 - mrSges(4,1) * t214 + mrSges(4,2) * t213) * t529 + t192 * t586 + t191 * t588 + t239 * t584 + (Ifges(7,1) * t192 + Ifges(7,4) * t191 + Ifges(7,5) * t238) * t583 + (Ifges(6,4) * t239 + Ifges(6,2) * t238) * t574 + (Ifges(6,1) * t239 + Ifges(6,4) * t238) * t575 + (Ifges(7,5) * t192 + Ifges(7,6) * t191 + Ifges(7,3) * t238) * t576 + t50 * t577 + t51 * t579 + (Ifges(7,4) * t192 + Ifges(7,2) * t191 + Ifges(7,6) * t238) * t582 + t469 * t564 + (Ifges(7,1) * t50 + Ifges(7,4) * t51 + Ifges(7,5) * t122) * t566 + (Ifges(7,4) * t50 + Ifges(7,2) * t51 + Ifges(7,6) * t122) * t567 + (Ifges(6,5) * t239 + Ifges(6,6) * t238) * t569 + t469 * t572 + (Ifges(5,1) * t468 - Ifges(5,4) * t469 + Ifges(5,5) * t249) * t549 + (Ifges(6,5) * t121 + Ifges(6,6) * t122 + Ifges(6,3) * t469) * t552 + (Ifges(6,1) * t121 + Ifges(6,4) * t122 + Ifges(6,5) * t469) * t557 + (Ifges(6,4) * t121 + Ifges(6,2) * t122 + Ifges(6,6) * t469) * t559 + (Ifges(7,5) * t50 + Ifges(7,6) * t51 + Ifges(7,3) * t122) * t561 + (Ifges(4,1) * t248 + Ifges(4,4) * t249) * t546 + (-Ifges(3,4) * t301 - Ifges(3,2) * t302 - Ifges(3,6) * qJDD(2)) * t339 + qJD(3) * (t505 + t509) / 0.2e1 + t309 * (-mrSges(4,1) * t249 + mrSges(4,2) * t248) + t288 * (Ifges(4,4) * t248 + Ifges(4,2) * t249) / 0.2e1 + Ifges(2,3) * qJDD(1) + t248 * t231 / 0.2e1 + t73 * (-mrSges(6,1) * t238 + mrSges(6,2) * t239) + t3 * (-mrSges(7,2) * t238 + mrSges(7,3) * t191) + t4 * (mrSges(7,1) * t238 - mrSges(7,3) * t192) + t38 * (-mrSges(7,1) * t191 + mrSges(7,2) * t192) + t121 * t120 / 0.2e1 + t61 * (mrSges(7,1) * t122 - mrSges(7,3) * t50) - t62 * (-mrSges(7,2) * t122 + mrSges(7,3) * t51) + t111 * t14 + t112 * t15 + t12 * t107 + t13 * t108 + m(7) * (t111 * t4 + t112 * t3 - t12 * t62 + t13 * t61 + t402) + m(6) * t402 + t344 * (Ifges(3,1) * t301 + Ifges(3,4) * t302 + Ifges(3,5) * qJDD(2)) - t469 * t497; t476 * t198 + t345 - t612 * t107 + t527 * t265 + 0.2e1 * (t349 * t590 + t365 * t589) * t320 + Ifges(3,5) * t301 + Ifges(3,6) * t302 + (-qJD(2) * (-t506 - t510) / 0.2e1 + t595) * qJD(1) + (mrSges(4,1) * t328 + t573 * t603 - t308 + t596) * g(3) + t611 * t108 + t266 * t44 + t211 * t14 - t370 * t15 + t197 * t162 + Ifges(3,3) * qJDD(2) + (mrSges(6,3) * t535 + (-qJD(4) * t216 + t475) * t320) * t337 + (mrSges(5,1) * t535 + t320 * t594 - t532) * t342 + t609 * t321 + (t211 * t4 - t370 * t3 + (-t270 - t603) * g(3) + t612 * t62 + t611 * t61 + t599) * m(7) + (t167 * t197 + t266 * t37 + t599) * m(6) + (t343 * (mrSges(4,1) * t332 - mrSges(4,3) * t213) + (-mrSges(4,2) * t332 + (-qJD(2) * t619 + t214) * mrSges(4,3)) * t338 + (t449 * t338 + (t342 * t216 + t337 * t450 + t258) * t343) * qJD(3) + 0.2e1 * (t429 * t589 + (t187 * t342 * t343 - t305 * t338 + t429) * t590) * qJD(3) - t404 * t447 + (g(3) * t339 + t294 * t343 + t295 * t338) * m(4)) * pkin(3) + t616 * t229; (t573 * t396 + (mrSges(6,3) * t337 + mrSges(4,1) + t524) * t328 + t596) * g(3) + t604 * t162 + t345 - t450 * t218 - m(5) * (t187 * t219 - t305 * t433 + t472) + t527 * t299 - t608 * t107 + t610 * t108 + t300 * t44 + (-t258 * t343 + (-t449 - t484) * t338) * t526 + t240 * t14 - t369 * t15 - t219 * t216 - mrSges(5,3) * t479 + t476 * t628 + t609 * pkin(4) + ((-t270 + t541) * g(3) + t240 * t4 - t369 * t3 + t608 * t62 + t610 * t61 + t600) * m(7) + (t167 * t604 + t300 * t37 - t472 + t600) * m(6) + (-m(5) * t349 - m(6) * t365 + m(7) * t534 - t342 * t106 - t475 * t337 + (t337 * t216 - t605) * qJD(4)) * pkin(5); -t476 * t471 - t450 * t187 + (t167 * t467 - t342 * t534 - t482 + (-t440 + t466) * t166 - t601) * mrSges(6,3) + (-m(6) * (t166 * t336 - t167 * t341 + t187) - t403) * t186 + (Ifges(5,1) * t550 + Ifges(5,5) * t548 + t376 * t553 + t381 * t560 + t386 * t558 - t362 + t597) * t254 + (Ifges(6,5) * t558 - Ifges(5,2) * t551 - Ifges(5,6) * t548 + Ifges(6,6) * t560 + Ifges(6,3) * t553 - t615) * t255 + (-t291 * mrSges(7,1) - t290 * mrSges(7,2) - (mrSges(7,3) + t587) * t460 + t391 - t337 * t390) * t534 + t625 * (-t442 / 0.2e1 + t467 / 0.2e1) + t593 + t341 * t626 + (-t519 + t118) * t550 - t61 * (-mrSges(7,1) * t467 - mrSges(7,3) * t184) + t62 * (mrSges(7,2) * t467 + mrSges(7,3) * t183) + (t120 + t483) * t407 + (-m(7) * t471 - mrSges(6,1) * t255 + (-t184 - t356) * mrSges(7,2) + (t183 - t355) * mrSges(7,1)) * t166 + t73 * t390 + (t107 * t355 + t14 * t459 + t15 * t462) * pkin(6) - t388 * t482 + t61 * (-mrSges(7,1) * t442 + mrSges(7,3) * t356) - t62 * (mrSges(7,2) * t442 + mrSges(7,3) * t355) + t3 * (-mrSges(7,2) * t341 + mrSges(7,3) * t462) + (t335 * t65 + t340 * t64) * t438 / 0.2e1 + (t205 * t381 + t206 * t386 + t253 * t376) * qJD(5) / 0.2e1 - m(7) * (t61 * t82 - t62 * t83) + (t172 + t252) * t551 + t4 * (mrSges(7,1) * t341 + mrSges(7,3) * t459) + (-t372 * t440 + (t3 * t335 + t340 * t4 + (-t335 * t61 - t340 * t62) * qJD(6)) * t336) * t587 + t462 * t588 + (Ifges(7,5) * t341 - t336 * t384) * t583 + t336 * t584 + t380 * t574 + t385 * t575 + (Ifges(7,3) * t341 - t336 * t374) * t576 + t184 * t578 + t183 * t580 + (Ifges(7,6) * t341 - t336 * t379) * t582 + (Ifges(7,1) * t184 + Ifges(7,4) * t183 - Ifges(7,5) * t467) * t565 + (t383 * t438 + (-Ifges(7,5) * t336 - t341 * t384) * qJD(5)) * t566 + qJD(5) * t362 + (t378 * t438 + (-Ifges(7,6) * t336 - t341 * t379) * qJD(5)) * t567 + (Ifges(7,4) * t184 + Ifges(7,2) * t183 - Ifges(7,6) * t467) * t568 + t375 * t569 + t466 * t571 + t171 * t549 + (t373 * t438 + (-Ifges(7,3) * t336 - t341 * t374) * qJD(5)) * t561 + (Ifges(7,5) * t184 + Ifges(7,6) * t183 - Ifges(7,3) * t467) * t562 - t10 * t459 / 0.2e1 + (t578 + t539) * t417 - t416 * t539 - t83 * t107 - t82 * t108; -m(7) * (t113 * t61 - t114 * t62 - t167 * t166) + t476 * t167 + t19 + (-t508 / 0.2e1 + t571 - t202 / 0.2e1 - t186 * mrSges(6,2) + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t206 + (mrSges(6,3) + t388) * t166 + t591) * t205 + (t166 * t388 + (m(7) * t372 - t335 * t107 - t340 * t108) * pkin(6) + t591) * qJD(6) + t383 * t583 + t378 * t582 + t373 * t576 + t360 * mrSges(7,3) + (m(7) * t360 - t335 * t14 + t340 * t15) * pkin(6) + (t61 * mrSges(7,1) + t62 * mrSges(7,2) - t186 * mrSges(6,1) + t504 / 0.2e1 - t507 / 0.2e1 + t502 / 0.2e1 + t503 / 0.2e1 + t119 / 0.2e1 + t63 / 0.2e1 + t167 * mrSges(6,3) + t516 / 0.2e1) * t206 + (mrSges(6,2) * t275 + t367 * (-t328 * t341 - t329 * t458)) * g(3) - t367 * t38 - t166 * t162 - t113 * t108 - t114 * t107 - t533 + t335 * t586 + t9 * t544; -t3 * mrSges(7,2) + t4 * mrSges(7,1) - t166 * (-mrSges(7,1) * t161 + mrSges(7,2) * t160) + (Ifges(7,1) * t160 + t513) * t565 + t64 * t566 + (Ifges(7,5) * t160 + Ifges(7,6) * t161) * t562 - t61 * t107 - t62 * t108 - g(3) * ((t275 * t335 + t329 * t457) * mrSges(7,1) + (t275 * t340 - t329 * t461) * mrSges(7,2)) + (t160 * t61 + t161 * t62) * mrSges(7,3) + t8 + (Ifges(7,2) * t161 + t155 + t65) * t568;];
tau = t1;
