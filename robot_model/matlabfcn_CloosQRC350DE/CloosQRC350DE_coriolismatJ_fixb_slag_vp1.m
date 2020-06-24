% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = CloosQRC350DE_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350DE_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:02:09
% EndTime: 2020-06-23 21:02:37
% DurationCPUTime: 21.74s
% Computational Cost: add. (9894->542), mult. (15016->811), div. (0->0), fcn. (10392->14), ass. (0->395)
t609 = qJD(2) + qJD(3);
t261 = sin(qJ(5));
t267 = cos(qJ(4));
t476 = t261 * t267;
t555 = Icges(7,3) / 0.2e1;
t363 = t476 * t555;
t260 = Icges(7,2) - Icges(7,1);
t246 = pkin(7) * qJ(5) - qJ(6);
t238 = sin(t246);
t239 = cos(t246);
t262 = sin(qJ(4));
t266 = cos(qJ(5));
t472 = t262 * t266;
t429 = t239 * t472;
t236 = t238 ^ 2;
t237 = t239 ^ 2;
t448 = t236 - t237;
t311 = -0.2e1 * t238 * t429 + t267 * t448;
t601 = t311 * t260 * t261;
t603 = t601 / 0.2e1 + t363;
t608 = qJD(4) * t603;
t604 = t363 - t601 / 0.2e1;
t607 = qJD(4) * t604;
t606 = qJD(6) * t603;
t605 = t604 * qJD(6);
t265 = sin(qJ(1));
t270 = cos(qJ(1));
t595 = t265 ^ 2 + t270 ^ 2;
t271 = pkin(6) + rSges(7,3);
t359 = m(6) * rSges(6,2) ^ 2 + m(7) * t271 ^ 2;
t591 = Icges(7,3) + Icges(6,2);
t570 = Icges(6,1) - t591;
t214 = Icges(7,2) + t359 + t570;
t474 = t262 * t237;
t478 = t261 * t262;
t602 = t214 * t478 + (pkin(7) * t311 - t261 * t474) * t260;
t559 = t266 ^ 2;
t247 = 0.1e1 + t559;
t257 = t267 ^ 2;
t458 = t266 * t267;
t489 = t238 * t239;
t600 = ((t247 * t257 - t559) * t489 + t262 * t236 * t458) * t260;
t576 = Icges(5,2) + Icges(6,3);
t599 = Icges(5,1) - t576;
t264 = sin(qJ(2));
t268 = cos(qJ(3));
t263 = sin(qJ(3));
t269 = cos(qJ(2));
t464 = t263 * t269;
t212 = t268 * t264 + t464;
t328 = t212 * t261;
t453 = t268 * t269;
t467 = t263 * t264;
t339 = -t453 + t467;
t437 = t339 * t458;
t388 = 0.2e1 * t437;
t598 = t328 + t388;
t401 = -t476 / 0.2e1;
t445 = t257 - 0.1e1 / 0.2e1;
t579 = t266 * t445;
t111 = t212 * t401 - t339 * t579;
t477 = t261 * t266;
t487 = t247 * t267;
t136 = -t212 * t477 - t339 * t487;
t431 = t262 * t489;
t491 = t236 * t260;
t88 = t111 * t491;
t302 = t88 + (-t111 * t237 - t136 * t431) * t260;
t259 = qJ(2) + qJ(3);
t248 = sin(t259);
t249 = cos(t259);
t331 = t248 * t476 - t249 * t266;
t320 = Icges(7,3) * t331;
t308 = t320 / 0.2e1;
t31 = t308 + t302;
t457 = t267 * t237;
t417 = t260 * t457;
t87 = t417 * t472 - t600;
t537 = t609 * t87;
t597 = qJD(1) * t31 - t537;
t166 = t263 * t401 + t268 * t579;
t545 = 0.2e1 * t257 - 0.1e1;
t410 = t545 * t266;
t455 = t268 * t261;
t168 = t263 * t410 + t267 * t455;
t549 = -t264 / 0.2e1;
t108 = t166 * t269 + t168 * t549;
t454 = t268 * t266;
t465 = t263 * t267;
t183 = t247 * t465 + t261 * t454;
t456 = t267 * t268;
t466 = t263 * t266;
t184 = t247 * t456 - t261 * t466;
t121 = -t183 * t264 + t184 * t269;
t81 = t108 * t491;
t303 = t81 + (-t108 * t237 - t121 * t431) * t260;
t28 = t308 + t303;
t596 = qJD(1) * t28 - t537;
t460 = t265 * t267;
t470 = t262 * t270;
t197 = t248 * t460 + t470;
t503 = t197 * t261;
t241 = pkin(3) * t264 + pkin(2);
t392 = t248 * pkin(4) + t241;
t592 = t265 * (t249 * (t266 * rSges(6,2) + pkin(5)) + t392);
t93 = rSges(6,2) * t503 - t592;
t452 = t270 * t267;
t473 = t262 * t265;
t198 = -t248 * t452 + t473;
t482 = t249 * t270;
t148 = t198 * t261 + t266 * t482;
t326 = (-pkin(5) * t249 - t392) * t270;
t94 = -rSges(6,2) * t148 + t326;
t594 = m(6) * (t265 * t93 + t270 * t94);
t323 = t214 * t559 - t359;
t154 = -Icges(7,1) + t323 + t591 + t599;
t490 = t237 * t260;
t407 = -t490 / 0.2e1;
t588 = t339 * t247 * t407 + (t467 / 0.2e1 - t453 / 0.2e1) * t154;
t553 = t248 / 0.2e1;
t256 = t262 ^ 2;
t585 = t256 / 0.2e1;
t556 = rSges(4,3) * m(4);
t244 = rSges(4,2) * t556 - Icges(4,6);
t372 = rSges(4,1) * t556 - Icges(4,5);
t349 = (-t244 * t268 - t263 * t372) * t269;
t584 = t349 / 0.2e1;
t560 = t261 ^ 2;
t583 = -t560 / 0.2e1;
t557 = m(6) * rSges(6,2);
t234 = m(7) * t271 + t557;
t500 = t214 * t266;
t355 = pkin(5) * t234 + t500;
t540 = t263 * pkin(4);
t128 = -t234 * t540 - t268 * t355;
t61 = -t128 * t261 + t154 * t465;
t582 = t269 * t61;
t475 = t261 * t271;
t220 = t267 * t475 - pkin(4);
t390 = t266 * t271 + pkin(5);
t316 = -t220 * t263 + t268 * t390;
t580 = t264 * t316;
t578 = (t266 * t457 + t431) * t260;
t395 = t237 / 0.2e1 - t236 / 0.2e1;
t577 = t395 * t260;
t550 = t262 / 0.2e1;
t400 = t212 * t550;
t405 = t489 / 0.2e1;
t419 = t212 * t458;
t479 = t261 * t339;
t574 = t260 * (-t419 + t479) * t405 + t400 * t490;
t343 = t405 * t458;
t399 = -t474 / 0.2e1;
t366 = t260 * t399;
t573 = t260 * t343 + t366;
t561 = -t220 * t268 - t263 * t390;
t309 = pkin(3) + t561;
t104 = t264 * t309 + t269 * t316 + pkin(2);
t469 = t262 * t271;
t423 = t265 * t469;
t65 = -t104 * t270 - t261 * t423;
t420 = t270 * t469;
t66 = -t104 * t265 + t261 * t420;
t572 = -t265 * t66 - t270 * t65;
t226 = -t268 * pkin(5) - t540;
t569 = m(4) * t595 * (rSges(4,1) * t248 + rSges(4,2) * t249 + t241);
t547 = pkin(5) + rSges(5,3);
t568 = m(5) * t595 * (t249 * t547 + t392);
t567 = m(7) * t572;
t333 = t248 * t266 + t249 * t476;
t542 = pkin(4) * t249;
t566 = pkin(5) * t248 + rSges(6,2) * t333 - t542;
t564 = t260 * t448 + Icges(7,3);
t558 = m(3) * rSges(3,1);
t554 = -t136 / 0.2e1;
t552 = -t261 / 0.2e1;
t551 = -t262 / 0.2e1;
t548 = t267 / 0.2e1;
t546 = t257 - 0.1e1;
t544 = pkin(3) * t263;
t253 = pkin(3) * t268;
t543 = pkin(3) * t269;
t227 = pkin(4) * t234;
t541 = pkin(5) * t263;
t538 = t253 + 0.2e1 * pkin(4);
t78 = t87 * qJD(6);
t416 = t261 * t465;
t90 = (-(m(4) * rSges(4,2) + m(5) * rSges(5,3)) * t268 - m(4) * rSges(4,1) * t263 + (t416 - t454) * t234 + t226 * (m(5) + m(6) + m(7))) * pkin(3);
t535 = t90 * qJD(3) - t78;
t534 = Icges(7,3) * pkin(7);
t532 = pkin(7) * t262;
t484 = t249 * t261;
t332 = t248 * t458 + t484;
t485 = t248 * t262;
t186 = -t248 * t261 + t249 * t458;
t483 = t249 * t262;
t130 = -t186 * t238 + t239 * t483;
t518 = Icges(7,2) * t130;
t131 = t186 * t239 + t238 * t483;
t521 = Icges(7,1) * t131;
t522 = Icges(6,1) * t186;
t274 = t332 * t522 - (t238 * t332 - t239 * t485) * t518 + (t238 * t485 + t239 * t332) * t521 + (Icges(4,4) * t249 + 0.2e1 * (Icges(4,1) - Icges(4,2)) * t553) * t249 + t591 * t333 * t331 + (-0.2e1 * Icges(4,4) * t553 + (t257 * Icges(5,1) + t576 * t256 - Icges(5,3)) * t249) * t248;
t336 = t248 * t547 - t542;
t354 = -rSges(4,1) * t249 + rSges(4,2) * t248;
t1 = t274 + (t336 - t543) * t568 + (t354 - t543) * t569 + t580 * t567 - (t566 - t543) * t594 + (-t595 * (rSges(3,1) * t264 + pkin(2)) * t558 - t309 * t567 + (-Icges(3,2) + Icges(3,1)) * t264) * t269;
t531 = t1 * qJD(1);
t2 = t274 + t354 * t569 + t336 * t568 - (t269 * t561 - t580) * t567 - t566 * t594;
t530 = t2 * qJD(1);
t413 = t267 * t454;
t415 = t263 * t458;
t468 = t263 * t261;
t119 = ((-t413 + t468) * t264 + (-t415 - t455) * t269) * t271;
t3 = m(7) * (t66 * (-t119 * t265 + t266 * t420) + t65 * (-t119 * t270 - t266 * t423)) + (t93 * (t197 * t266 + t265 * t484) - t94 * (t198 * t266 - t261 * t482)) * t557 - (-pkin(7) * t130 + t239 * t333) * t521 + (-pkin(7) * t131 + t238 * t333) * t518 - t570 * t333 * t186;
t525 = t3 * qJD(1);
t427 = t249 * t478;
t377 = Icges(7,3) * t427;
t471 = t262 * t267;
t488 = t239 * t267;
t8 = -m(7) * ((-t265 * t65 + t270 * t66) * t267 + t572 * t262 * t212) * t475 - (-t93 * (t248 * t473 - t452) - t94 * (t248 * t470 + t460)) * t261 * t557 + (Icges(6,2) * t427 + t377) * t333 + ((-t238 * t267 + t429) * t521 - (t238 * t472 + t488) * t518 + t472 * t522 + t599 * t471 * t249) * t249;
t524 = t8 * qJD(1);
t375 = t261 * t437;
t254 = t266 + 0.1e1;
t255 = t266 - 0.1e1;
t481 = t254 * t255;
t125 = t212 * t481 - t375;
t432 = t260 * t489;
t381 = t125 * t432;
t435 = t214 * t477;
t523 = pkin(7) * t381 - t212 * t435;
t517 = Icges(7,3) * t186;
t516 = pkin(7) * qJD(5);
t39 = t260 * t131 * t130;
t513 = qJD(1) * t39;
t512 = t128 * t269;
t391 = t268 * pkin(4) - t541;
t222 = pkin(3) + t391;
t162 = t222 * t264 - t226 * t269 + pkin(2);
t507 = t162 * t234;
t506 = t183 * t269;
t505 = t184 * t264;
t502 = t339 * t267;
t501 = t214 * t261;
t221 = Icges(7,1) + Icges(6,3) + t359;
t497 = t221 * t212;
t496 = t226 * t264;
t242 = -pkin(5) + t544;
t495 = t234 * t242;
t494 = t234 * t261;
t493 = t234 * t266;
t492 = t234 * t268;
t463 = t264 * t128;
t190 = -t244 * t263 + t372 * t268;
t462 = t264 * t190;
t459 = t266 * t339;
t430 = t266 * t489;
t378 = t260 * t430;
t449 = t154 * t471 + t257 * t378;
t447 = Icges(7,3) * t532 + t227;
t446 = qJD(4) * t261;
t443 = pkin(7) * t490;
t442 = -t253 / 0.2e1;
t441 = ((-t263 * t579 + t268 * t401) * t269 + (t268 * t410 - t416) * t549) * t489;
t439 = t339 * t478;
t438 = t339 * t476;
t436 = t212 * t500;
t434 = t234 * t478;
t433 = t234 * t476;
t426 = t262 * t507;
t422 = Icges(7,3) * t472;
t418 = t238 * t488;
t414 = t154 * t456;
t295 = (-Icges(6,1) - Icges(7,2) - Icges(5,3) + t323) * t339;
t376 = t260 * t431;
t352 = t261 * t376;
t406 = t490 / 0.2e1;
t412 = (-t261 * t419 - t339 * t481) * t406 + t295 / 0.2e1 - t212 * t352 / 0.2e1;
t411 = pkin(7) * t600 + t267 * t352;
t409 = t328 / 0.2e1;
t408 = t494 / 0.2e1;
t404 = t247 * t585;
t403 = -t481 / 0.2e1;
t402 = t479 / 0.2e1;
t397 = -t458 / 0.2e1;
t389 = (-t254 - t255) * t212;
t387 = -0.2e1 * t413;
t386 = t261 * t445;
t383 = t262 * t443;
t382 = pkin(3) * t408;
t380 = t481 * t489;
t69 = (-t212 * t579 - t339 * t401) * t432;
t373 = t69 + t584;
t371 = t546 * t490;
t370 = qJD(5) * t432;
t369 = t266 * t409;
t368 = -t439 / 0.2e1;
t367 = t234 * t397;
t365 = t212 * t397;
t364 = t267 * t407;
t362 = t583 + t559 / 0.2e1;
t361 = t234 * t538 / 0.2e1;
t356 = pkin(7) * t376;
t353 = t260 * t380;
t351 = t401 * t534;
t243 = t253 + pkin(4);
t304 = (-t214 * t257 + t214 + t371) * t261;
t32 = -t242 * t494 + ((t234 * t243 - t383) * t267 + t304) * t266 + t411;
t107 = (t263 * t397 - t268 * t386) * t269 + (-t468 * t545 + t413) * t549;
t218 = pkin(4) * t492;
t334 = t218 + t234 * (pkin(3) - t541);
t273 = -pkin(7) * t81 + ((t269 * t334 + t463) * t266 + (t598 * t261 - t559 * t464) * t214) * t550 + ((-t263 * t559 + (t387 + t468) * t261) * t269 - t264 * (t559 * t268 + (-0.2e1 * t415 - t455) * t261)) * t366;
t306 = -(-t212 * t476 - t459) * t534 / 0.2e1;
t329 = (t222 * t269 + t496) * t234;
t279 = (-t266 * t329 + t497) * t550 + t306;
t4 = -t107 * t432 - t108 * t443 - t121 * t356 - t273 + t279 + t574;
t350 = -t4 * qJD(1) - t32 * qJD(2);
t348 = t405 * t328;
t347 = t260 * t262 * t402;
t346 = t260 * t368;
t345 = -t377 / 0.2e1;
t290 = (t329 - t436) * t476;
t292 = -t295 / 0.2e1;
t57 = t414 + t261 * (-t214 * t466 + t334);
t313 = ((t387 + t468 / 0.2e1) * t269 + (-0.4e1 * t415 - t455) * t549) * t376 + (-t264 * t61 + t269 * t57) * t548;
t10 = t290 / 0.2e1 + t292 + (-t339 * t403 + (-t121 / 0.2e1 + t369) * t267) * t490 + (t260 * t348 + t588 * t262) * t262 + t313;
t291 = (-t247 * t262 * t457 - t256 * t430) * t260 + t449;
t51 = t243 * t434 + t291;
t342 = t10 * qJD(1) + t51 * qJD(2);
t124 = t212 * t154 * t267;
t135 = -t212 * t487 + t261 * t459;
t12 = (t261 * t512 - t124) * t550 + (-t264 * t414 - t582) * t551 - t349 / 0.2e1 + (-t441 + (-t135 / 0.2e1 - t506 / 0.2e1 - t505 / 0.2e1) * t474) * t260 + t373;
t341 = qJD(1) * t12 - qJD(2) * t90;
t129 = -t263 * t355 + t218;
t340 = t129 * t269 + t463;
t338 = t395 * t472;
t330 = pkin(5) * t494 + t411;
t322 = t402 + t365;
t23 = ((t361 - t383) * t267 + t304) * t266 + t367 * t253 + t330;
t34 = ((-t383 + t227) * t267 + t304) * t266 + t330;
t109 = t339 * t386 + t365;
t276 = -pkin(7) * t88 + (t266 * t340 + t598 * t501) * t550 + (0.2e1 * t375 + (-t559 + t560) * t212) * t366;
t307 = (t269 * t391 + t496) * t234;
t278 = (-t266 * t307 + t497) * t550 + t306;
t6 = -t109 * t432 - t111 * t443 - t136 * t356 - t276 + t278 + t574;
t321 = -t6 * qJD(1) - t23 * qJD(2) - t34 * qJD(3);
t318 = (-t505 - t506) * t474;
t282 = (t388 + t409) * t376 + (-t154 * t502 + t261 * t340) * t548 + t339 * t154 * t585;
t286 = (t307 - t436) * t476;
t13 = t286 / 0.2e1 + t292 + (t262 * t348 + ((t554 + t369) * t267 - (t404 + t403) * t339) * t237) * t260 + t282;
t46 = -t256 * t378 + (-t247 * t417 + (t253 / 0.2e1 + pkin(4) + t442) * t494) * t262 + t449;
t52 = pkin(4) * t434 + t291;
t317 = t13 * qJD(1) + t46 * qJD(2) + t52 * qJD(3);
t315 = -t124 + t261 * (-t129 * t264 + t512);
t312 = (-t237 * t472 + t418) * t260;
t305 = (t399 + t343) * t260;
t15 = t261 * t389 * t407 - ((-t378 + (t555 - t577) * t261 * pkin(7)) * t262 + (-t221 / 0.2e1 + (-0.1e1 / 0.2e1 - t362) * t490 + t362 * t214) * t267) * t339 + t523;
t188 = t260 * t418;
t277 = (-t188 + (t495 + (-t214 + t490) * t266) * t262) * t266 / 0.2e1 + t602 * t261 / 0.2e1;
t284 = (-t242 * t493 + t221) * t551 + t351;
t20 = -t277 + t284 + t573;
t275 = -(t262 * t355 + t312) * t266 / 0.2e1 - t602 * t552;
t285 = (pkin(5) * t493 + t221) * t551 + t351;
t25 = -t275 + t285 + t573;
t62 = -t435 + (pkin(7) * t380 + (t254 / 0.2e1 + t255 / 0.2e1) * t261 * t237) * t260;
t301 = -t15 * qJD(1) - t20 * qJD(2) - t25 * qJD(3) - t62 * qJD(4);
t297 = (-t338 + t418) * t260;
t142 = -t437 - t328;
t296 = (-t142 * t395 + t339 * t431) * t260;
t48 = t345 + (-t125 * t489 - t395 * t439) * t260;
t289 = qJD(1) * t48 - qJD(4) * t353 + t609 * t603;
t101 = -t422 / 0.2e1 + t297;
t103 = t188 + (-Icges(7,3) / 0.2e1 - t577) * t472;
t143 = (t555 + t577) * t261;
t42 = -t517 / 0.2e1 + t296;
t287 = qJD(1) * t42 + qJD(2) * t103 + qJD(3) * t101 + qJD(4) * t143 + t370;
t283 = (-t261 * t431 + (-t448 * t472 - 0.2e1 * t418) * pkin(7)) * qJD(5) * t260;
t231 = t422 / 0.2e1;
t173 = -t320 / 0.2e1;
t144 = Icges(7,3) * t552 + t261 * t577;
t102 = -t260 * t338 + t188 + t231;
t100 = t231 + t297;
t49 = t236 * t346 + t237 * t347 + t345 + t381;
t45 = t262 * t268 * t382 + t361 * t478 + t291;
t41 = t517 / 0.2e1 + t296;
t30 = t173 + t302;
t27 = t173 + t303;
t26 = t305 + t275 + t285;
t24 = -t371 * t477 + t538 * t367 + (-0.2e1 * pkin(5) + t544) * t408 + t263 * t382 + ((t234 * t442 + t383) * t267 + t546 * t501) * t266 - t411;
t21 = t305 + t277 + t284;
t16 = (-t559 * t502 + (t389 + t438) * t261) * t406 + t214 * t502 * t583 + (-t214 * t459 + t507) * t397 - t339 * t364 + (-t162 * t493 + t221 * t339) * t548 + t368 * t534 - t523 + (t236 * t347 + t237 * t346) * pkin(7);
t14 = -t286 / 0.2e1 + (t267 * t554 - t339 * t404) * t490 + t282 + t412;
t11 = t121 * t364 - t290 / 0.2e1 + t588 * t256 + t313 + t412;
t9 = t584 - t462 + t373 + (-t582 - t264 * (t129 * t261 + t414) + t315) * t550 + (-t318 / 0.2e1 + t135 * t399 + t441) * t260;
t7 = ((t111 * pkin(7) + t400) * t237 + (t136 * t532 + t109 + t322) * t489) * t260 + t276 + t278;
t5 = ((t108 * pkin(7) + t400) * t237 + (t121 * t532 + t107 + t322) * t489) * t260 + t273 + t279;
t17 = [m(6) * ((-t93 * t148 - t503 * t94) * rSges(6,2) + t93 * t326 + t94 * t592) * qJD(1) - t1 * qJD(2) - t2 * qJD(3) - t8 * qJD(4) + t3 * qJD(5) + t39 * qJD(6), -t531 + ((-t57 * t264 - t582) * t262 - (pkin(3) * t556 + rSges(3,3) * t558 - Icges(3,5) + t190) * t264 + t349 + (-t318 + 0.2e1 * (-t166 * t264 - t269 * t168 / 0.2e1) * t489) * t260) * qJD(2) + t9 * qJD(3) + t11 * qJD(4) + t5 * qJD(5) + t27 * qJD(6), -t530 + t9 * qJD(2) + (t349 - t462 + 0.2e1 * t69 + (-t135 * t490 + t315) * t262) * qJD(3) + t14 * qJD(4) + t7 * qJD(5) + t30 * qJD(6), -t524 + t11 * qJD(2) + t14 * qJD(3) + (t426 - (t214 * t472 + t312) * t339) * t446 + t16 * qJD(5) + t49 * qJD(6), t525 + t5 * qJD(2) + t7 * qJD(3) + t16 * qJD(4) + ((-t212 * t266 + t438) * t432 + t261 * t426) * qJD(5) + t41 * qJD(6) + (-t142 * t564 - 0.2e1 * t339 * t376) * t516, qJD(2) * t27 + qJD(3) * t30 + qJD(4) * t49 + qJD(5) * t41 + t513; -qJD(3) * t12 + qJD(4) * t10 - qJD(5) * t4 + qJD(6) * t28 + t531, qJD(4) * t51 - qJD(5) * t32 + t535, qJD(4) * t45 + qJD(5) * t24 - t341 + t535, t45 * qJD(3) + ((t495 - t500) * t267 + t578) * t446 + t21 * qJD(5) + t605 + t342, t24 * qJD(3) + t21 * qJD(4) + (t242 * t433 - (pkin(3) * t492 + t447) * t266) * qJD(5) + t102 * qJD(6) + t283 + t350, qJD(5) * t102 + t596 + t607; qJD(2) * t12 + qJD(4) * t13 - qJD(5) * t6 + qJD(6) * t31 + t530, qJD(4) * t46 - qJD(5) * t23 + t341 - t78, qJD(4) * t52 - qJD(5) * t34 - t78, -(t267 * t355 - t578) * t446 + t26 * qJD(5) + t605 + t317, t26 * qJD(4) + (-pkin(5) * t433 - t266 * t447) * qJD(5) + t100 * qJD(6) + t283 + t321, qJD(5) * t100 + t597 + t607; -qJD(2) * t10 - qJD(3) * t13 - qJD(5) * t15 - qJD(6) * t48 + t524, -qJD(3) * t46 - qJD(5) * t20 - t342 - t606, -qJD(5) * t25 - t317 - t606, -qJD(5) * t62 + qJD(6) * t353, t261 * t516 * t564 + t144 * qJD(6) - t266 * t370 + t301, qJD(5) * t144 - t289; qJD(2) * t4 + qJD(3) * t6 + qJD(4) * t15 + qJD(6) * t42 - t525, qJD(3) * t23 + qJD(4) * t20 + qJD(6) * t103 - t350, qJD(4) * t25 + qJD(6) * t101 - t321, qJD(6) * t143 - t301, (qJD(6) - t516) * t432, t287; -qJD(2) * t28 - qJD(3) * t31 + qJD(4) * t48 - qJD(5) * t42 - t513, -qJD(5) * t103 - t596 + t608, -qJD(5) * t101 - t597 + t608, -qJD(5) * t143 + t289, -t287, 0;];
Cq = t17;
