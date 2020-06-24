% Calculate kinetic energy for
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = CloosQRC350DE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:01:29
% EndTime: 2020-06-23 21:01:35
% DurationCPUTime: 6.52s
% Computational Cost: add. (1329->232), mult. (1934->451), div. (0->0), fcn. (1912->14), ass. (0->158)
t399 = pkin(6) + rSges(7,3);
t389 = sin(qJ(5));
t395 = cos(qJ(4));
t426 = t395 * t389;
t359 = t399 * t426 - pkin(4);
t394 = cos(qJ(5));
t369 = t399 * t394 + pkin(5);
t391 = sin(qJ(3));
t396 = cos(qJ(3));
t316 = t359 * t396 + t369 * t391;
t315 = -pkin(3) + t316;
t392 = sin(qJ(2));
t317 = -t359 * t391 + t369 * t396;
t397 = cos(qJ(2));
t439 = t317 * t397;
t461 = t315 * t392 - t439;
t398 = cos(qJ(1));
t372 = t398 * t395;
t388 = qJ(2) + qJ(3);
t381 = sin(t388);
t390 = sin(qJ(4));
t393 = sin(qJ(1));
t430 = t393 * t390;
t343 = t381 * t372 - t430;
t382 = cos(t388);
t433 = t382 * t398;
t415 = -t343 * t389 + t394 * t433;
t460 = t415 ^ 2 / 0.2e1;
t450 = Icges(4,4) * t381;
t410 = Icges(4,2) * t382 + t450;
t321 = Icges(4,6) * t398 - t410 * t393;
t322 = Icges(4,6) * t393 + t410 * t398;
t449 = Icges(4,4) * t382;
t411 = Icges(4,1) * t381 + t449;
t323 = Icges(4,5) * t398 - t411 * t393;
t324 = Icges(4,5) * t393 + t411 * t398;
t347 = -Icges(4,2) * t381 + t449;
t348 = Icges(4,1) * t382 - t450;
t379 = qJD(2) * t393;
t360 = qJD(3) * t393 + t379;
t380 = qJD(2) * t398;
t361 = qJD(3) * t398 + t380;
t459 = -(t322 * t382 + t324 * t381) * t360 - (t321 * t382 + t323 * t381) * t361 + (t347 * t382 + t348 * t381) * qJD(1);
t338 = t381 * t430 - t372;
t458 = t338 ^ 2;
t425 = t398 * t390;
t428 = t393 * t395;
t341 = t381 * t425 + t428;
t457 = t341 ^ 2;
t387 = t398 ^ 2;
t455 = t392 * pkin(3);
t453 = Icges(5,1) * t395;
t427 = t394 * t395;
t334 = -t381 * t389 + t382 * t427;
t452 = Icges(6,1) * t334;
t376 = pkin(7) * qJ(5) - qJ(6);
t370 = sin(t376);
t371 = cos(t376);
t435 = t382 * t390;
t300 = -t334 * t371 - t370 * t435;
t451 = Icges(7,1) * t300;
t448 = Icges(3,5) * t393;
t386 = t397 ^ 2;
t447 = Icges(3,2) * t386;
t446 = Icges(5,2) * t390;
t299 = -t334 * t370 + t371 * t435;
t445 = Icges(7,2) * t299;
t375 = t382 ^ 2;
t444 = Icges(5,3) * t375;
t443 = Icges(5,3) * t398;
t340 = -t381 * t428 - t425;
t429 = t393 * t394;
t307 = -t340 * t389 - t382 * t429;
t442 = t307 * t415;
t333 = -t381 * t394 - t382 * t426;
t441 = t307 * t333;
t440 = t415 * t333;
t438 = t338 * t341;
t432 = t392 * t393;
t355 = -Icges(3,1) * t432 + Icges(3,5) * t398;
t437 = t355 * t392;
t373 = pkin(2) + t455;
t436 = t373 * t398;
t434 = t382 * t393;
t431 = t392 * t398;
t424 = qJD(1) * t393;
t423 = qJD(4) * t382;
t422 = t397 * qJD(2);
t421 = t398 * qJD(1);
t420 = qJD(2) * t455;
t419 = Icges(6,3) * t435;
t329 = t398 * t423 + t360;
t418 = t393 * t422;
t417 = t398 * t422;
t416 = -rSges(3,1) * t392 - pkin(2);
t414 = pkin(3) * t418;
t367 = -qJD(4) * t381 - qJD(1);
t301 = qJD(5) * t341 + t329;
t413 = pkin(4) * t381 + pkin(5) * t382;
t412 = rSges(4,1) * t381 + rSges(4,2) * t382;
t337 = qJD(5) * t435 + t367;
t409 = Icges(4,5) * t381 + Icges(4,6) * t382;
t314 = t317 * t392;
t405 = (t315 * t397 + t314) * qJD(2) + (t316 * t397 + t314) * qJD(3);
t330 = -t393 * t423 + t361;
t335 = t413 * t393;
t351 = t382 * pkin(4) - t381 * pkin(5);
t368 = pkin(3) * t417;
t404 = -qJD(1) * t335 + t361 * t351 - t373 * t424 + t368;
t302 = -qJD(5) * t338 + t330;
t336 = t413 * t398;
t403 = -t360 * t335 - t361 * t336 - t420;
t402 = -(Icges(4,5) * t382 - Icges(4,6) * t381) * qJD(1) + (Icges(4,3) * t398 - t409 * t393) * t361 + (Icges(4,3) * t393 + t409 * t398) * t360;
t401 = -t414 - t360 * t351 + (-t336 - t436) * qJD(1);
t385 = t393 ^ 2;
t384 = t392 ^ 2;
t383 = t390 ^ 2;
t374 = -pkin(7) * qJD(5) + qJD(6);
t364 = t398 * rSges(2,1) + t393 * rSges(2,2);
t363 = t393 * rSges(2,1) - t398 * rSges(2,2);
t362 = t397 * Icges(3,2) * t432;
t356 = Icges(3,1) * t431 + t448;
t354 = Icges(3,5) * t431 + Icges(3,3) * t393;
t353 = -Icges(3,5) * t432 + Icges(3,3) * t398;
t352 = t397 * t391 + t392 * t396;
t350 = t382 * rSges(4,1) - t381 * rSges(4,2);
t345 = -t391 * t389 + t396 * t427;
t344 = t396 * t389 + t391 * t427;
t328 = t333 ^ 2;
t326 = t393 * rSges(4,3) + t412 * t398;
t325 = t398 * rSges(4,3) - t412 * t393;
t313 = -rSges(3,1) * t418 + (-t393 * rSges(3,3) + t416 * t398) * qJD(1);
t312 = rSges(3,1) * t417 + (t398 * rSges(3,3) + t416 * t393) * qJD(1);
t310 = t343 * t394 + t389 * t433;
t308 = t340 * t394 - t389 * t434;
t304 = t307 ^ 2;
t303 = t344 * t397 + t392 * t345;
t298 = t374 * t333 + t337;
t295 = pkin(2) - t461;
t294 = -t310 * t371 - t341 * t370;
t293 = -t310 * t370 + t341 * t371;
t292 = -t308 * t371 + t338 * t370;
t291 = -t308 * t370 - t338 * t371;
t288 = -t414 - t360 * t350 + (-t326 - t436) * qJD(1);
t287 = t361 * t350 + t368 + (-t373 * t393 + t325) * qJD(1);
t286 = t360 * t325 - t361 * t326 - t420;
t285 = t374 * t307 + t302;
t284 = t374 * t415 + t301;
t282 = (t329 * t381 + t367 * t433) * rSges(5,3) + t401;
t281 = (-t330 * t381 + t367 * t434) * rSges(5,3) + t404;
t280 = (-t329 * t393 - t330 * t398) * t382 * rSges(5,3) + t403;
t276 = (-t301 * t333 + t337 * t415) * rSges(6,2) + t401;
t275 = (t302 * t333 - t307 * t337) * rSges(6,2) + t404;
t274 = t461 * qJD(2) + (t316 * t392 - t439) * qJD(3) + (-(-t344 * t392 + t397 * t345) * qJD(5) + t390 * t389 * (-t392 * t391 + t397 * t396) * qJD(4)) * t399;
t273 = (t301 * t307 - t302 * t415) * rSges(6,2) + t403;
t270 = -t295 * t424 - t405 * t398 + (-(t303 * t398 - t390 * t429) * qJD(5) + (t390 * t421 + (t352 * t425 + t428) * qJD(4)) * t389) * t399;
t269 = -t295 * t421 + t405 * t393 + ((t303 * t393 + t394 * t425) * qJD(5) + (-t390 * t424 + (-t352 * t430 + t372) * qJD(4)) * t389) * t399;
t1 = t360 * (t402 * t393 - t459 * t398) / 0.2e1 + t361 * (t459 * t393 + t402 * t398) / 0.2e1 + (Icges(6,1) * t308 ^ 2 + Icges(6,2) * t304 + Icges(6,3) * t458) * t302 ^ 2 / 0.2e1 + (Icges(6,1) * t334 ^ 2 + t375 * t383 * Icges(6,3) + Icges(6,2) * t328) * t337 ^ 2 / 0.2e1 + (Icges(7,1) * t292 ^ 2 + Icges(7,2) * t291 ^ 2 + Icges(7,3) * t304) * t285 ^ 2 / 0.2e1 + (Icges(7,1) * t300 ^ 2 + Icges(7,2) * t299 ^ 2 + Icges(7,3) * t328) * t298 ^ 2 / 0.2e1 + ((t387 * t437 + (t387 * t447 + t393 * t354 + (t356 * t392 - t398 * t447 + t353) * t398) * t393) * qJD(2) - (t448 + (Icges(3,1) - Icges(3,2)) * t431) * qJD(1) * t397) * t379 / 0.2e1 + (-(t355 * t397 + t362) * qJD(1) + ((t398 * t353 + t385 * t447) * t398 + (-t356 * t432 + (-t393 * t447 + t354 - t437) * t398) * t393) * qJD(2)) * t380 / 0.2e1 + (Icges(5,1) * t343 ^ 2 + Icges(5,2) * t457 + t387 * t444) * t329 ^ 2 / 0.2e1 + (Icges(5,1) * t340 ^ 2 + Icges(5,2) * t458 + t385 * t444) * t330 ^ 2 / 0.2e1 + m(4) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(5) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + m(7) * (t269 ^ 2 + t270 ^ 2 + t274 ^ 2) / 0.2e1 + m(3) * (t384 * qJD(2) ^ 2 * rSges(3,1) ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + m(6) * (t273 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + t329 * t330 * (t343 * Icges(5,1) * t340 - t375 * t393 * t443 - Icges(5,2) * t438) + t302 * t337 * (Icges(6,2) * t441 + t308 * t452 - t338 * t419) + t285 * t298 * (Icges(7,3) * t441 + t291 * t445 + t292 * t451) + (Icges(2,3) + m(2) * (t363 ^ 2 + t364 ^ 2)) * qJD(1) ^ 2 / 0.2e1 - ((-t381 * t322 + t382 * t324) * t360 + (-t381 * t321 + t382 * t323) * t361 + (t362 * t398 + ((-Icges(3,2) * t431 + t356) * t393 + t355 * t398) * t397) * qJD(2) + (-Icges(3,1) * t386 - t384 * Icges(3,2) + t381 * t347 - t382 * t348) * qJD(1)) * qJD(1) / 0.2e1 + ((t310 * Icges(6,1) * t308 + Icges(6,2) * t442 - Icges(6,3) * t438) * t302 + (Icges(6,2) * t440 + t310 * t452 + t341 * t419) * t337 + (Icges(6,1) * t310 ^ 2 / 0.2e1 + Icges(6,2) * t460 + Icges(6,3) * t457 / 0.2e1) * t301) * t301 + ((t294 * Icges(7,1) * t292 + t293 * Icges(7,2) * t291 + Icges(7,3) * t442) * t285 + (Icges(7,3) * t440 + t293 * t445 + t294 * t451) * t298 + (Icges(7,1) * t294 ^ 2 / 0.2e1 + Icges(7,2) * t293 ^ 2 / 0.2e1 + Icges(7,3) * t460) * t284) * t284 + ((t381 ^ 2 * Icges(5,3) / 0.2e1 + (Icges(5,1) * t395 ^ 2 + Icges(5,2) * t383) * t375 / 0.2e1) * t367 + ((Icges(5,3) * t381 * t393 - t338 * t446 + t340 * t453) * t330 + (t341 * t446 + t343 * t453 - t381 * t443) * t329) * t382) * t367;
T = t1;
