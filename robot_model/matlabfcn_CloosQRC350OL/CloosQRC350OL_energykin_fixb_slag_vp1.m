% Calculate kinetic energy for
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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = CloosQRC350OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:00:06
% EndTime: 2020-06-20 08:00:11
% DurationCPUTime: 5.44s
% Computational Cost: add. (2331->319), mult. (3295->517), div. (0->0), fcn. (3738->12), ass. (0->160)
t372 = sin(qJ(1));
t416 = pkin(2) * t372;
t371 = sin(qJ(2));
t414 = pkin(3) * t371;
t413 = Icges(3,4) * t371;
t376 = cos(qJ(2));
t412 = Icges(3,4) * t376;
t367 = qJ(2) + qJ(3);
t365 = sin(t367);
t411 = Icges(4,4) * t365;
t366 = cos(t367);
t410 = Icges(4,4) * t366;
t370 = sin(qJ(4));
t409 = t366 * t370;
t408 = t366 * t372;
t375 = cos(qJ(4));
t407 = t366 * t375;
t377 = cos(qJ(1));
t406 = t366 * t377;
t405 = t370 * t372;
t404 = t370 * t377;
t403 = t372 * t375;
t402 = t375 * t377;
t364 = qJD(2) * t377;
t347 = qJD(3) * t377 + t364;
t401 = qJD(2) * t372;
t400 = qJD(4) * t366;
t399 = pkin(3) * qJD(2) * t376;
t340 = t414 * t377;
t362 = qJD(1) * t377 * pkin(2);
t398 = qJD(1) * t340 + t372 * t399 + t362;
t324 = t372 * t400 + t347;
t339 = t414 * t372;
t397 = -t339 - t416;
t346 = (-qJD(2) - qJD(3)) * t372;
t357 = -qJD(4) * t365 + qJD(1);
t334 = t365 * t405 + t402;
t295 = qJD(5) * t334 + t324;
t396 = pkin(4) * t365 + pkin(5) * t366;
t395 = rSges(3,1) * t371 + rSges(3,2) * t376;
t394 = rSges(4,1) * t365 + rSges(4,2) * t366;
t333 = qJD(5) * t409 + t357;
t393 = Icges(3,1) * t371 + t412;
t392 = Icges(4,1) * t365 + t410;
t391 = Icges(3,2) * t376 + t413;
t390 = Icges(4,2) * t366 + t411;
t389 = Icges(3,5) * t371 + Icges(3,6) * t376;
t388 = Icges(4,5) * t365 + Icges(4,6) * t366;
t319 = Icges(3,6) * t377 + t372 * t391;
t321 = Icges(3,5) * t377 + t372 * t393;
t387 = t319 * t376 + t321 * t371;
t320 = -Icges(3,6) * t372 + t377 * t391;
t322 = -Icges(3,5) * t372 + t377 * t393;
t386 = -t320 * t376 - t322 * t371;
t349 = -Icges(3,2) * t371 + t412;
t350 = Icges(3,1) * t376 - t413;
t385 = t349 * t376 + t350 * t371;
t323 = t377 * t400 + t346;
t336 = t365 * t404 - t403;
t294 = qJD(5) * t336 + t323;
t332 = t396 * t377;
t345 = pkin(4) * t366 - pkin(5) * t365;
t384 = qJD(1) * t332 - t345 * t346 + t398;
t383 = (-t339 * t372 - t340 * t377) * qJD(2);
t382 = (Icges(4,5) * t366 - Icges(4,6) * t365) * qJD(1) + (Icges(4,3) * t377 + t372 * t388) * t347 + (-Icges(4,3) * t372 + t377 * t388) * t346;
t331 = t396 * t372;
t359 = t377 * t399;
t381 = t347 * t345 + t359 + (-t331 + t397) * qJD(1);
t380 = t346 * t331 - t332 * t347 + t383;
t310 = Icges(4,6) * t377 + t372 * t390;
t311 = -Icges(4,6) * t372 + t377 * t390;
t312 = Icges(4,5) * t377 + t372 * t392;
t313 = -Icges(4,5) * t372 + t377 * t392;
t342 = -Icges(4,2) * t365 + t410;
t343 = Icges(4,1) * t366 - t411;
t379 = (t311 * t366 + t313 * t365) * t346 + (t310 * t366 + t312 * t365) * t347 + (t342 * t366 + t343 * t365) * qJD(1);
t374 = cos(qJ(5));
t373 = cos(qJ(6));
t369 = sin(qJ(5));
t368 = sin(qJ(6));
t353 = rSges(2,1) * t377 - rSges(2,2) * t372;
t352 = rSges(3,1) * t376 - rSges(3,2) * t371;
t351 = rSges(2,1) * t372 + rSges(2,2) * t377;
t348 = Icges(3,5) * t376 - Icges(3,6) * t371;
t344 = rSges(4,1) * t366 - rSges(4,2) * t365;
t337 = t365 * t402 + t405;
t335 = t365 * t403 - t404;
t330 = -t365 * t369 + t374 * t407;
t329 = -t365 * t374 - t369 * t407;
t328 = -rSges(3,3) * t372 + t377 * t395;
t327 = rSges(3,3) * t377 + t372 * t395;
t318 = -Icges(3,3) * t372 + t377 * t389;
t317 = Icges(3,3) * t377 + t372 * t389;
t315 = -rSges(4,3) * t372 + t377 * t394;
t314 = rSges(4,3) * t377 + t372 * t394;
t306 = -rSges(5,3) * t365 + (rSges(5,1) * t375 - rSges(5,2) * t370) * t366;
t305 = -Icges(5,5) * t365 + (Icges(5,1) * t375 - Icges(5,4) * t370) * t366;
t304 = -Icges(5,6) * t365 + (Icges(5,4) * t375 - Icges(5,2) * t370) * t366;
t303 = -Icges(5,3) * t365 + (Icges(5,5) * t375 - Icges(5,6) * t370) * t366;
t301 = t337 * t374 + t369 * t406;
t300 = -t337 * t369 + t374 * t406;
t299 = t335 * t374 + t369 * t408;
t298 = -t335 * t369 + t374 * t408;
t297 = -t330 * t373 + t368 * t409;
t296 = t330 * t368 + t373 * t409;
t293 = qJD(6) * t329 + t333;
t292 = qJD(1) * t328 + t352 * t401 + t362;
t291 = t352 * t364 + (-t327 - t416) * qJD(1);
t290 = rSges(5,1) * t337 - rSges(5,2) * t336 + rSges(5,3) * t406;
t289 = rSges(5,1) * t335 - rSges(5,2) * t334 + rSges(5,3) * t408;
t288 = Icges(5,1) * t337 - Icges(5,4) * t336 + Icges(5,5) * t406;
t287 = Icges(5,1) * t335 - Icges(5,4) * t334 + Icges(5,5) * t408;
t286 = Icges(5,4) * t337 - Icges(5,2) * t336 + Icges(5,6) * t406;
t285 = Icges(5,4) * t335 - Icges(5,2) * t334 + Icges(5,6) * t408;
t284 = Icges(5,5) * t337 - Icges(5,6) * t336 + Icges(5,3) * t406;
t283 = Icges(5,5) * t335 - Icges(5,6) * t334 + Icges(5,3) * t408;
t282 = (-t327 * t372 - t328 * t377) * qJD(2);
t281 = rSges(6,1) * t330 + rSges(6,2) * t329 + rSges(6,3) * t409;
t280 = Icges(6,1) * t330 + Icges(6,4) * t329 + Icges(6,5) * t409;
t279 = Icges(6,4) * t330 + Icges(6,2) * t329 + Icges(6,6) * t409;
t278 = Icges(6,5) * t330 + Icges(6,6) * t329 + Icges(6,3) * t409;
t277 = -t301 * t373 + t336 * t368;
t276 = t301 * t368 + t336 * t373;
t275 = -t299 * t373 + t334 * t368;
t274 = t299 * t368 + t334 * t373;
t273 = qJD(6) * t298 + t295;
t272 = qJD(6) * t300 + t294;
t271 = qJD(1) * t315 - t344 * t346 + t398;
t270 = t344 * t347 + t359 + (-t314 + t397) * qJD(1);
t269 = rSges(6,1) * t301 + rSges(6,2) * t300 + rSges(6,3) * t336;
t268 = rSges(6,1) * t299 + rSges(6,2) * t298 + rSges(6,3) * t334;
t267 = Icges(6,1) * t301 + Icges(6,4) * t300 + Icges(6,5) * t336;
t266 = Icges(6,1) * t299 + Icges(6,4) * t298 + Icges(6,5) * t334;
t265 = Icges(6,4) * t301 + Icges(6,2) * t300 + Icges(6,6) * t336;
t264 = Icges(6,4) * t299 + Icges(6,2) * t298 + Icges(6,6) * t334;
t263 = Icges(6,5) * t301 + Icges(6,6) * t300 + Icges(6,3) * t336;
t262 = Icges(6,5) * t299 + Icges(6,6) * t298 + Icges(6,3) * t334;
t261 = rSges(7,1) * t297 + rSges(7,2) * t296 + rSges(7,3) * t329;
t260 = Icges(7,1) * t297 + Icges(7,4) * t296 + Icges(7,5) * t329;
t259 = Icges(7,4) * t297 + Icges(7,2) * t296 + Icges(7,6) * t329;
t258 = Icges(7,5) * t297 + Icges(7,6) * t296 + Icges(7,3) * t329;
t257 = t314 * t346 - t315 * t347 + t383;
t256 = rSges(7,1) * t277 + rSges(7,2) * t276 + rSges(7,3) * t300;
t255 = rSges(7,1) * t275 + rSges(7,2) * t274 + rSges(7,3) * t298;
t254 = Icges(7,1) * t277 + Icges(7,4) * t276 + Icges(7,5) * t300;
t253 = Icges(7,1) * t275 + Icges(7,4) * t274 + Icges(7,5) * t298;
t252 = Icges(7,4) * t277 + Icges(7,2) * t276 + Icges(7,6) * t300;
t251 = Icges(7,4) * t275 + Icges(7,2) * t274 + Icges(7,6) * t298;
t250 = Icges(7,5) * t277 + Icges(7,6) * t276 + Icges(7,3) * t300;
t249 = Icges(7,5) * t275 + Icges(7,6) * t274 + Icges(7,3) * t298;
t248 = t290 * t357 - t306 * t323 + t384;
t247 = -t289 * t357 + t306 * t324 + t381;
t246 = t289 * t323 - t290 * t324 + t380;
t245 = t269 * t333 - t281 * t294 + t384;
t244 = -t268 * t333 + t281 * t295 + t381;
t243 = t268 * t294 - t269 * t295 + t380;
t242 = t256 * t293 - t261 * t272 + (-t329 * t294 + t300 * t333) * pkin(6) + t384;
t241 = -t255 * t293 + t261 * t273 + (t329 * t295 - t298 * t333) * pkin(6) + t381;
t240 = t255 * t272 - t256 * t273 + (t298 * t294 - t300 * t295) * pkin(6) + t380;
t1 = -((-t372 * t348 + t377 * t385) * qJD(1) + (t372 ^ 2 * t318 + (t387 * t377 + (-t317 + t386) * t372) * t377) * qJD(2)) * t401 / 0.2e1 + ((t377 * t348 + t372 * t385) * qJD(1) + (t377 ^ 2 * t317 + (t386 * t372 + (-t318 + t387) * t377) * t372) * qJD(2)) * t364 / 0.2e1 + m(3) * (t282 ^ 2 + t291 ^ 2 + t292 ^ 2) / 0.2e1 + t272 * ((t300 * t250 + t276 * t252 + t277 * t254) * t272 + (t249 * t300 + t251 * t276 + t253 * t277) * t273 + (t258 * t300 + t259 * t276 + t260 * t277) * t293) / 0.2e1 + t273 * ((t250 * t298 + t252 * t274 + t254 * t275) * t272 + (t298 * t249 + t274 * t251 + t275 * t253) * t273 + (t258 * t298 + t259 * t274 + t260 * t275) * t293) / 0.2e1 + t293 * ((t250 * t329 + t252 * t296 + t254 * t297) * t272 + (t249 * t329 + t251 * t296 + t253 * t297) * t273 + (t329 * t258 + t296 * t259 + t297 * t260) * t293) / 0.2e1 + t295 * ((t263 * t334 + t265 * t298 + t267 * t299) * t294 + (t334 * t262 + t298 * t264 + t299 * t266) * t295 + (t278 * t334 + t279 * t298 + t280 * t299) * t333) / 0.2e1 + t294 * ((t336 * t263 + t300 * t265 + t301 * t267) * t294 + (t262 * t336 + t264 * t300 + t266 * t301) * t295 + (t278 * t336 + t279 * t300 + t280 * t301) * t333) / 0.2e1 + t333 * ((t263 * t409 + t265 * t329 + t267 * t330) * t294 + (t262 * t409 + t264 * t329 + t266 * t330) * t295 + (t278 * t409 + t329 * t279 + t330 * t280) * t333) / 0.2e1 + t357 * ((-t283 * t324 - t284 * t323 - t303 * t357) * t365 + ((-t286 * t370 + t288 * t375) * t323 + (-t285 * t370 + t287 * t375) * t324 + (-t304 * t370 + t305 * t375) * t357) * t366) / 0.2e1 + t323 * ((t284 * t406 - t336 * t286 + t337 * t288) * t323 + (t283 * t406 - t285 * t336 + t287 * t337) * t324 + (t303 * t406 - t304 * t336 + t305 * t337) * t357) / 0.2e1 + t324 * ((t284 * t408 - t286 * t334 + t288 * t335) * t323 + (t283 * t408 - t334 * t285 + t335 * t287) * t324 + (t303 * t408 - t304 * t334 + t305 * t335) * t357) / 0.2e1 + t346 * (-t372 * t382 + t377 * t379) / 0.2e1 + t347 * (t372 * t379 + t377 * t382) / 0.2e1 + m(7) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + m(6) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + m(5) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + m(4) * (t257 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t351 ^ 2 + t353 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((-t311 * t365 + t313 * t366) * t346 + (-t310 * t365 + t312 * t366) * t347 + (-(-t320 * t371 + t322 * t376) * t372 + (-t319 * t371 + t321 * t376) * t377) * qJD(2) + (-t365 * t342 + t366 * t343 - t349 * t371 + t376 * t350) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
