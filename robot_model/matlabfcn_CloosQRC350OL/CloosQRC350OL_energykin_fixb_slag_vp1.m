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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 21:55:49
% EndTime: 2020-06-23 21:55:56
% DurationCPUTime: 6.87s
% Computational Cost: add. (1344->207), mult. (1994->403), div. (0->0), fcn. (2046->12), ass. (0->144)
t319 = qJ(2) + qJ(3);
t313 = sin(t319);
t327 = cos(qJ(4));
t329 = cos(qJ(1));
t354 = t327 * t329;
t322 = sin(qJ(4));
t324 = sin(qJ(1));
t360 = t322 * t324;
t281 = t313 * t354 + t360;
t321 = sin(qJ(5));
t326 = cos(qJ(5));
t314 = cos(t319);
t361 = t314 * t329;
t251 = -t281 * t321 + t326 * t361;
t387 = t251 ^ 2 / 0.2e1;
t276 = t313 * t360 + t354;
t386 = t276 ^ 2;
t356 = t324 * t327;
t359 = t322 * t329;
t279 = t313 * t359 - t356;
t385 = t279 ^ 2;
t383 = pkin(2) * t324;
t323 = sin(qJ(2));
t381 = pkin(3) * t323;
t380 = Icges(5,1) * t327;
t362 = t314 * t327;
t272 = -t313 * t321 + t326 * t362;
t379 = Icges(6,1) * t272;
t320 = sin(qJ(6));
t325 = cos(qJ(6));
t364 = t314 * t322;
t246 = -t272 * t325 + t320 * t364;
t378 = Icges(7,1) * t246;
t377 = Icges(4,4) * t313;
t376 = Icges(4,4) * t314;
t375 = Icges(3,5) * t324;
t374 = Icges(3,5) * t329;
t328 = cos(qJ(2));
t317 = t328 ^ 2;
t373 = Icges(3,2) * t317;
t372 = Icges(5,2) * t322;
t245 = t272 * t320 + t325 * t364;
t371 = Icges(7,2) * t245;
t310 = t314 ^ 2;
t370 = Icges(5,3) * t310;
t369 = Icges(5,3) * t313;
t278 = t313 * t356 - t359;
t363 = t314 * t324;
t249 = -t278 * t321 + t326 * t363;
t368 = t249 * t251;
t271 = -t313 * t326 - t321 * t362;
t367 = t249 * t271;
t366 = t251 * t271;
t365 = t276 * t279;
t358 = t323 * t324;
t357 = t323 * t329;
t355 = t324 * t329;
t312 = qJD(2) * t329;
t297 = qJD(3) * t329 + t312;
t353 = qJD(1) * t328;
t352 = qJD(2) * t328;
t351 = qJD(4) * t314;
t350 = Icges(6,3) * t364;
t284 = t381 * t329;
t309 = qJD(1) * t329 * pkin(2);
t348 = t324 * t352;
t349 = pkin(3) * t348 + qJD(1) * t284 + t309;
t268 = t324 * t351 + t297;
t347 = t328 * t312;
t283 = t381 * t324;
t346 = -t283 - t383;
t345 = t323 * (Icges(3,1) - Icges(3,2));
t296 = (-qJD(2) - qJD(3)) * t324;
t304 = -qJD(4) * t313 + qJD(1);
t243 = qJD(5) * t276 + t268;
t344 = pkin(4) * t313 + pkin(5) * t314;
t343 = rSges(4,1) * t313 + rSges(4,2) * t314;
t275 = qJD(5) * t364 + t304;
t342 = Icges(4,1) * t313 + t376;
t341 = Icges(4,2) * t314 + t377;
t340 = Icges(4,5) * t313 + Icges(4,6) * t314;
t267 = t329 * t351 + t296;
t242 = qJD(5) * t279 + t267;
t339 = t242 * t249 - t243 * t251;
t338 = -t242 * t271 + t251 * t275;
t337 = t243 * t271 - t249 * t275;
t274 = t344 * t329;
t289 = pkin(4) * t314 - pkin(5) * t313;
t336 = qJD(1) * t274 - t289 * t296 + t349;
t335 = (-t283 * t324 - t284 * t329) * qJD(2);
t334 = (Icges(4,5) * t314 - Icges(4,6) * t313) * qJD(1) + (Icges(4,3) * t329 + t324 * t340) * t297 + (-Icges(4,3) * t324 + t329 * t340) * t296;
t273 = t344 * t324;
t306 = pkin(3) * t347;
t333 = t297 * t289 + t306 + (-t273 + t346) * qJD(1);
t332 = t296 * t273 - t274 * t297 + t335;
t259 = Icges(4,6) * t329 + t324 * t341;
t260 = -Icges(4,6) * t324 + t329 * t341;
t261 = Icges(4,5) * t329 + t324 * t342;
t262 = -Icges(4,5) * t324 + t329 * t342;
t286 = -Icges(4,2) * t313 + t376;
t287 = Icges(4,1) * t314 - t377;
t331 = (t260 * t314 + t262 * t313) * t296 + (t259 * t314 + t261 * t313) * t297 + (t286 * t314 + t287 * t313) * qJD(1);
t318 = t329 ^ 2;
t316 = t324 ^ 2;
t315 = t322 ^ 2;
t300 = rSges(2,1) * t329 - rSges(2,2) * t324;
t299 = rSges(2,1) * t324 + rSges(2,2) * t329;
t298 = t355 * t373;
t295 = rSges(3,1) * t357 - rSges(3,3) * t324;
t294 = rSges(3,1) * t358 + rSges(3,3) * t329;
t293 = Icges(3,1) * t357 - t375;
t292 = Icges(3,1) * t358 + t374;
t291 = Icges(3,5) * t357 - Icges(3,3) * t324;
t290 = Icges(3,5) * t358 + Icges(3,3) * t329;
t288 = rSges(4,1) * t314 - rSges(4,2) * t313;
t266 = t271 ^ 2;
t264 = -rSges(4,3) * t324 + t329 * t343;
t263 = rSges(4,3) * t329 + t324 * t343;
t255 = rSges(3,1) * t348 + qJD(1) * t295 + t309;
t254 = rSges(3,1) * t347 + (-t294 - t383) * qJD(1);
t252 = t281 * t326 + t321 * t361;
t250 = t278 * t326 + t321 * t363;
t247 = t249 ^ 2;
t244 = (-t294 * t324 - t295 * t329) * qJD(2);
t241 = qJD(6) * t271 + t275;
t240 = -t252 * t325 + t279 * t320;
t239 = t252 * t320 + t279 * t325;
t238 = -t250 * t325 + t276 * t320;
t237 = t250 * t320 + t276 * t325;
t234 = qJD(6) * t249 + t243;
t233 = qJD(6) * t251 + t242;
t232 = qJD(1) * t264 - t288 * t296 + t349;
t231 = t288 * t297 + t306 + (-t263 + t346) * qJD(1);
t229 = t263 * t296 - t264 * t297 + t335;
t228 = (t267 * t313 + t304 * t361) * rSges(5,3) + t336;
t227 = (-t268 * t313 - t304 * t363) * rSges(5,3) + t333;
t224 = (t267 * t324 - t268 * t329) * t314 * rSges(5,3) + t332;
t222 = rSges(6,2) * t338 + t336;
t221 = rSges(6,2) * t337 + t333;
t218 = rSges(6,2) * t339 + t332;
t216 = (-t233 * t271 + t241 * t251) * rSges(7,3) + t338 * pkin(6) + t336;
t215 = (t234 * t271 - t241 * t249) * rSges(7,3) + t337 * pkin(6) + t333;
t214 = (t233 * t249 - t234 * t251) * rSges(7,3) + t339 * pkin(6) + t332;
t1 = (Icges(6,1) * t250 ^ 2 + Icges(6,2) * t247 + Icges(6,3) * t386) * t243 ^ 2 / 0.2e1 + (Icges(7,1) * t238 ^ 2 + Icges(7,2) * t237 ^ 2 + Icges(7,3) * t247) * t234 ^ 2 / 0.2e1 + (Icges(7,1) * t246 ^ 2 + Icges(7,2) * t245 ^ 2 + Icges(7,3) * t266) * t241 ^ 2 / 0.2e1 + (Icges(6,1) * t272 ^ 2 + Icges(6,3) * t310 * t315 + Icges(6,2) * t266) * t275 ^ 2 / 0.2e1 + (Icges(5,1) * t281 ^ 2 + Icges(5,2) * t385 + t318 * t370) * t267 ^ 2 / 0.2e1 + (Icges(5,1) * t278 ^ 2 + Icges(5,2) * t386 + t316 * t370) * t268 ^ 2 / 0.2e1 + t296 * (-t334 * t324 + t331 * t329) / 0.2e1 + t297 * (t331 * t324 + t334 * t329) / 0.2e1 - qJD(2) * t324 * ((-(-t291 * t324 + t293 * t357 + t318 * t373) * t324 + (-t290 * t324 + t292 * t357 + t298) * t329) * qJD(2) + (t329 * t345 - t375) * t353) / 0.2e1 + ((-(t291 * t329 + t293 * t358 + t298) * t324 + (t290 * t329 + t292 * t358 + t316 * t373) * t329) * qJD(2) + (t324 * t345 + t374) * t353) * t312 / 0.2e1 + m(6) * (t218 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + m(7) * (t214 ^ 2 + t215 ^ 2 + t216 ^ 2) / 0.2e1 + m(5) * (t224 ^ 2 + t227 ^ 2 + t228 ^ 2) / 0.2e1 + m(4) * (t229 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + m(3) * (t244 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t267 * t268 * (Icges(5,1) * t278 * t281 + Icges(5,2) * t365 + t355 * t370) + t243 * t275 * (Icges(6,2) * t367 + t250 * t379 + t276 * t350) + t234 * t241 * (Icges(7,3) * t367 + t237 * t371 + t238 * t378) + (Icges(2,3) + m(2) * (t299 ^ 2 + t300 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((-t260 * t313 + t262 * t314) * t296 + (-t259 * t313 + t261 * t314) * t297 + (t292 * t329 - t293 * t324) * t352 + (t323 ^ 2 * Icges(3,2) + Icges(3,1) * t317 - t313 * t286 + t314 * t287) * qJD(1)) * qJD(1) / 0.2e1 + ((Icges(6,1) * t250 * t252 + Icges(6,2) * t368 + Icges(6,3) * t365) * t243 + (Icges(6,2) * t366 + t252 * t379 + t279 * t350) * t275 + (Icges(6,1) * t252 ^ 2 / 0.2e1 + Icges(6,2) * t387 + Icges(6,3) * t385 / 0.2e1) * t242) * t242 + ((Icges(7,1) * t238 * t240 + Icges(7,2) * t237 * t239 + Icges(7,3) * t368) * t234 + (Icges(7,3) * t366 + t239 * t371 + t240 * t378) * t241 + (Icges(7,1) * t240 ^ 2 / 0.2e1 + Icges(7,2) * t239 ^ 2 / 0.2e1 + Icges(7,3) * t387) * t233) * t233 + ((Icges(5,3) * t313 ^ 2 / 0.2e1 + (Icges(5,1) * t327 ^ 2 + Icges(5,2) * t315) * t310 / 0.2e1) * t304 + ((t276 * t372 + t278 * t380 - t324 * t369) * t268 + (t279 * t372 + t281 * t380 - t329 * t369) * t267) * t314) * t304;
T = t1;
