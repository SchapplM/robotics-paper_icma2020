% Calculate vector of inverse dynamics joint torques for with Newton-Euler
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = CloosQRC350OL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:56:44
% EndTime: 2020-06-23 21:56:51
% DurationCPUTime: 6.02s
% Computational Cost: add. (28614->209), mult. (59815->266), div. (0->0), fcn. (43635->10), ass. (0->104)
t395 = sin(qJ(3));
t396 = sin(qJ(2));
t400 = cos(qJ(3));
t401 = cos(qJ(2));
t377 = (-t396 * t395 + t401 * t400) * qJD(1);
t417 = qJD(1) * qJD(2);
t415 = t401 * t417;
t379 = -qJDD(1) * t396 - t415;
t416 = t396 * t417;
t380 = qJDD(1) * t401 - t416;
t349 = -qJD(3) * t377 + t400 * t379 - t380 * t395;
t376 = (-t401 * t395 - t396 * t400) * qJD(1);
t350 = qJD(3) * t376 + t379 * t395 + t380 * t400;
t391 = qJDD(1) * pkin(2);
t368 = t391 + (-t379 + t415) * pkin(3);
t389 = qJD(2) + qJD(3);
t325 = (t376 * t389 + t350) * pkin(5) + (t377 * t389 - t349) * pkin(4) + t368;
t402 = qJD(1) ^ 2;
t420 = t401 * t402;
t381 = -pkin(2) * t420 + t396 * g(3);
t372 = (-t396 * t420 + qJDD(2)) * pkin(3) + t381;
t408 = -pkin(2) * t396 * t402 - g(3) * t401;
t421 = t396 ^ 2 * t402;
t374 = (-qJD(2) ^ 2 - t421) * pkin(3) + t408;
t352 = t395 * t372 + t400 * t374;
t358 = -pkin(4) * t376 + pkin(5) * t377;
t387 = t389 ^ 2;
t388 = qJDD(2) + qJDD(3);
t330 = -pkin(4) * t387 - pkin(5) * t388 + t358 * t376 + t352;
t394 = sin(qJ(4));
t399 = cos(qJ(4));
t321 = -t325 * t394 + t330 * t399;
t351 = t400 * t372 - t374 * t395;
t329 = pkin(4) * t388 - pkin(5) * t387 - t358 * t377 + t351;
t393 = sin(qJ(5));
t398 = cos(qJ(5));
t316 = t398 * t321 + t393 * t329;
t364 = t377 * t399 - t389 * t394;
t332 = -qJD(4) * t364 - t350 * t394 - t388 * t399;
t331 = qJDD(5) - t332;
t375 = qJD(4) + t376;
t342 = -t364 * t393 + t398 * t375;
t343 = t364 * t398 + t375 * t393;
t426 = t342 * t343;
t412 = -t331 + t426;
t314 = t412 * pkin(6) + t316;
t320 = t325 * t399 + t330 * t394;
t363 = -t377 * t394 - t389 * t399;
t333 = qJD(4) * t363 + t350 * t399 - t388 * t394;
t348 = qJDD(4) + t349;
t323 = qJD(5) * t342 + t333 * t398 + t348 * t393;
t360 = qJD(5) - t363;
t425 = t342 * t360;
t413 = t323 + t425;
t315 = t413 * pkin(6) + t320;
t392 = sin(qJ(6));
t397 = cos(qJ(6));
t308 = t314 * t392 + t315 * t397;
t334 = t343 * t392 + t360 * t397;
t318 = qJD(6) * t334 - t323 * t397 + t331 * t392;
t339 = qJD(6) + t342;
t428 = t334 * t339;
t306 = m(7) * t308 + (-t318 + t428) * mrSges(7,3);
t309 = -t314 * t397 + t315 * t392;
t335 = -t343 * t397 + t360 * t392;
t317 = -qJD(6) * t335 + t323 * t392 + t331 * t397;
t427 = t335 * t339;
t307 = m(7) * t309 + (t317 + t427) * mrSges(7,3);
t302 = t392 * t306 - t307 * t397;
t303 = mrSges(7,3) * t309 + Ifges(7,2) * t317 + (Ifges(7,1) - Ifges(7,3)) * t427;
t304 = -mrSges(7,3) * t308 + Ifges(7,1) * t318 + (-Ifges(7,2) + Ifges(7,3)) * t428;
t431 = -pkin(6) * t302 - mrSges(6,2) * t316 + Ifges(6,3) * t331 + t303 * t397 + t304 * t392 - (Ifges(6,1) - Ifges(6,2)) * t426;
t423 = t363 * t375;
t422 = t364 * t375;
t419 = -t343 ^ 2 - t360 ^ 2;
t418 = t401 * qJD(1);
t301 = m(6) * t316 + t412 * mrSges(6,2) + t302;
t414 = -t321 * t393 + t398 * t329;
t311 = m(6) * t414 + m(7) * (t419 * pkin(6) + t414) + (-t334 ^ 2 - t335 ^ 2) * mrSges(7,3) + t419 * mrSges(6,2);
t297 = m(5) * t321 + t301 * t398 - t311 * t393 + (t332 + t422) * mrSges(5,3);
t410 = t306 * t397 + t392 * t307;
t299 = (-m(5) - m(6)) * t320 + (-t333 + t423) * mrSges(5,3) - t413 * mrSges(6,2) - t410;
t291 = t399 * t297 - t299 * t394;
t411 = -t333 * t393 + t398 * t348;
t290 = -t394 * t297 - t399 * t299;
t409 = Ifges(7,3) * (-qJD(5) * t343 + qJDD(6) + t411) + (-Ifges(7,1) + Ifges(7,2)) * t334 * t335;
t407 = t393 * t301 + t398 * t311 + m(5) * t329 + (-t363 ^ 2 - t364 ^ 2) * mrSges(5,3);
t295 = Ifges(6,1) * t323 + mrSges(6,2) * t320 - t397 * t304 + t392 * t303 + pkin(6) * t410 + (-Ifges(6,2) + Ifges(6,3)) * t425;
t313 = Ifges(6,2) * t411 + (-Ifges(6,2) * qJD(5) + (Ifges(6,1) - Ifges(6,3)) * t360) * t343 + t409;
t406 = Ifges(5,3) * t348 + t393 * t295 + t398 * t313 + (-Ifges(5,1) + Ifges(5,2)) * t363 * t364;
t369 = -mrSges(4,2) * t389 + mrSges(4,3) * t376;
t370 = mrSges(4,1) * t389 - mrSges(4,3) * t377;
t404 = m(4) * t368 - t349 * mrSges(4,1) + t350 * mrSges(4,2) - t376 * t369 + t377 * t370 + t290;
t292 = mrSges(5,3) * t320 + Ifges(5,1) * t333 + t295 * t398 - t313 * t393 + (-Ifges(5,2) + Ifges(5,3)) * t423;
t293 = mrSges(5,3) * t321 + Ifges(5,2) * t332 + (Ifges(5,1) - Ifges(5,3)) * t422 - t431;
t354 = Ifges(4,4) * t377 + Ifges(4,2) * t376 + Ifges(4,6) * t389;
t355 = Ifges(4,1) * t377 + Ifges(4,4) * t376 + Ifges(4,5) * t389;
t403 = pkin(4) * t407 - pkin(5) * t291 + mrSges(4,1) * t351 - mrSges(4,2) * t352 + Ifges(4,5) * t350 + Ifges(4,6) * t349 + Ifges(4,3) * t388 - t394 * t292 - t399 * t293 + t377 * t354 - t376 * t355;
t383 = Ifges(3,1) * t418 + Ifges(3,5) * qJD(2);
t357 = -mrSges(4,1) * t376 + mrSges(4,2) * t377;
t353 = Ifges(4,5) * t377 + Ifges(4,6) * t376 + Ifges(4,3) * t389;
t289 = -pkin(4) * t290 - mrSges(4,1) * t368 + mrSges(4,3) * t352 + Ifges(4,4) * t350 + Ifges(4,2) * t349 + Ifges(4,6) * t388 - t353 * t377 + t355 * t389 + t406;
t288 = pkin(5) * t290 + mrSges(4,2) * t368 - mrSges(4,3) * t351 + Ifges(4,1) * t350 + Ifges(4,4) * t349 + Ifges(4,5) * t388 + t292 * t399 - t293 * t394 + t353 * t376 - t354 * t389;
t1 = [Ifges(2,3) * qJDD(1) + t401 * (-mrSges(3,3) * t381 + Ifges(3,1) * t380 + Ifges(3,5) * qJDD(2) + Ifges(3,2) * t416 + t288 * t400 - t395 * t289) - t396 * (-pkin(3) * t404 - mrSges(3,1) * t391 + mrSges(3,3) * t408 + Ifges(3,2) * t379 + qJD(2) * t383 + t395 * t288 + t400 * t289) + pkin(2) * (m(3) * t391 - t379 * mrSges(3,1) + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t418) * t418 - mrSges(3,3) * t421 + t404); t403 + pkin(3) * (t395 * (m(4) * t352 - mrSges(4,2) * t388 + mrSges(4,3) * t349 + t357 * t376 - t370 * t389 + t291) + t400 * (m(4) * t351 + mrSges(4,1) * t388 - mrSges(4,3) * t350 - t357 * t377 + t369 * t389 + t407)) + Ifges(3,3) * qJDD(2) + Ifges(3,5) * t380 + mrSges(3,1) * t381 + (-Ifges(3,2) * t420 + qJD(1) * t383) * t396; t403; t406; t431; t409;];
tauJ = t1;
