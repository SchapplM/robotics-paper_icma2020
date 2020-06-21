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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 08:01:49
% EndTime: 2020-06-20 08:02:06
% DurationCPUTime: 10.03s
% Computational Cost: add. (48831->298), mult. (99296->382), div. (0->0), fcn. (75069->10), ass. (0->115)
t440 = qJDD(1) * pkin(2);
t450 = cos(qJ(2));
t451 = qJD(1) ^ 2;
t463 = t450 * t451;
t444 = sin(qJ(3));
t445 = sin(qJ(2));
t449 = cos(qJ(3));
t428 = (-t445 * t444 + t450 * t449) * qJD(1);
t460 = qJD(1) * qJD(2);
t459 = t450 * t460;
t430 = -t445 * qJDD(1) - t459;
t431 = t450 * qJDD(1) - t445 * t460;
t398 = -t428 * qJD(3) + t449 * t430 - t444 * t431;
t427 = (-t450 * t444 - t445 * t449) * qJD(1);
t399 = t427 * qJD(3) + t444 * t430 + t449 * t431;
t416 = t440 + (-t430 + t459) * pkin(3);
t439 = qJD(2) + qJD(3);
t359 = (t427 * t439 + t399) * pkin(5) + (t428 * t439 - t398) * pkin(4) + t416;
t433 = -pkin(2) * t463 + t445 * g(3);
t420 = (-t445 * t463 + qJDD(2)) * pkin(3) + t433;
t432 = -t445 * t451 * pkin(2) - t450 * g(3);
t422 = (-t445 ^ 2 * t451 - qJD(2) ^ 2) * pkin(3) + t432;
t401 = t444 * t420 + t449 * t422;
t409 = -t427 * pkin(4) + t428 * pkin(5);
t437 = t439 ^ 2;
t438 = qJDD(2) + qJDD(3);
t369 = -t437 * pkin(4) - t438 * pkin(5) + t427 * t409 + t401;
t443 = sin(qJ(4));
t448 = cos(qJ(4));
t350 = -t443 * t359 + t448 * t369;
t400 = t449 * t420 - t444 * t422;
t368 = t438 * pkin(4) - t437 * pkin(5) - t428 * t409 + t400;
t442 = sin(qJ(5));
t447 = cos(qJ(5));
t341 = t447 * t350 + t442 * t368;
t462 = t445 * qJD(1);
t461 = t450 * qJD(1);
t412 = t448 * t428 - t443 * t439;
t376 = -t412 * qJD(4) - t443 * t399 - t448 * t438;
t375 = qJDD(5) - t376;
t423 = qJD(4) + t427;
t391 = -t442 * t412 + t447 * t423;
t392 = t447 * t412 + t442 * t423;
t337 = (t391 * t392 - t375) * pkin(6) + t341;
t349 = t448 * t359 + t443 * t369;
t411 = -t443 * t428 - t448 * t439;
t377 = t411 * qJD(4) + t448 * t399 - t443 * t438;
t397 = qJDD(4) + t398;
t357 = t391 * qJD(5) + t447 * t377 + t442 * t397;
t410 = qJD(5) - t411;
t338 = (t391 * t410 + t357) * pkin(6) + t349;
t441 = sin(qJ(6));
t446 = cos(qJ(6));
t335 = t441 * t337 + t446 * t338;
t378 = t441 * t392 + t446 * t410;
t345 = t378 * qJD(6) - t446 * t357 + t441 * t375;
t356 = -t392 * qJD(5) - t442 * t377 + t447 * t397;
t355 = qJDD(6) + t356;
t379 = -t446 * t392 + t441 * t410;
t358 = -t378 * mrSges(7,1) + t379 * mrSges(7,2);
t389 = qJD(6) + t391;
t360 = -t389 * mrSges(7,2) + t378 * mrSges(7,3);
t332 = m(7) * t335 + t355 * mrSges(7,1) - t345 * mrSges(7,3) - t379 * t358 + t389 * t360;
t336 = -t446 * t337 + t441 * t338;
t344 = -t379 * qJD(6) + t441 * t357 + t446 * t375;
t361 = t389 * mrSges(7,1) - t379 * mrSges(7,3);
t333 = m(7) * t336 - t355 * mrSges(7,2) + t344 * mrSges(7,3) + t378 * t358 - t389 * t361;
t326 = t441 * t332 - t446 * t333;
t370 = -t391 * mrSges(6,1) + t392 * mrSges(6,2);
t381 = t410 * mrSges(6,1) - t392 * mrSges(6,3);
t325 = m(6) * t341 - t375 * mrSges(6,2) + t356 * mrSges(6,3) + t391 * t370 - t410 * t381 + t326;
t340 = -t442 * t350 + t447 * t368;
t339 = (-t392 ^ 2 - t410 ^ 2) * pkin(6) + t340;
t380 = -t410 * mrSges(6,2) + t391 * mrSges(6,3);
t330 = m(6) * t340 + m(7) * t339 + t375 * mrSges(6,1) - t344 * mrSges(7,1) + t345 * mrSges(7,2) - t357 * mrSges(6,3) - t378 * t360 + t379 * t361 - t392 * t370 + t410 * t380;
t386 = -t411 * mrSges(5,1) + t412 * mrSges(5,2);
t403 = t423 * mrSges(5,1) - t412 * mrSges(5,3);
t319 = m(5) * t350 - t397 * mrSges(5,2) + t376 * mrSges(5,3) + t447 * t325 - t442 * t330 + t411 * t386 - t423 * t403;
t402 = -t423 * mrSges(5,2) + t411 * mrSges(5,3);
t458 = t446 * t332 + t441 * t333;
t321 = t397 * mrSges(5,1) + t356 * mrSges(6,1) - t357 * mrSges(6,2) - t377 * mrSges(5,3) + t391 * t380 - t392 * t381 - t412 * t386 + t423 * t402 + (-m(5) - m(6)) * t349 - t458;
t314 = t448 * t319 - t443 * t321;
t313 = -t443 * t319 - t448 * t321;
t457 = m(5) * t368 - t376 * mrSges(5,1) + t377 * mrSges(5,2) + t442 * t325 + t447 * t330 - t411 * t402 + t412 * t403;
t352 = Ifges(7,4) * t379 + Ifges(7,2) * t378 + Ifges(7,6) * t389;
t353 = Ifges(7,1) * t379 + Ifges(7,4) * t378 + Ifges(7,5) * t389;
t456 = mrSges(7,1) * t335 - mrSges(7,2) * t336 + Ifges(7,5) * t345 + Ifges(7,6) * t344 + Ifges(7,3) * t355 + t379 * t352 - t378 * t353;
t351 = Ifges(7,5) * t379 + Ifges(7,6) * t378 + Ifges(7,3) * t389;
t327 = -mrSges(7,1) * t339 + mrSges(7,3) * t336 + Ifges(7,4) * t345 + Ifges(7,2) * t344 + Ifges(7,6) * t355 - t379 * t351 + t389 * t353;
t328 = mrSges(7,2) * t339 - mrSges(7,3) * t335 + Ifges(7,1) * t345 + Ifges(7,4) * t344 + Ifges(7,5) * t355 + t378 * t351 - t389 * t352;
t362 = Ifges(6,5) * t392 + Ifges(6,6) * t391 + Ifges(6,3) * t410;
t363 = Ifges(6,4) * t392 + Ifges(6,2) * t391 + Ifges(6,6) * t410;
t317 = pkin(6) * t458 + mrSges(6,2) * t349 - mrSges(6,3) * t340 + Ifges(6,1) * t357 + Ifges(6,4) * t356 + Ifges(6,5) * t375 + t441 * t327 - t446 * t328 + t391 * t362 - t410 * t363;
t364 = Ifges(6,1) * t392 + Ifges(6,4) * t391 + Ifges(6,5) * t410;
t323 = -mrSges(6,1) * t349 + mrSges(6,3) * t341 + Ifges(6,4) * t357 + Ifges(6,2) * t356 + Ifges(6,6) * t375 - t392 * t362 + t410 * t364 + t456;
t383 = Ifges(5,4) * t412 + Ifges(5,2) * t411 + Ifges(5,6) * t423;
t384 = Ifges(5,1) * t412 + Ifges(5,4) * t411 + Ifges(5,5) * t423;
t455 = -mrSges(5,1) * t349 - mrSges(5,2) * t350 + Ifges(5,5) * t377 + Ifges(5,6) * t376 + Ifges(5,3) * t397 + t442 * t317 + t447 * t323 + t412 * t383 - t411 * t384;
t417 = -t439 * mrSges(4,2) + t427 * mrSges(4,3);
t418 = t439 * mrSges(4,1) - t428 * mrSges(4,3);
t454 = m(4) * t416 - t398 * mrSges(4,1) + t399 * mrSges(4,2) - t427 * t417 + t428 * t418 + t313;
t382 = Ifges(5,5) * t412 + Ifges(5,6) * t411 + Ifges(5,3) * t423;
t312 = mrSges(5,2) * t368 + mrSges(5,3) * t349 + Ifges(5,1) * t377 + Ifges(5,4) * t376 + Ifges(5,5) * t397 + t447 * t317 - t442 * t323 + t411 * t382 - t423 * t383;
t452 = pkin(6) * t326 - mrSges(6,1) * t340 + mrSges(6,2) * t341 - Ifges(6,5) * t357 - Ifges(6,6) * t356 - Ifges(6,3) * t375 - t446 * t327 - t441 * t328 - t392 * t363 + t391 * t364;
t315 = -mrSges(5,1) * t368 + mrSges(5,3) * t350 + Ifges(5,4) * t377 + Ifges(5,2) * t376 + Ifges(5,6) * t397 - t412 * t382 + t423 * t384 + t452;
t405 = Ifges(4,4) * t428 + Ifges(4,2) * t427 + Ifges(4,6) * t439;
t406 = Ifges(4,1) * t428 + Ifges(4,4) * t427 + Ifges(4,5) * t439;
t453 = pkin(4) * t457 - pkin(5) * t314 + mrSges(4,1) * t400 - mrSges(4,2) * t401 + Ifges(4,5) * t399 + Ifges(4,6) * t398 + Ifges(4,3) * t438 - t443 * t312 - t448 * t315 + t428 * t405 - t427 * t406;
t426 = Ifges(3,5) * qJD(2) + (t450 * Ifges(3,1) - t445 * Ifges(3,4)) * qJD(1);
t425 = Ifges(3,6) * qJD(2) + (t450 * Ifges(3,4) - t445 * Ifges(3,2)) * qJD(1);
t408 = -t427 * mrSges(4,1) + t428 * mrSges(4,2);
t404 = Ifges(4,5) * t428 + Ifges(4,6) * t427 + Ifges(4,3) * t439;
t311 = -pkin(4) * t313 - mrSges(4,1) * t416 + mrSges(4,3) * t401 + Ifges(4,4) * t399 + Ifges(4,2) * t398 + Ifges(4,6) * t438 - t428 * t404 + t439 * t406 + t455;
t310 = pkin(5) * t313 + mrSges(4,2) * t416 - mrSges(4,3) * t400 + Ifges(4,1) * t399 + Ifges(4,4) * t398 + Ifges(4,5) * t438 + t448 * t312 - t443 * t315 + t427 * t404 - t439 * t405;
t1 = [Ifges(2,3) * qJDD(1) + t450 * (mrSges(3,2) * t440 - mrSges(3,3) * t433 + Ifges(3,1) * t431 + Ifges(3,4) * t430 + Ifges(3,5) * qJDD(2) - qJD(2) * t425 + t449 * t310 - t444 * t311) - t445 * (-pkin(3) * t454 - mrSges(3,1) * t440 + mrSges(3,3) * t432 + Ifges(3,4) * t431 + Ifges(3,2) * t430 + Ifges(3,6) * qJDD(2) + qJD(2) * t426 + t444 * t310 + t449 * t311) + pkin(2) * (m(3) * t440 + t431 * mrSges(3,2) - t430 * mrSges(3,1) + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t461) * t461 + (-qJD(2) * mrSges(3,2) - mrSges(3,3) * t462) * t462 + t454); t453 + pkin(3) * (t444 * (m(4) * t401 - t438 * mrSges(4,2) + t398 * mrSges(4,3) + t427 * t408 - t439 * t418 + t314) + t449 * (m(4) * t400 + t438 * mrSges(4,1) - t399 * mrSges(4,3) - t428 * t408 + t439 * t417 + t457)) + Ifges(3,6) * t430 + Ifges(3,5) * t431 - mrSges(3,2) * t432 + mrSges(3,1) * t433 + Ifges(3,3) * qJDD(2) + (t450 * t425 + t445 * t426) * qJD(1); t453; t455; -t452; t456;];
tauJ = t1;
