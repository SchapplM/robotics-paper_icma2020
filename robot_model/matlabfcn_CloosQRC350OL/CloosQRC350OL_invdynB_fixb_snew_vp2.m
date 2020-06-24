% Calculate vector of inverse dynamics base forces with Newton-Euler for
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = CloosQRC350OL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:56:51
% EndTime: 2020-06-23 21:57:09
% DurationCPUTime: 17.58s
% Computational Cost: add. (91620->245), mult. (192630->310), div. (0->0), fcn. (138840->12), ass. (0->121)
t540 = qJDD(1) * pkin(2);
t544 = sin(qJ(3));
t545 = sin(qJ(2));
t550 = cos(qJ(3));
t551 = cos(qJ(2));
t523 = (-t544 * t545 + t550 * t551) * qJD(1);
t538 = qJD(2) + qJD(3);
t543 = sin(qJ(4));
t549 = cos(qJ(4));
t510 = t523 * t549 - t538 * t543;
t522 = (-t544 * t551 - t545 * t550) * qJD(1);
t521 = qJD(4) + t522;
t542 = sin(qJ(5));
t548 = cos(qJ(5));
t492 = t510 * t548 + t521 * t542;
t509 = -t523 * t543 - t538 * t549;
t506 = qJD(5) - t509;
t541 = sin(qJ(6));
t547 = cos(qJ(6));
t486 = t492 * t541 + t506 * t547;
t491 = -t510 * t542 + t548 * t521;
t489 = qJD(6) + t491;
t578 = t486 * t489;
t487 = -t492 * t547 + t506 * t541;
t577 = t487 * t489;
t576 = t491 * t492;
t575 = t491 * t506;
t574 = t509 * t521;
t573 = t510 * t521;
t539 = t545 ^ 2;
t553 = qJD(1) ^ 2;
t572 = t539 * t553;
t571 = t551 * t553;
t567 = qJD(1) * qJD(2);
t565 = t551 * t567;
t526 = -qJDD(1) * t545 - t565;
t566 = t545 * t567;
t527 = qJDD(1) * t551 - t566;
t495 = -qJD(3) * t523 + t550 * t526 - t527 * t544;
t529 = -pkin(2) * t571 + t545 * g(3);
t558 = -t545 * t571 + qJDD(2);
t518 = t558 * pkin(3) + t529;
t528 = -t545 * t553 * pkin(2) - t551 * g(3);
t520 = (-qJD(2) ^ 2 - t572) * pkin(3) + t528;
t498 = t544 * t518 + t550 * t520;
t503 = -mrSges(4,1) * t522 + mrSges(4,2) * t523;
t516 = mrSges(4,1) * t538 - mrSges(4,3) * t523;
t537 = qJDD(2) + qJDD(3);
t496 = qJD(3) * t522 + t526 * t544 + t527 * t550;
t514 = t540 + (-t526 + t565) * pkin(3);
t477 = (t522 * t538 + t496) * pkin(5) + (t523 * t538 - t495) * pkin(4) + t514;
t504 = -pkin(4) * t522 + pkin(5) * t523;
t535 = t538 ^ 2;
t482 = -pkin(4) * t535 - pkin(5) * t537 + t504 * t522 + t498;
t475 = -t477 * t543 + t482 * t549;
t497 = t550 * t518 - t520 * t544;
t481 = pkin(4) * t537 - pkin(5) * t535 - t504 * t523 + t497;
t470 = t548 * t475 + t542 * t481;
t484 = -qJD(4) * t510 - t496 * t543 - t537 * t549;
t483 = qJDD(5) - t484;
t560 = -t483 + t576;
t468 = t560 * pkin(6) + t470;
t474 = t549 * t477 + t543 * t482;
t485 = qJD(4) * t509 + t496 * t549 - t537 * t543;
t494 = qJDD(4) + t495;
t476 = qJD(5) * t491 + t485 * t548 + t494 * t542;
t561 = t476 + t575;
t469 = t561 * pkin(6) + t474;
t463 = t468 * t541 + t469 * t547;
t472 = qJD(6) * t486 - t476 * t547 + t483 * t541;
t461 = m(7) * t463 + (-t472 + t578) * mrSges(7,3);
t464 = -t468 * t547 + t469 * t541;
t471 = -qJD(6) * t487 + t476 * t541 + t483 * t547;
t462 = m(7) * t464 + (t471 + t577) * mrSges(7,3);
t563 = t541 * t461 - t547 * t462;
t457 = m(6) * t470 + t560 * mrSges(6,2) + t563;
t562 = -t542 * t475 + t548 * t481;
t569 = -t492 ^ 2 - t506 ^ 2;
t466 = m(6) * t562 + m(7) * (t569 * pkin(6) + t562) + (-t486 ^ 2 - t487 ^ 2) * mrSges(7,3) + t569 * mrSges(6,2);
t452 = m(5) * t475 + t548 * t457 - t542 * t466 + (t484 + t573) * mrSges(5,3);
t557 = t547 * t461 + t541 * t462;
t455 = (-m(5) - m(6)) * t474 + (-t485 + t574) * mrSges(5,3) - t561 * mrSges(6,2) - t557;
t564 = t549 * t452 - t455 * t543;
t446 = m(4) * t498 - mrSges(4,2) * t537 + mrSges(4,3) * t495 + t503 * t522 - t516 * t538 + t564;
t515 = -mrSges(4,2) * t538 + mrSges(4,3) * t522;
t556 = t542 * t457 + t548 * t466 + m(5) * t481 + (-t509 ^ 2 - t510 ^ 2) * mrSges(5,3);
t454 = m(4) * t497 + t537 * mrSges(4,1) - t496 * mrSges(4,3) - t523 * t503 + t538 * t515 + t556;
t570 = t544 * t446 + t550 * t454;
t568 = qJD(1) * t551;
t559 = -t485 * t542 + t548 * t494;
t447 = -t543 * t452 - t549 * t455;
t555 = m(4) * t514 - t495 * mrSges(4,1) + t496 * mrSges(4,2) - t522 * t515 + t523 * t516 + t447;
t530 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t568;
t554 = m(3) * t540 - t526 * mrSges(3,1) + t530 * t568 + t555;
t552 = cos(qJ(1));
t546 = sin(qJ(1));
t532 = Ifges(3,1) * t568 + Ifges(3,5) * qJD(2);
t531 = Ifges(3,5) * t568 + Ifges(3,3) * qJD(2);
t501 = Ifges(4,1) * t523 + Ifges(4,4) * t522 + Ifges(4,5) * t538;
t500 = Ifges(4,4) * t523 + Ifges(4,2) * t522 + Ifges(4,6) * t538;
t499 = Ifges(4,5) * t523 + Ifges(4,6) * t522 + Ifges(4,3) * t538;
t467 = Ifges(6,2) * t559 + Ifges(7,3) * (qJDD(6) + t559) + ((Ifges(6,1) - Ifges(6,3)) * t506 + (-Ifges(6,2) - Ifges(7,3)) * qJD(5)) * t492 + (-Ifges(7,1) + Ifges(7,2)) * t487 * t486;
t459 = -mrSges(7,3) * t463 + Ifges(7,1) * t472 + (-Ifges(7,2) + Ifges(7,3)) * t578;
t458 = mrSges(7,3) * t464 + Ifges(7,2) * t471 + (Ifges(7,1) - Ifges(7,3)) * t577;
t450 = Ifges(6,1) * t476 + mrSges(6,2) * t474 - t547 * t459 + t541 * t458 + pkin(6) * t557 + (-Ifges(6,2) + Ifges(6,3)) * t575;
t449 = Ifges(5,2) * t484 + mrSges(5,3) * t475 - Ifges(6,3) * t483 + mrSges(6,2) * t470 - t541 * t459 - t547 * t458 + pkin(6) * t563 + (Ifges(5,1) - Ifges(5,3)) * t573 + (Ifges(6,1) - Ifges(6,2)) * t576;
t448 = mrSges(5,3) * t474 + Ifges(5,1) * t485 + t548 * t450 - t542 * t467 + (-Ifges(5,2) + Ifges(5,3)) * t574;
t444 = qJDD(1) * mrSges(2,1) + (-mrSges(3,3) * t539 - mrSges(2,2)) * t553 + t554;
t443 = m(3) * t528 - mrSges(3,1) * t572 + mrSges(3,3) * t526 - qJD(2) * t530 + t446 * t550 - t454 * t544;
t442 = m(3) * t529 + (-t527 - t566) * mrSges(3,3) + t558 * mrSges(3,1) + t570;
t441 = t551 * t443;
t440 = -pkin(4) * t447 - mrSges(4,1) * t514 + mrSges(4,3) * t498 + Ifges(4,4) * t496 + Ifges(4,2) * t495 + Ifges(4,6) * t537 + Ifges(5,3) * t494 + t542 * t450 + t548 * t467 - t523 * t499 + t538 * t501 + (-Ifges(5,1) + Ifges(5,2)) * t510 * t509;
t439 = pkin(5) * t447 + mrSges(4,2) * t514 - mrSges(4,3) * t497 + Ifges(4,1) * t496 + Ifges(4,4) * t495 + Ifges(4,5) * t537 + t448 * t549 - t449 * t543 + t499 * t522 - t500 * t538;
t438 = -mrSges(2,1) * t553 - qJDD(1) * mrSges(2,2) + t442 * t551 + t443 * t545;
t437 = t438 * t552 - t444 * t546;
t436 = t438 * t546 + t444 * t552;
t435 = -mrSges(3,3) * t529 + Ifges(3,1) * t527 + Ifges(3,5) * qJDD(2) + t550 * t439 - t544 * t440 + (Ifges(3,2) * qJD(2) - t531) * t545 * qJD(1);
t434 = -pkin(3) * t555 - mrSges(3,1) * t540 + mrSges(3,3) * t528 + Ifges(3,2) * t526 + qJD(2) * t532 + t544 * t439 + t550 * t440 - t531 * t568;
t433 = t553 * Ifges(2,5) + mrSges(2,1) * g(3) - pkin(2) * t441 + pkin(3) * t570 + Ifges(3,5) * t527 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t529 - t543 * t448 - t549 * t449 + pkin(4) * t556 - pkin(5) * t564 + Ifges(4,5) * t496 + Ifges(4,6) * t495 + Ifges(4,3) * t537 + t523 * t500 - t522 * t501 + mrSges(4,1) * t497 - mrSges(4,2) * t498 + Ifges(2,6) * qJDD(1) + (pkin(2) * t442 - Ifges(3,2) * t571 + qJD(1) * t532) * t545;
t432 = -mrSges(2,2) * g(3) + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t553 + t434 * t551 + t435 * t545;
t1 = [t437; t436; -t545 * t442 + t441 + (-m(1) - m(2)) * g(3); -pkin(1) * t436 + t432 * t552 - t433 * t546; pkin(1) * t437 + t432 * t546 + t433 * t552; Ifges(2,3) * qJDD(1) + t551 * t435 - t545 * t434 + pkin(2) * (-mrSges(3,3) * t572 + t554);];
tauB = t1;
