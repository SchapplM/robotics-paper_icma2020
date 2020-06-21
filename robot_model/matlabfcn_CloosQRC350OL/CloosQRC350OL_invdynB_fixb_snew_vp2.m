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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 08:02:06
% EndTime: 2020-06-20 08:02:35
% DurationCPUTime: 22.38s
% Computational Cost: add. (147520->337), mult. (299733->433), div. (0->0), fcn. (226880->12), ass. (0->132)
t596 = qJDD(1) * pkin(2);
t607 = cos(qJ(2));
t609 = qJD(1) ^ 2;
t622 = t607 * t609;
t600 = sin(qJ(3));
t601 = sin(qJ(2));
t606 = cos(qJ(3));
t578 = (-t600 * t601 + t606 * t607) * qJD(1);
t618 = qJD(1) * qJD(2);
t617 = t607 * t618;
t584 = -t601 * qJDD(1) - t617;
t585 = t607 * qJDD(1) - t601 * t618;
t548 = -t578 * qJD(3) + t606 * t584 - t600 * t585;
t587 = -pkin(2) * t622 + t601 * g(3);
t570 = (-t601 * t622 + qJDD(2)) * pkin(3) + t587;
t586 = -t601 * t609 * pkin(2) - t607 * g(3);
t572 = (-t601 ^ 2 * t609 - qJD(2) ^ 2) * pkin(3) + t586;
t551 = t600 * t570 + t606 * t572;
t577 = (-t600 * t607 - t601 * t606) * qJD(1);
t558 = -t577 * mrSges(4,1) + t578 * mrSges(4,2);
t595 = qJD(2) + qJD(3);
t568 = t595 * mrSges(4,1) - t578 * mrSges(4,3);
t594 = qJDD(2) + qJDD(3);
t549 = t577 * qJD(3) + t600 * t584 + t606 * t585;
t566 = t596 + (-t584 + t617) * pkin(3);
t517 = (t577 * t595 + t549) * pkin(5) + (t578 * t595 - t548) * pkin(4) + t566;
t559 = -t577 * pkin(4) + t578 * pkin(5);
t592 = t595 ^ 2;
t527 = -t592 * pkin(4) - t594 * pkin(5) + t577 * t559 + t551;
t599 = sin(qJ(4));
t605 = cos(qJ(4));
t509 = -t599 * t517 + t605 * t527;
t550 = t606 * t570 - t600 * t572;
t526 = t594 * pkin(4) - t592 * pkin(5) - t578 * t559 + t550;
t598 = sin(qJ(5));
t604 = cos(qJ(5));
t504 = t604 * t509 + t598 * t526;
t561 = -t599 * t578 - t605 * t595;
t532 = t561 * qJD(4) + t605 * t549 - t599 * t594;
t562 = t605 * t578 - t599 * t595;
t573 = qJD(4) + t577;
t545 = t604 * t562 + t598 * t573;
t547 = qJDD(4) + t548;
t514 = -t545 * qJD(5) - t598 * t532 + t604 * t547;
t544 = -t598 * t562 + t604 * t573;
t528 = -t544 * mrSges(6,1) + t545 * mrSges(6,2);
t531 = -t562 * qJD(4) - t599 * t549 - t605 * t594;
t530 = qJDD(5) - t531;
t560 = qJD(5) - t561;
t536 = t560 * mrSges(6,1) - t545 * mrSges(6,3);
t500 = (t544 * t545 - t530) * pkin(6) + t504;
t508 = t605 * t517 + t599 * t527;
t515 = t544 * qJD(5) + t604 * t532 + t598 * t547;
t501 = (t544 * t560 + t515) * pkin(6) + t508;
t597 = sin(qJ(6));
t603 = cos(qJ(6));
t498 = t597 * t500 + t603 * t501;
t533 = t597 * t545 + t603 * t560;
t506 = t533 * qJD(6) - t603 * t515 + t597 * t530;
t513 = qJDD(6) + t514;
t534 = -t603 * t545 + t597 * t560;
t516 = -t533 * mrSges(7,1) + t534 * mrSges(7,2);
t543 = qJD(6) + t544;
t518 = -t543 * mrSges(7,2) + t533 * mrSges(7,3);
t496 = m(7) * t498 + t513 * mrSges(7,1) - t506 * mrSges(7,3) - t534 * t516 + t543 * t518;
t499 = -t603 * t500 + t597 * t501;
t505 = -t534 * qJD(6) + t597 * t515 + t603 * t530;
t519 = t543 * mrSges(7,1) - t534 * mrSges(7,3);
t497 = m(7) * t499 - t513 * mrSges(7,2) + t505 * mrSges(7,3) + t533 * t516 - t543 * t519;
t614 = t597 * t496 - t603 * t497;
t490 = m(6) * t504 - t530 * mrSges(6,2) + t514 * mrSges(6,3) + t544 * t528 - t560 * t536 + t614;
t503 = -t598 * t509 + t604 * t526;
t502 = (-t545 ^ 2 - t560 ^ 2) * pkin(6) + t503;
t535 = -t560 * mrSges(6,2) + t544 * mrSges(6,3);
t494 = m(6) * t503 + m(7) * t502 + t530 * mrSges(6,1) - t505 * mrSges(7,1) + t506 * mrSges(7,2) - t515 * mrSges(6,3) - t533 * t518 + t534 * t519 - t545 * t528 + t560 * t535;
t541 = -t561 * mrSges(5,1) + t562 * mrSges(5,2);
t553 = t573 * mrSges(5,1) - t562 * mrSges(5,3);
t486 = m(5) * t509 - t547 * mrSges(5,2) + t531 * mrSges(5,3) + t604 * t490 - t598 * t494 + t561 * t541 - t573 * t553;
t552 = -t573 * mrSges(5,2) + t561 * mrSges(5,3);
t613 = t603 * t496 + t597 * t497;
t487 = t547 * mrSges(5,1) + t514 * mrSges(6,1) - t515 * mrSges(6,2) - t532 * mrSges(5,3) + t544 * t535 - t545 * t536 - t562 * t541 + t573 * t552 + (-m(5) - m(6)) * t508 - t613;
t615 = t605 * t486 - t599 * t487;
t479 = m(4) * t551 - t594 * mrSges(4,2) + t548 * mrSges(4,3) + t577 * t558 - t595 * t568 + t615;
t567 = -t595 * mrSges(4,2) + t577 * mrSges(4,3);
t612 = m(5) * t526 - t531 * mrSges(5,1) + t532 * mrSges(5,2) + t598 * t490 + t604 * t494 - t561 * t552 + t562 * t553;
t484 = m(4) * t550 + t594 * mrSges(4,1) - t549 * mrSges(4,3) - t578 * t558 + t595 * t567 + t612;
t621 = t600 * t479 + t606 * t484;
t620 = qJD(1) * t601;
t619 = qJD(1) * t607;
t583 = (mrSges(3,1) * t601 + mrSges(3,2) * t607) * qJD(1);
t589 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t619;
t474 = m(3) * t586 - qJDD(2) * mrSges(3,2) + t584 * mrSges(3,3) - qJD(2) * t589 + t606 * t479 - t600 * t484 - t583 * t620;
t588 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t620;
t475 = m(3) * t587 + qJDD(2) * mrSges(3,1) - t585 * mrSges(3,3) + qJD(2) * t588 - t583 * t619 + t621;
t616 = t607 * t474 - t601 * t475;
t480 = -t599 * t486 - t605 * t487;
t611 = m(4) * t566 - t548 * mrSges(4,1) + t549 * mrSges(4,2) - t577 * t567 + t578 * t568 + t480;
t610 = m(3) * t596 - t584 * mrSges(3,1) + t585 * mrSges(3,2) + t588 * t620 + t589 * t619 + t611;
t608 = cos(qJ(1));
t602 = sin(qJ(1));
t576 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t607 - Ifges(3,4) * t601) * qJD(1);
t575 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t607 - Ifges(3,2) * t601) * qJD(1);
t574 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t607 - Ifges(3,6) * t601) * qJD(1);
t556 = Ifges(4,1) * t578 + Ifges(4,4) * t577 + Ifges(4,5) * t595;
t555 = Ifges(4,4) * t578 + Ifges(4,2) * t577 + Ifges(4,6) * t595;
t554 = Ifges(4,5) * t578 + Ifges(4,6) * t577 + Ifges(4,3) * t595;
t539 = Ifges(5,1) * t562 + Ifges(5,4) * t561 + Ifges(5,5) * t573;
t538 = Ifges(5,4) * t562 + Ifges(5,2) * t561 + Ifges(5,6) * t573;
t537 = Ifges(5,5) * t562 + Ifges(5,6) * t561 + Ifges(5,3) * t573;
t522 = Ifges(6,1) * t545 + Ifges(6,4) * t544 + Ifges(6,5) * t560;
t521 = Ifges(6,4) * t545 + Ifges(6,2) * t544 + Ifges(6,6) * t560;
t520 = Ifges(6,5) * t545 + Ifges(6,6) * t544 + Ifges(6,3) * t560;
t512 = Ifges(7,1) * t534 + Ifges(7,4) * t533 + Ifges(7,5) * t543;
t511 = Ifges(7,4) * t534 + Ifges(7,2) * t533 + Ifges(7,6) * t543;
t510 = Ifges(7,5) * t534 + Ifges(7,6) * t533 + Ifges(7,3) * t543;
t492 = mrSges(7,2) * t502 - mrSges(7,3) * t498 + Ifges(7,1) * t506 + Ifges(7,4) * t505 + Ifges(7,5) * t513 + t533 * t510 - t543 * t511;
t491 = -mrSges(7,1) * t502 + mrSges(7,3) * t499 + Ifges(7,4) * t506 + Ifges(7,2) * t505 + Ifges(7,6) * t513 - t534 * t510 + t543 * t512;
t488 = -mrSges(6,1) * t508 + mrSges(7,1) * t498 - mrSges(7,2) * t499 + mrSges(6,3) * t504 + Ifges(6,4) * t515 + Ifges(7,5) * t506 + Ifges(6,2) * t514 + Ifges(6,6) * t530 + Ifges(7,6) * t505 + Ifges(7,3) * t513 + t534 * t511 - t533 * t512 - t545 * t520 + t560 * t522;
t482 = pkin(6) * t613 + mrSges(6,2) * t508 - mrSges(6,3) * t503 + Ifges(6,1) * t515 + Ifges(6,4) * t514 + Ifges(6,5) * t530 + t597 * t491 - t603 * t492 + t544 * t520 - t560 * t521;
t481 = Ifges(5,4) * t532 + Ifges(5,2) * t531 + Ifges(5,6) * t547 - t562 * t537 + t573 * t539 - mrSges(5,1) * t526 + mrSges(5,3) * t509 - Ifges(6,5) * t515 - Ifges(6,6) * t514 - Ifges(6,3) * t530 - t545 * t521 + t544 * t522 - mrSges(6,1) * t503 + mrSges(6,2) * t504 - t597 * t492 - t603 * t491 + pkin(6) * t614;
t477 = qJDD(1) * mrSges(2,1) - t609 * mrSges(2,2) + t610;
t476 = mrSges(5,2) * t526 + mrSges(5,3) * t508 + Ifges(5,1) * t532 + Ifges(5,4) * t531 + Ifges(5,5) * t547 + t604 * t482 - t598 * t488 + t561 * t537 - t573 * t538;
t472 = Ifges(4,4) * t549 + Ifges(4,2) * t548 + Ifges(4,6) * t594 - t578 * t554 + t595 * t556 - mrSges(4,1) * t566 + mrSges(4,3) * t551 + Ifges(5,5) * t532 + Ifges(5,6) * t531 + Ifges(5,3) * t547 + t562 * t538 - t561 * t539 - mrSges(5,1) * t508 - mrSges(5,2) * t509 + t598 * t482 + t604 * t488 - pkin(4) * t480;
t471 = pkin(5) * t480 + mrSges(4,2) * t566 - mrSges(4,3) * t550 + Ifges(4,1) * t549 + Ifges(4,4) * t548 + Ifges(4,5) * t594 + t605 * t476 - t599 * t481 + t577 * t554 - t595 * t555;
t470 = -t609 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t601 * t474 + t607 * t475;
t469 = t608 * t470 - t602 * t477;
t468 = t602 * t470 + t608 * t477;
t467 = mrSges(3,2) * t596 - mrSges(3,3) * t587 + Ifges(3,1) * t585 + Ifges(3,4) * t584 + Ifges(3,5) * qJDD(2) - qJD(2) * t575 + t606 * t471 - t600 * t472 - t574 * t620;
t466 = -pkin(3) * t611 - mrSges(3,1) * t596 + mrSges(3,3) * t586 + Ifges(3,4) * t585 + Ifges(3,2) * t584 + Ifges(3,6) * qJDD(2) + qJD(2) * t576 + t600 * t471 + t606 * t472 - t574 * t619;
t465 = mrSges(2,1) * g(3) - pkin(2) * t616 + pkin(3) * t621 + Ifges(3,5) * t585 + Ifges(3,6) * t584 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t587 - mrSges(3,2) * t586 + t609 * Ifges(2,5) - pkin(5) * t615 - t599 * t476 - t605 * t481 + pkin(4) * t612 + Ifges(4,6) * t548 + Ifges(4,3) * t594 + t578 * t555 - t577 * t556 + mrSges(4,1) * t550 - mrSges(4,2) * t551 + Ifges(4,5) * t549 + Ifges(2,6) * qJDD(1) + (t607 * t575 + t601 * t576) * qJD(1);
t464 = -mrSges(2,2) * g(3) + Ifges(2,5) * qJDD(1) - t609 * Ifges(2,6) + t607 * t466 + t601 * t467;
t1 = [t469; t468; (-m(1) - m(2)) * g(3) + t616; -pkin(1) * t468 - mrSges(1,2) * g(3) + t608 * t464 - t602 * t465; pkin(1) * t469 + mrSges(1,1) * g(3) + t602 * t464 + t608 * t465; pkin(2) * t610 + Ifges(2,3) * qJDD(1) - t601 * t466 + t607 * t467;];
tauB = t1;
