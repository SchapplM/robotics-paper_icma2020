% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = CloosQRC350OL_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:57:09
% EndTime: 2020-06-23 21:57:27
% DurationCPUTime: 18.16s
% Computational Cost: add. (120234->245), mult. (252445->311), div. (0->0), fcn. (182475->12), ass. (0->129)
t711 = sin(qJ(3));
t712 = sin(qJ(2));
t717 = cos(qJ(3));
t718 = cos(qJ(2));
t685 = (-t711 * t712 + t717 * t718) * qJD(1);
t739 = qJD(1) * qJD(2);
t736 = t718 * t739;
t691 = -t712 * qJDD(1) - t736;
t737 = t712 * t739;
t692 = t718 * qJDD(1) - t737;
t657 = -t685 * qJD(3) + t717 * t691 - t711 * t692;
t684 = (-t711 * t718 - t712 * t717) * qJD(1);
t658 = t684 * qJD(3) + t711 * t691 + t717 * t692;
t707 = qJDD(1) * pkin(2);
t676 = t707 + (-t691 + t736) * pkin(3);
t704 = qJD(2) + qJD(3);
t633 = (t684 * t704 + t658) * pkin(5) + (t685 * t704 - t657) * pkin(4) + t676;
t720 = qJD(1) ^ 2;
t744 = t718 * t720;
t694 = -pkin(2) * t744 + t712 * g(3);
t738 = t712 * t744;
t730 = qJDD(2) - t738;
t680 = t730 * pkin(3) + t694;
t693 = -t712 * t720 * pkin(2) - t718 * g(3);
t706 = t712 ^ 2;
t745 = t706 * t720;
t682 = (-qJD(2) ^ 2 - t745) * pkin(3) + t693;
t660 = t711 * t680 + t717 * t682;
t666 = -t684 * pkin(4) + t685 * pkin(5);
t701 = t704 ^ 2;
t703 = qJDD(2) + qJDD(3);
t638 = -t701 * pkin(4) - t703 * pkin(5) + t684 * t666 + t660;
t710 = sin(qJ(4));
t716 = cos(qJ(4));
t629 = -t710 * t633 + t716 * t638;
t659 = t717 * t680 - t711 * t682;
t637 = t703 * pkin(4) - t701 * pkin(5) - t685 * t666 + t659;
t709 = sin(qJ(5));
t715 = cos(qJ(5));
t624 = t715 * t629 + t709 * t637;
t672 = t716 * t685 - t710 * t704;
t640 = -t672 * qJD(4) - t710 * t658 - t716 * t703;
t639 = qJDD(5) - t640;
t683 = qJD(4) + t684;
t650 = -t709 * t672 + t715 * t683;
t651 = t715 * t672 + t709 * t683;
t750 = t650 * t651;
t732 = -t639 + t750;
t622 = t732 * pkin(6) + t624;
t628 = t716 * t633 + t710 * t638;
t671 = -t710 * t685 - t716 * t704;
t641 = t671 * qJD(4) + t716 * t658 - t710 * t703;
t656 = qJDD(4) + t657;
t631 = t650 * qJD(5) + t715 * t641 + t709 * t656;
t668 = qJD(5) - t671;
t749 = t650 * t668;
t733 = t631 + t749;
t623 = t733 * pkin(6) + t628;
t708 = sin(qJ(6));
t714 = cos(qJ(6));
t616 = t708 * t622 + t714 * t623;
t642 = t708 * t651 + t714 * t668;
t626 = t642 * qJD(6) - t714 * t631 + t708 * t639;
t647 = qJD(6) + t650;
t752 = t642 * t647;
t614 = m(7) * t616 + (-t626 + t752) * mrSges(7,3);
t617 = -t714 * t622 + t708 * t623;
t643 = -t714 * t651 + t708 * t668;
t625 = -t643 * qJD(6) + t708 * t631 + t714 * t639;
t751 = t643 * t647;
t615 = m(7) * t617 + (t625 + t751) * mrSges(7,3);
t610 = t708 * t614 - t714 * t615;
t611 = mrSges(7,3) * t617 + Ifges(7,2) * t625 + (Ifges(7,1) - Ifges(7,3)) * t751;
t612 = -mrSges(7,3) * t616 + Ifges(7,1) * t626 + (-Ifges(7,2) + Ifges(7,3)) * t752;
t755 = -pkin(6) * t610 - mrSges(6,2) * t624 + Ifges(6,3) * t639 + t714 * t611 + t708 * t612 - (Ifges(6,1) - Ifges(6,2)) * t750;
t747 = t671 * t683;
t746 = t672 * t683;
t609 = m(6) * t624 + t732 * mrSges(6,2) + t610;
t734 = -t709 * t629 + t715 * t637;
t742 = -t651 ^ 2 - t668 ^ 2;
t619 = m(6) * t734 + m(7) * (t742 * pkin(6) + t734) + (-t642 ^ 2 - t643 ^ 2) * mrSges(7,3) + t742 * mrSges(6,2);
t603 = m(5) * t629 + t715 * t609 - t709 * t619 + (t640 + t746) * mrSges(5,3);
t729 = t714 * t614 + t708 * t615;
t607 = (-m(5) - m(6)) * t628 + (-t641 + t747) * mrSges(5,3) - t733 * mrSges(6,2) - t729;
t597 = t716 * t603 - t710 * t607;
t665 = -t684 * mrSges(4,1) + t685 * mrSges(4,2);
t678 = t704 * mrSges(4,1) - t685 * mrSges(4,3);
t595 = m(4) * t660 - t703 * mrSges(4,2) + t657 * mrSges(4,3) + t684 * t665 - t704 * t678 + t597;
t677 = -t704 * mrSges(4,2) + t684 * mrSges(4,3);
t727 = t709 * t609 + t715 * t619 + m(5) * t637 + (-t671 ^ 2 - t672 ^ 2) * mrSges(5,3);
t605 = m(4) * t659 + t703 * mrSges(4,1) - t658 * mrSges(4,3) - t685 * t665 + t704 * t677 + t727;
t743 = t711 * t595 + t717 * t605;
t741 = qJD(1) * t712;
t740 = qJD(1) * t718;
t590 = m(3) * t694 + (-t692 - t737) * mrSges(3,3) + t730 * mrSges(3,1) + t743;
t695 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t740;
t591 = m(3) * t693 - mrSges(3,1) * t745 + t691 * mrSges(3,3) - qJD(2) * t695 + t717 * t595 - t711 * t605;
t735 = -t712 * t590 + t718 * t591;
t731 = -t709 * t641 + t715 * t656;
t596 = -t710 * t603 - t716 * t607;
t728 = Ifges(7,3) * (-t651 * qJD(5) + qJDD(6) + t731) + (-Ifges(7,1) + Ifges(7,2)) * t642 * t643;
t601 = Ifges(6,1) * t631 + mrSges(6,2) * t628 - t714 * t612 + t708 * t611 + pkin(6) * t729 + (-Ifges(6,2) + Ifges(6,3)) * t749;
t621 = Ifges(6,2) * t731 + (-Ifges(6,2) * qJD(5) + (Ifges(6,1) - Ifges(6,3)) * t668) * t651 + t728;
t726 = Ifges(5,3) * t656 + t709 * t601 + t715 * t621 + (-Ifges(5,1) + Ifges(5,2)) * t671 * t672;
t724 = m(4) * t676 - t657 * mrSges(4,1) + t658 * mrSges(4,2) - t684 * t677 + t685 * t678 + t596;
t723 = m(3) * t707 - t691 * mrSges(3,1) + t695 * t740 + t724;
t598 = mrSges(5,3) * t628 + Ifges(5,1) * t641 + t715 * t601 - t709 * t621 + (-Ifges(5,2) + Ifges(5,3)) * t747;
t599 = mrSges(5,3) * t629 + Ifges(5,2) * t640 + (Ifges(5,1) - Ifges(5,3)) * t746 - t755;
t662 = Ifges(4,4) * t685 + Ifges(4,2) * t684 + Ifges(4,6) * t704;
t663 = Ifges(4,1) * t685 + Ifges(4,4) * t684 + Ifges(4,5) * t704;
t722 = pkin(4) * t727 - pkin(5) * t597 + mrSges(4,1) * t659 - mrSges(4,2) * t660 + Ifges(4,5) * t658 + Ifges(4,6) * t657 + Ifges(4,3) * t703 - t710 * t598 - t716 * t599 + t685 * t662 - t684 * t663;
t697 = Ifges(3,1) * t740 + Ifges(3,5) * qJD(2);
t721 = pkin(3) * t743 + mrSges(3,1) * t694 + Ifges(3,5) * t692 + Ifges(3,3) * qJDD(2) + t697 * t741 + t722;
t719 = cos(qJ(1));
t713 = sin(qJ(1));
t696 = Ifges(3,5) * t740 + Ifges(3,3) * qJD(2);
t661 = Ifges(4,5) * t685 + Ifges(4,6) * t684 + Ifges(4,3) * t704;
t593 = qJDD(1) * mrSges(2,1) + (-mrSges(3,3) * t706 - mrSges(2,2)) * t720 + t723;
t588 = -pkin(4) * t596 - mrSges(4,1) * t676 + mrSges(4,3) * t660 + Ifges(4,4) * t658 + Ifges(4,2) * t657 + Ifges(4,6) * t703 - t685 * t661 + t704 * t663 + t726;
t587 = pkin(5) * t596 + mrSges(4,2) * t676 - mrSges(4,3) * t659 + Ifges(4,1) * t658 + Ifges(4,4) * t657 + Ifges(4,5) * t703 + t716 * t598 - t710 * t599 + t684 * t661 - t704 * t662;
t586 = -t720 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t718 * t590 + t712 * t591;
t585 = t719 * t586 - t713 * t593;
t584 = t713 * t586 + t719 * t593;
t583 = -mrSges(3,3) * t694 + Ifges(3,1) * t692 + Ifges(3,5) * qJDD(2) + t717 * t587 - t711 * t588 + (Ifges(3,2) * qJD(2) - t696) * t741;
t582 = -pkin(3) * t724 - mrSges(3,1) * t707 + mrSges(3,3) * t693 + Ifges(3,2) * t691 + qJD(2) * t697 + t711 * t587 + t717 * t588 - t696 * t740;
t581 = t721 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + (-t718 * Ifges(3,2) * t712 + Ifges(2,5)) * t720 - pkin(2) * t735;
t580 = -mrSges(2,2) * g(3) + Ifges(2,5) * qJDD(1) - t720 * Ifges(2,6) + t718 * t582 + t712 * t583;
t579 = Ifges(2,3) * qJDD(1) + t718 * t583 - t712 * t582 + pkin(2) * (-mrSges(3,3) * t745 + t723);
t1 = [t585; t584; (-m(1) - m(2)) * g(3) + t735; -pkin(1) * t584 + t719 * t580 - t713 * t581; pkin(1) * t585 + t713 * t580 + t719 * t581; t579; t579; -Ifges(3,2) * t738 + t721; t722; t726; t755; t728;];
tauJB = t1;
