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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 08:02:35
% EndTime: 2020-06-20 08:03:07
% DurationCPUTime: 28.74s
% Computational Cost: add. (196351->338), mult. (399029->432), div. (0->0), fcn. (301949->12), ass. (0->138)
t783 = qJDD(1) * pkin(2);
t794 = cos(qJ(2));
t796 = qJD(1) ^ 2;
t812 = t794 * t796;
t787 = sin(qJ(3));
t788 = sin(qJ(2));
t793 = cos(qJ(3));
t760 = (-t787 * t788 + t793 * t794) * qJD(1);
t808 = qJD(1) * qJD(2);
t807 = t794 * t808;
t769 = -t788 * qJDD(1) - t807;
t770 = t794 * qJDD(1) - t788 * t808;
t728 = -qJD(3) * t760 + t793 * t769 - t770 * t787;
t759 = (-t787 * t794 - t788 * t793) * qJD(1);
t729 = qJD(3) * t759 + t769 * t787 + t770 * t793;
t746 = t783 + (-t769 + t807) * pkin(3);
t781 = qJD(2) + qJD(3);
t689 = (t759 * t781 + t729) * pkin(5) + (t760 * t781 - t728) * pkin(4) + t746;
t772 = -pkin(2) * t812 + t788 * g(3);
t750 = (-t788 * t812 + qJDD(2)) * pkin(3) + t772;
t771 = -t788 * t796 * pkin(2) - t794 * g(3);
t754 = (-t788 ^ 2 * t796 - qJD(2) ^ 2) * pkin(3) + t771;
t731 = t787 * t750 + t793 * t754;
t739 = -pkin(4) * t759 + pkin(5) * t760;
t778 = t781 ^ 2;
t780 = qJDD(2) + qJDD(3);
t699 = -pkin(4) * t778 - pkin(5) * t780 + t739 * t759 + t731;
t786 = sin(qJ(4));
t792 = cos(qJ(4));
t680 = -t689 * t786 + t699 * t792;
t730 = t793 * t750 - t754 * t787;
t698 = pkin(4) * t780 - pkin(5) * t778 - t739 * t760 + t730;
t785 = sin(qJ(5));
t791 = cos(qJ(5));
t671 = t791 * t680 + t785 * t698;
t742 = t760 * t792 - t781 * t786;
t706 = -qJD(4) * t742 - t729 * t786 - t780 * t792;
t705 = qJDD(5) - t706;
t755 = qJD(4) + t759;
t721 = -t742 * t785 + t791 * t755;
t722 = t742 * t791 + t755 * t785;
t667 = (t721 * t722 - t705) * pkin(6) + t671;
t679 = t792 * t689 + t786 * t699;
t741 = -t760 * t786 - t781 * t792;
t707 = qJD(4) * t741 + t729 * t792 - t780 * t786;
t727 = qJDD(4) + t728;
t687 = qJD(5) * t721 + t707 * t791 + t727 * t785;
t740 = qJD(5) - t741;
t668 = (t721 * t740 + t687) * pkin(6) + t679;
t784 = sin(qJ(6));
t790 = cos(qJ(6));
t665 = t667 * t784 + t668 * t790;
t708 = t722 * t784 + t740 * t790;
t675 = qJD(6) * t708 - t687 * t790 + t705 * t784;
t686 = -qJD(5) * t722 - t707 * t785 + t791 * t727;
t685 = qJDD(6) + t686;
t709 = -t722 * t790 + t740 * t784;
t688 = -mrSges(7,1) * t708 + mrSges(7,2) * t709;
t719 = qJD(6) + t721;
t690 = -mrSges(7,2) * t719 + mrSges(7,3) * t708;
t662 = m(7) * t665 + mrSges(7,1) * t685 - mrSges(7,3) * t675 - t688 * t709 + t690 * t719;
t666 = -t667 * t790 + t668 * t784;
t674 = -qJD(6) * t709 + t687 * t784 + t705 * t790;
t691 = mrSges(7,1) * t719 - mrSges(7,3) * t709;
t663 = m(7) * t666 - mrSges(7,2) * t685 + mrSges(7,3) * t674 + t688 * t708 - t691 * t719;
t656 = t784 * t662 - t663 * t790;
t700 = -mrSges(6,1) * t721 + mrSges(6,2) * t722;
t711 = mrSges(6,1) * t740 - mrSges(6,3) * t722;
t655 = m(6) * t671 - mrSges(6,2) * t705 + mrSges(6,3) * t686 + t700 * t721 - t711 * t740 + t656;
t670 = -t785 * t680 + t791 * t698;
t669 = (-t722 ^ 2 - t740 ^ 2) * pkin(6) + t670;
t710 = -mrSges(6,2) * t740 + mrSges(6,3) * t721;
t660 = m(6) * t670 + m(7) * t669 + mrSges(6,1) * t705 - mrSges(7,1) * t674 + mrSges(7,2) * t675 - mrSges(6,3) * t687 - t690 * t708 + t691 * t709 - t700 * t722 + t710 * t740;
t716 = -mrSges(5,1) * t741 + mrSges(5,2) * t742;
t733 = mrSges(5,1) * t755 - mrSges(5,3) * t742;
t649 = m(5) * t680 - mrSges(5,2) * t727 + mrSges(5,3) * t706 + t655 * t791 - t660 * t785 + t716 * t741 - t733 * t755;
t732 = -mrSges(5,2) * t755 + mrSges(5,3) * t741;
t805 = t790 * t662 + t784 * t663;
t651 = t727 * mrSges(5,1) + t686 * mrSges(6,1) - t687 * mrSges(6,2) - t707 * mrSges(5,3) + t721 * t710 - t722 * t711 - t742 * t716 + t755 * t732 + (-m(5) - m(6)) * t679 - t805;
t642 = t792 * t649 - t651 * t786;
t738 = -mrSges(4,1) * t759 + mrSges(4,2) * t760;
t748 = mrSges(4,1) * t781 - mrSges(4,3) * t760;
t640 = m(4) * t731 - mrSges(4,2) * t780 + mrSges(4,3) * t728 + t738 * t759 - t748 * t781 + t642;
t747 = -mrSges(4,2) * t781 + mrSges(4,3) * t759;
t804 = m(5) * t698 - mrSges(5,1) * t706 + t707 * mrSges(5,2) + t785 * t655 + t791 * t660 - t741 * t732 + t742 * t733;
t647 = m(4) * t730 + mrSges(4,1) * t780 - mrSges(4,3) * t729 - t738 * t760 + t747 * t781 + t804;
t811 = t787 * t640 + t793 * t647;
t810 = qJD(1) * t788;
t809 = qJD(1) * t794;
t767 = (mrSges(3,1) * t788 + mrSges(3,2) * t794) * qJD(1);
t774 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t809;
t634 = m(3) * t771 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t769 - qJD(2) * t774 + t640 * t793 - t647 * t787 - t767 * t810;
t773 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t810;
t635 = m(3) * t772 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t770 + qJD(2) * t773 - t767 * t809 + t811;
t806 = t794 * t634 - t788 * t635;
t641 = -t786 * t649 - t792 * t651;
t682 = Ifges(7,4) * t709 + Ifges(7,2) * t708 + Ifges(7,6) * t719;
t683 = Ifges(7,1) * t709 + Ifges(7,4) * t708 + Ifges(7,5) * t719;
t803 = mrSges(7,1) * t665 - mrSges(7,2) * t666 + Ifges(7,5) * t675 + Ifges(7,6) * t674 + Ifges(7,3) * t685 + t709 * t682 - t683 * t708;
t681 = Ifges(7,5) * t709 + Ifges(7,6) * t708 + Ifges(7,3) * t719;
t657 = -mrSges(7,1) * t669 + mrSges(7,3) * t666 + Ifges(7,4) * t675 + Ifges(7,2) * t674 + Ifges(7,6) * t685 - t681 * t709 + t683 * t719;
t658 = mrSges(7,2) * t669 - mrSges(7,3) * t665 + Ifges(7,1) * t675 + Ifges(7,4) * t674 + Ifges(7,5) * t685 + t681 * t708 - t682 * t719;
t692 = Ifges(6,5) * t722 + Ifges(6,6) * t721 + Ifges(6,3) * t740;
t693 = Ifges(6,4) * t722 + Ifges(6,2) * t721 + Ifges(6,6) * t740;
t645 = pkin(6) * t805 + mrSges(6,2) * t679 - mrSges(6,3) * t670 + Ifges(6,1) * t687 + Ifges(6,4) * t686 + Ifges(6,5) * t705 + t784 * t657 - t790 * t658 + t721 * t692 - t740 * t693;
t694 = Ifges(6,1) * t722 + Ifges(6,4) * t721 + Ifges(6,5) * t740;
t653 = -mrSges(6,1) * t679 + mrSges(6,3) * t671 + Ifges(6,4) * t687 + Ifges(6,2) * t686 + Ifges(6,6) * t705 - t692 * t722 + t694 * t740 + t803;
t713 = Ifges(5,4) * t742 + Ifges(5,2) * t741 + Ifges(5,6) * t755;
t714 = Ifges(5,1) * t742 + Ifges(5,4) * t741 + Ifges(5,5) * t755;
t802 = -mrSges(5,1) * t679 - mrSges(5,2) * t680 + Ifges(5,5) * t707 + Ifges(5,6) * t706 + Ifges(5,3) * t727 + t785 * t645 + t791 * t653 + t742 * t713 - t714 * t741;
t801 = m(4) * t746 - t728 * mrSges(4,1) + t729 * mrSges(4,2) - t759 * t747 + t760 * t748 + t641;
t712 = Ifges(5,5) * t742 + Ifges(5,6) * t741 + Ifges(5,3) * t755;
t637 = mrSges(5,2) * t698 + mrSges(5,3) * t679 + Ifges(5,1) * t707 + Ifges(5,4) * t706 + Ifges(5,5) * t727 + t645 * t791 - t653 * t785 + t712 * t741 - t713 * t755;
t798 = pkin(6) * t656 - mrSges(6,1) * t670 + mrSges(6,2) * t671 - Ifges(6,5) * t687 - Ifges(6,6) * t686 - Ifges(6,3) * t705 - t790 * t657 - t784 * t658 - t722 * t693 + t721 * t694;
t643 = -mrSges(5,1) * t698 + mrSges(5,3) * t680 + Ifges(5,4) * t707 + Ifges(5,2) * t706 + Ifges(5,6) * t727 - t742 * t712 + t755 * t714 + t798;
t735 = Ifges(4,4) * t760 + Ifges(4,2) * t759 + Ifges(4,6) * t781;
t736 = Ifges(4,1) * t760 + Ifges(4,4) * t759 + Ifges(4,5) * t781;
t800 = pkin(4) * t804 - pkin(5) * t642 + mrSges(4,1) * t730 - mrSges(4,2) * t731 + Ifges(4,5) * t729 + Ifges(4,6) * t728 + Ifges(4,3) * t780 - t786 * t637 - t792 * t643 + t760 * t735 - t759 * t736;
t799 = m(3) * t783 - t769 * mrSges(3,1) + t770 * mrSges(3,2) + t773 * t810 + t774 * t809 + t801;
t757 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t794 - Ifges(3,2) * t788) * qJD(1);
t758 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t794 - Ifges(3,4) * t788) * qJD(1);
t797 = pkin(3) * t811 + mrSges(3,1) * t772 - mrSges(3,2) * t771 + Ifges(3,5) * t770 + Ifges(3,6) * t769 + Ifges(3,3) * qJDD(2) + t757 * t809 + t758 * t810 + t800;
t795 = cos(qJ(1));
t789 = sin(qJ(1));
t756 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t794 - Ifges(3,6) * t788) * qJD(1);
t734 = Ifges(4,5) * t760 + Ifges(4,6) * t759 + Ifges(4,3) * t781;
t638 = qJDD(1) * mrSges(2,1) - mrSges(2,2) * t796 + t799;
t632 = -pkin(4) * t641 - mrSges(4,1) * t746 + mrSges(4,3) * t731 + Ifges(4,4) * t729 + Ifges(4,2) * t728 + Ifges(4,6) * t780 - t734 * t760 + t736 * t781 + t802;
t631 = pkin(5) * t641 + mrSges(4,2) * t746 - mrSges(4,3) * t730 + Ifges(4,1) * t729 + Ifges(4,4) * t728 + Ifges(4,5) * t780 + t637 * t792 - t643 * t786 + t734 * t759 - t735 * t781;
t630 = -mrSges(2,1) * t796 - qJDD(1) * mrSges(2,2) + t634 * t788 + t635 * t794;
t629 = t630 * t795 - t638 * t789;
t628 = t630 * t789 + t638 * t795;
t627 = mrSges(3,2) * t783 - mrSges(3,3) * t772 + Ifges(3,1) * t770 + Ifges(3,4) * t769 + Ifges(3,5) * qJDD(2) - qJD(2) * t757 + t631 * t793 - t632 * t787 - t756 * t810;
t626 = -pkin(3) * t801 - mrSges(3,1) * t783 + mrSges(3,3) * t771 + Ifges(3,4) * t770 + Ifges(3,2) * t769 + Ifges(3,6) * qJDD(2) + qJD(2) * t758 + t787 * t631 + t793 * t632 - t756 * t809;
t625 = -pkin(2) * t806 + mrSges(2,1) * g(3) + t796 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + t797;
t624 = -mrSges(2,2) * g(3) + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t796 + t626 * t794 + t627 * t788;
t623 = pkin(2) * t799 + Ifges(2,3) * qJDD(1) - t788 * t626 + t794 * t627;
t1 = [t629; t628; (-m(1) - m(2)) * g(3) + t806; -pkin(1) * t628 - mrSges(1,2) * g(3) + t624 * t795 - t625 * t789; pkin(1) * t629 + mrSges(1,1) * g(3) + t624 * t789 + t625 * t795; t623; t623; t797; t800; t802; -t798; t803;];
tauJB = t1;
