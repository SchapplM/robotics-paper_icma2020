% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = CloosQRC350OL_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350OL_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
JRD_rot=NaN(9,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(2));
	t37 = qJD(2) * t31;
	t33 = cos(qJ(2));
	t36 = qJD(2) * t33;
	t35 = qJD(2) * t34;
	t30 = -t32 * t37 + t33 * t38;
	t29 = -t31 * t38 - t32 * t36;
	t28 = -t31 * t35 - t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t29, t28, 0, 0, 0, 0; -t27, t30, 0, 0, 0, 0; 0, -t36, 0, 0, 0, 0; -t30, t27, 0, 0, 0, 0; t28, t29, 0, 0, 0, 0; 0, t37, 0, 0, 0, 0; t39, 0, 0, 0, 0, 0; -t38, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (59->11), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t72 = qJ(2) + qJ(3);
	t70 = cos(t72);
	t71 = qJD(2) + qJD(3);
	t79 = t71 * t70;
	t73 = sin(qJ(1));
	t78 = t71 * t73;
	t74 = cos(qJ(1));
	t77 = t71 * t74;
	t76 = qJD(1) * t73;
	t75 = qJD(1) * t74;
	t69 = sin(t72);
	t68 = t71 * t69;
	t67 = -t69 * t78 + t70 * t75;
	t66 = -t69 * t75 - t70 * t78;
	t65 = -t69 * t77 - t70 * t76;
	t64 = t69 * t76 - t70 * t77;
	t1 = [t66, t65, t65, 0, 0, 0; -t64, t67, t67, 0, 0, 0; 0, -t79, -t79, 0, 0, 0; -t67, t64, t64, 0, 0, 0; t65, t66, t66, 0, 0, 0; 0, t68, t68, 0, 0, 0; t76, 0, 0, 0, 0, 0; -t75, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:52
	% EndTime: 2020-06-23 22:04:53
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (164->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t313 = qJD(2) + qJD(3);
	t315 = sin(qJ(4));
	t337 = t313 * t315;
	t316 = sin(qJ(1));
	t336 = t313 * t316;
	t317 = cos(qJ(4));
	t335 = t313 * t317;
	t318 = cos(qJ(1));
	t334 = t313 * t318;
	t333 = t317 * t318;
	t332 = qJD(1) * t316;
	t331 = qJD(1) * t318;
	t330 = qJD(4) * t315;
	t329 = qJD(4) * t317;
	t328 = qJD(4) * t318;
	t314 = qJ(2) + qJ(3);
	t312 = cos(t314);
	t327 = t312 * t335;
	t311 = sin(t314);
	t326 = t311 * t336;
	t325 = t312 * t336;
	t324 = t311 * t334;
	t323 = t312 * t334;
	t322 = qJD(4) * t311 - qJD(1);
	t321 = qJD(1) * t311 - qJD(4);
	t320 = t322 * t315;
	t319 = t321 * t316 - t323;
	t310 = t313 * t311;
	t309 = -t311 * t331 - t325;
	t308 = t311 * t332 - t323;
	t307 = t311 * t330 - t327;
	t306 = t311 * t329 + t312 * t337;
	t305 = -t317 * t326 + (-t316 * t330 + t317 * t331) * t312;
	t304 = t315 * t326 + (-t315 * t331 - t316 * t329) * t312;
	t303 = -t317 * t324 + (-t315 * t328 - t317 * t332) * t312;
	t302 = t315 * t324 + (t315 * t332 - t317 * t328) * t312;
	t301 = -t321 * t333 + (t320 - t327) * t316;
	t300 = t322 * t317 * t316 + (t321 * t318 + t325) * t315;
	t299 = t319 * t317 + t318 * t320;
	t298 = t319 * t315 - t322 * t333;
	t1 = [t301, t303, t303, t298, 0, 0; -t299, t305, t305, -t300, 0, 0; 0, t307, t307, t311 * t337 - t312 * t329, 0, 0; t300, t302, t302, t299, 0, 0; t298, t304, t304, t301, 0, 0; 0, t306, t306, t311 * t335 + t312 * t330, 0, 0; -t312 * t331 + t326, t308, t308, 0, 0, 0; -t312 * t332 - t324, t309, t309, 0, 0, 0; 0, t310, t310, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:55
	% EndTime: 2020-06-23 22:04:56
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (454->61), mult. (696->114), div. (0->0), fcn. (708->8), ass. (0->73)
	t464 = qJ(2) + qJ(3);
	t462 = cos(t464);
	t465 = sin(qJ(5));
	t468 = cos(qJ(5));
	t463 = qJD(2) + qJD(3);
	t469 = cos(qJ(4));
	t486 = qJD(5) * t469 + t463;
	t466 = sin(qJ(4));
	t501 = qJD(4) * t466;
	t476 = t465 * t486 + t468 * t501;
	t518 = t462 * t476;
	t467 = sin(qJ(1));
	t461 = sin(t464);
	t484 = qJD(1) * t461 - qJD(4);
	t470 = cos(qJ(1));
	t508 = t463 * t470;
	t517 = -t462 * t508 + t467 * t484;
	t516 = -t465 * t501 + t468 * t486;
	t504 = t470 * t469;
	t507 = t467 * t466;
	t456 = t461 * t504 + t507;
	t500 = qJD(4) * t469;
	t488 = t470 * t500;
	t493 = t461 * t507;
	t510 = t463 * t467;
	t496 = t462 * t510;
	t449 = qJD(1) * t456 - qJD(4) * t493 + t469 * t496 - t488;
	t505 = t470 * t466;
	t506 = t467 * t469;
	t454 = t461 * t506 - t505;
	t502 = qJD(1) * t470;
	t491 = t462 * t502;
	t499 = qJD(5) * t462;
	t515 = t465 * (t467 * t499 + t449) - t468 * (-qJD(5) * t454 - t461 * t510 + t491);
	t514 = t461 * t465;
	t513 = t461 * t468;
	t512 = t462 * t469;
	t511 = t463 * t466;
	t509 = t463 * t469;
	t503 = qJD(1) * t467;
	t498 = qJD(5) * t465;
	t497 = qJD(5) * t468;
	t494 = t463 * t514;
	t492 = t461 * t505;
	t487 = qJD(5) + t509;
	t485 = -qJD(4) * t461 + qJD(1);
	t479 = t485 * t470;
	t447 = t466 * t479 - t517 * t469;
	t483 = t470 * t499 + t447;
	t481 = t487 * t468;
	t480 = t487 * t465;
	t478 = t468 * t512 - t514;
	t477 = t465 * t512 + t513;
	t474 = -qJD(5) * t456 - t461 * t508 - t462 * t503;
	t473 = -t465 * t491 - t449 * t468 + t467 * t494 + (-t462 * t467 * t468 + t454 * t465) * qJD(5);
	t472 = -t461 * t481 - t518;
	t471 = t461 * t480 - t516 * t462;
	t455 = -t492 + t506;
	t453 = -t493 - t504;
	t452 = -t461 * t500 - t462 * t511;
	t451 = -t463 * t493 + (t466 * t502 + t467 * t500) * t462;
	t450 = -t463 * t492 + (-t466 * t503 + t488) * t462;
	t448 = t485 * t506 + (-t470 * t484 - t496) * t466;
	t446 = t517 * t466 + t469 * t479;
	t445 = t461 * t476 - t462 * t481;
	t444 = t516 * t461 + t462 * t480;
	t443 = t467 * t472 + t478 * t502;
	t442 = t467 * t471 - t477 * t502;
	t441 = t470 * t472 - t478 * t503;
	t440 = t470 * t471 + t477 * t503;
	t439 = t465 * t474 + t468 * t483;
	t438 = -t465 * t483 + t468 * t474;
	t1 = [t473, t441, t441, t446 * t468 - t455 * t498, t438, 0; t439, t443, t443, t448 * t468 - t453 * t498, -t515, 0; 0, t445, t445, t511 * t513 + (t466 * t498 - t468 * t500) * t462, t471, 0; t515, t440, t440, -t446 * t465 - t455 * t497, -t439, 0; t438, t442, t442, -t448 * t465 - t453 * t497, t473, 0; 0, t444, t444, -t466 * t494 + (t465 * t500 + t466 * t497) * t462, t487 * t513 + t518, 0; t448, t450, t450, t447, 0, 0; -t446, t451, t451, t449, 0, 0; 0, t452, t452, -t461 * t509 - t462 * t501, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:59
	% EndTime: 2020-06-23 22:05:01
	% DurationCPUTime: 1.76s
	% Computational Cost: add. (1085->118), mult. (1726->211), div. (0->0), fcn. (1804->10), ass. (0->111)
	t685 = qJ(2) + qJ(3);
	t682 = sin(t685);
	t692 = cos(qJ(4));
	t693 = cos(qJ(1));
	t745 = t693 * t692;
	t688 = sin(qJ(4));
	t689 = sin(qJ(1));
	t750 = t689 * t688;
	t671 = t682 * t745 + t750;
	t739 = qJD(4) * t693;
	t721 = t692 * t739;
	t742 = qJD(4) * t688;
	t722 = t689 * t742;
	t683 = cos(t685);
	t684 = qJD(2) + qJD(3);
	t753 = t684 * t689;
	t731 = t683 * t753;
	t657 = t671 * qJD(1) - t682 * t722 + t692 * t731 - t721;
	t746 = t693 * t688;
	t749 = t689 * t692;
	t669 = t682 * t749 - t746;
	t687 = sin(qJ(5));
	t691 = cos(qJ(5));
	t756 = t683 * t689;
	t658 = -t669 * t687 + t691 * t756;
	t743 = qJD(1) * t693;
	t726 = t683 * t743;
	t759 = t682 * t687;
	t729 = t684 * t759;
	t644 = t658 * qJD(5) + t657 * t691 + t687 * t726 - t689 * t729;
	t740 = qJD(4) * t692;
	t705 = t688 * t743 + t689 * t740;
	t744 = qJD(1) * t689;
	t656 = t705 * t682 - t692 * t744 + (t731 - t739) * t688;
	t659 = t669 * t691 + t687 * t756;
	t728 = t682 * t750;
	t668 = t728 + t745;
	t686 = sin(qJ(6));
	t690 = cos(qJ(6));
	t773 = t644 * t690 - t656 * t686 + (-t659 * t686 - t668 * t690) * qJD(6);
	t772 = t644 * t686 + t656 * t690 + (t659 * t690 - t668 * t686) * qJD(6);
	t720 = t684 * t692 + qJD(5);
	t710 = t720 * t691;
	t734 = qJD(6) * t688;
	t723 = t683 * t734;
	t719 = qJD(5) * t692 + t684;
	t741 = qJD(4) * t691;
	t701 = t719 * t687 + t688 * t741;
	t768 = t701 * t683;
	t771 = t682 * t710 + t723 + t768;
	t747 = t691 * t692;
	t667 = t683 * t747 - t759;
	t767 = qJD(6) * t667;
	t766 = -t687 * t742 + t719 * t691;
	t738 = qJD(5) * t683;
	t765 = (t689 * t738 + t657) * t687 - (-qJD(5) * t669 - t682 * t753 + t726) * t691;
	t760 = t682 * t684;
	t758 = t682 * t691;
	t757 = t683 * t687;
	t755 = t683 * t693;
	t754 = t684 * t688;
	t752 = t684 * t693;
	t751 = t688 * t691;
	t748 = t690 * t692;
	t737 = qJD(5) * t687;
	t736 = qJD(5) * t691;
	t735 = qJD(6) * t686;
	t733 = qJD(6) * t690;
	t732 = qJD(6) * t691;
	t730 = t683 * t752;
	t727 = t682 * t746;
	t718 = -qJD(4) - t732;
	t717 = qJD(6) + t741;
	t655 = (-qJD(4) * t682 + qJD(1)) * t746 + (t730 + (-qJD(1) * t682 + qJD(4)) * t689) * t692;
	t716 = t693 * t738 + t655;
	t651 = t720 * t758 + t768;
	t714 = -t651 - t723;
	t713 = -t683 * t710 + (t701 + t734) * t682;
	t670 = t727 - t749;
	t712 = -t670 * t732 + t655;
	t711 = -t668 * t732 + t657;
	t709 = t720 * t687;
	t708 = t667 * t744 + t771 * t693;
	t664 = t667 * t693;
	t707 = -qJD(1) * t664 + t771 * t689;
	t706 = t692 * t757 + t758;
	t654 = t668 * qJD(1) - t682 * t721 - t688 * t730 - t722;
	t703 = -qJD(6) * t671 + t654 * t691 + t670 * t737;
	t702 = -qJD(6) * t669 - t656 * t691 + t668 * t737;
	t699 = -qJD(5) * t671 - t682 * t752 - t683 * t744;
	t698 = qJD(6) * (-t682 * t747 - t757) - t682 * t740 - t683 * t754;
	t697 = -t682 * t754 + t683 * t740 + t767;
	t696 = t705 * t683 - t684 * t728 + t689 * t767;
	t695 = -t684 * t727 + qJD(6) * t664 + (-t688 * t744 + t721) * t683;
	t650 = t682 * t709 - t766 * t683;
	t662 = t671 * t691 + t687 * t755;
	t661 = -t671 * t687 + t691 * t755;
	t652 = t766 * t682 + t683 * t709;
	t648 = t650 * t689 - t706 * t743;
	t646 = t650 * t693 + t706 * t744;
	t642 = t699 * t687 + t716 * t691;
	t641 = -t716 * t687 + t699 * t691;
	t640 = t698 * t686 - t713 * t690;
	t639 = t713 * t686 + t698 * t690;
	t638 = t686 * t696 + t690 * t707;
	t637 = -t686 * t707 + t690 * t696;
	t636 = t695 * t686 + t690 * t708;
	t635 = -t686 * t708 + t695 * t690;
	t634 = t642 * t690 + t654 * t686 + (-t662 * t686 - t670 * t690) * qJD(6);
	t633 = t642 * t686 - t654 * t690 + (t662 * t690 - t670 * t686) * qJD(6);
	t1 = [t773, t636, t636, t712 * t686 - t703 * t690, -t641 * t690 + t661 * t735, t633; -t634, t638, t638, t711 * t686 - t702 * t690, t658 * t735 + t690 * t765, t772; 0, t640, t640, (-t686 * t692 - t690 * t751) * t760 + (t717 * t748 + (t718 * t686 - t690 * t737) * t688) * t683, -t650 * t690 - t706 * t735, t714 * t686 + t697 * t690; -t772, t635, t635, t703 * t686 + t712 * t690, t641 * t686 + t661 * t733, t634; t633, t637, t637, t702 * t686 + t711 * t690, t658 * t733 - t686 * t765, t773; 0, t639, t639, (t686 * t751 - t748) * t760 + (t718 * t690 * t688 + (t688 * t737 - t717 * t692) * t686) * t683, t650 * t686 - t706 * t733, -t697 * t686 + t714 * t690; t765, t646, t646, -t654 * t687 + t670 * t736, -t642, 0; t641, t648, t648, t656 * t687 + t668 * t736, -t644, 0; 0, t652, t652, -t688 * t729 + (t687 * t740 + t688 * t736) * t683, t651, 0;];
	JRD_rot = t1;
end