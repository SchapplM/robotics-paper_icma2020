% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% CloosQRC350DE
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = CloosQRC350DE_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
JRD_rot=NaN(9,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t32 = qJD(1) * sin(qJ(1));
	t31 = qJD(1) * cos(qJ(1));
	t1 = [-t31, 0, 0, 0, 0, 0; t32, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t32, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t31 = sin(qJ(1));
	t38 = qJD(1) * t31;
	t33 = cos(qJ(1));
	t37 = qJD(1) * t33;
	t30 = sin(qJ(2));
	t36 = qJD(2) * t30;
	t32 = cos(qJ(2));
	t35 = qJD(2) * t32;
	t34 = qJD(2) * t33;
	t29 = t31 * t36 - t32 * t37;
	t28 = t30 * t37 + t31 * t35;
	t27 = t30 * t34 + t32 * t38;
	t26 = t30 * t38 - t32 * t34;
	t1 = [-t28, -t27, 0, 0, 0, 0; t26, t29, 0, 0, 0, 0; 0, -t35, 0, 0, 0, 0; t29, t26, 0, 0, 0, 0; t27, t28, 0, 0, 0, 0; 0, t36, 0, 0, 0, 0; -t38, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (60->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t70 = qJ(2) + qJ(3);
	t68 = cos(t70);
	t69 = qJD(2) + qJD(3);
	t77 = t69 * t68;
	t71 = sin(qJ(1));
	t76 = t69 * t71;
	t72 = cos(qJ(1));
	t75 = t69 * t72;
	t74 = qJD(1) * t71;
	t73 = qJD(1) * t72;
	t67 = sin(t70);
	t66 = t69 * t67;
	t65 = t67 * t76 - t68 * t73;
	t64 = t67 * t73 + t68 * t76;
	t63 = t67 * t75 + t68 * t74;
	t62 = t67 * t74 - t68 * t75;
	t1 = [-t64, -t63, -t63, 0, 0, 0; t62, t65, t65, 0, 0, 0; 0, -t77, -t77, 0, 0, 0; t65, t62, t62, 0, 0, 0; t63, t64, t64, 0, 0, 0; 0, t66, t66, 0, 0, 0; -t74, 0, 0, 0, 0, 0; -t73, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:01
	% EndTime: 2020-06-23 21:15:02
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (164->29), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t326 = cos(qJ(4));
	t323 = qJ(2) + qJ(3);
	t320 = sin(t323);
	t331 = qJD(4) * t320 + qJD(1);
	t347 = t326 * t331;
	t324 = sin(qJ(4));
	t346 = t331 * t324;
	t322 = qJD(2) + qJD(3);
	t345 = t322 * t324;
	t325 = sin(qJ(1));
	t344 = t322 * t325;
	t343 = t322 * t326;
	t327 = cos(qJ(1));
	t342 = t322 * t327;
	t341 = qJD(1) * t325;
	t340 = qJD(1) * t327;
	t339 = qJD(4) * t324;
	t338 = qJD(4) * t326;
	t337 = qJD(4) * t327;
	t321 = cos(t323);
	t336 = t321 * t343;
	t335 = t320 * t344;
	t334 = t321 * t344;
	t333 = t320 * t342;
	t332 = t321 * t342;
	t330 = qJD(1) * t320 + qJD(4);
	t329 = t330 * t327;
	t328 = t325 * t330 - t332;
	t319 = t322 * t320;
	t318 = t320 * t340 + t334;
	t317 = t320 * t341 - t332;
	t316 = t320 * t339 - t336;
	t315 = t320 * t338 + t321 * t345;
	t314 = t326 * t335 + (t325 * t339 - t326 * t340) * t321;
	t313 = -t324 * t335 + (t324 * t340 + t325 * t338) * t321;
	t312 = -t326 * t333 + (-t324 * t337 - t326 * t341) * t321;
	t311 = t324 * t333 + (t324 * t341 - t326 * t337) * t321;
	t310 = t326 * t329 + (t336 - t346) * t325;
	t309 = t325 * t347 + (t329 + t334) * t324;
	t308 = t326 * t328 + t327 * t346;
	t307 = t324 * t328 - t327 * t347;
	t1 = [-t310, t312, t312, t307, 0, 0; t308, t314, t314, t309, 0, 0; 0, t316, t316, t320 * t345 - t321 * t338, 0, 0; t309, t311, t311, t308, 0, 0; -t307, t313, t313, t310, 0, 0; 0, t315, t315, t320 * t343 + t321 * t339, 0, 0; -t321 * t340 + t335, t317, t317, 0, 0, 0; t321 * t341 + t333, t318, t318, 0, 0, 0; 0, t319, t319, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:04
	% EndTime: 2020-06-23 21:15:05
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (454->63), mult. (696->118), div. (0->0), fcn. (708->8), ass. (0->71)
	t462 = qJ(2) + qJ(3);
	t460 = cos(t462);
	t463 = sin(qJ(5));
	t466 = cos(qJ(5));
	t461 = qJD(2) + qJD(3);
	t467 = cos(qJ(4));
	t486 = qJD(5) * t467 + t461;
	t464 = sin(qJ(4));
	t499 = qJD(4) * t464;
	t474 = t486 * t463 + t466 * t499;
	t515 = t474 * t460;
	t475 = -t463 * t499 + t486 * t466;
	t514 = t475 * t460;
	t459 = sin(t462);
	t468 = cos(qJ(1));
	t503 = t468 * t464;
	t491 = t459 * t503;
	t465 = sin(qJ(1));
	t504 = t465 * t467;
	t478 = t491 + t504;
	t479 = t459 * t504 + t503;
	t506 = t461 * t468;
	t494 = t460 * t506;
	t445 = t479 * qJD(1) + t478 * qJD(4) - t467 * t494;
	t502 = t468 * t467;
	t505 = t465 * t464;
	t454 = t459 * t502 - t505;
	t501 = qJD(1) * t465;
	t490 = t460 * t501;
	t497 = qJD(5) * t460;
	t513 = (t468 * t497 - t445) * t463 + (qJD(5) * t454 + t459 * t506 + t490) * t466;
	t512 = t459 * t463;
	t511 = t459 * t466;
	t510 = t460 * t467;
	t509 = t461 * t464;
	t508 = t461 * t465;
	t507 = t461 * t467;
	t500 = qJD(1) * t468;
	t498 = qJD(4) * t467;
	t496 = qJD(5) * t463;
	t495 = qJD(5) * t466;
	t493 = t461 * t512;
	t492 = t459 * t505;
	t487 = qJD(5) + t507;
	t485 = qJD(4) * t459 + qJD(1);
	t484 = qJD(1) * t459 + qJD(4);
	t447 = -t484 * t502 + (-t460 * t507 + t485 * t464) * t465;
	t482 = t465 * t497 - t447;
	t481 = t487 * t466;
	t480 = t487 * t463;
	t477 = qJD(1) * (-t466 * t510 + t512);
	t476 = qJD(1) * (t463 * t510 + t511);
	t473 = qJD(5) * t479 + t459 * t508 - t460 * t500;
	t471 = t445 * t466 + t468 * t493 + t463 * t490 + (-t460 * t466 * t468 + t454 * t463) * qJD(5);
	t470 = t459 * t480 - t514;
	t469 = t487 * t511 + t515;
	t451 = t492 - t502;
	t450 = -t459 * t498 - t460 * t509;
	t449 = t461 * t492 + (-t464 * t500 - t465 * t498) * t460;
	t448 = -t461 * t491 + (-t464 * t501 + t468 * t498) * t460;
	t446 = t485 * t504 + (t460 * t508 + t484 * t468) * t464;
	t444 = -t485 * t502 + (t484 * t465 - t494) * t464;
	t443 = t474 * t459 - t460 * t481;
	t442 = t475 * t459 + t460 * t480;
	t441 = t469 * t465 + t468 * t477;
	t440 = t468 * t476 + (-t487 * t512 + t514) * t465;
	t439 = t465 * t477 + (-t459 * t481 - t515) * t468;
	t438 = t465 * t476 + t470 * t468;
	t437 = t473 * t463 - t482 * t466;
	t436 = t482 * t463 + t473 * t466;
	t1 = [t437, t439, t439, t444 * t466 + t478 * t496, -t513, 0; t471, t441, t441, t446 * t466 - t451 * t496, t436, 0; 0, t443, t443, t509 * t511 + (t464 * t496 - t466 * t498) * t460, t470, 0; t436, t438, t438, -t444 * t463 + t478 * t495, t471, 0; t513, t440, t440, -t446 * t463 - t451 * t495, -t437, 0; 0, t442, t442, -t464 * t493 + (t463 * t498 + t464 * t495) * t460, t469, 0; -t446, t448, t448, -t445, 0, 0; t444, t449, t449, t447, 0, 0; 0, t450, t450, -t459 * t507 - t460 * t499, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:09
	% EndTime: 2020-06-23 21:15:11
	% DurationCPUTime: 2.01s
	% Computational Cost: add. (1531->123), mult. (2288->211), div. (0->0), fcn. (2102->10), ass. (0->114)
	t724 = pkin(7) * qJ(5) - qJ(6);
	t719 = sin(t724);
	t720 = cos(t724);
	t728 = qJ(2) + qJ(3);
	t725 = sin(t728);
	t727 = qJD(2) + qJD(3);
	t733 = cos(qJ(4));
	t763 = qJD(1) * t725 + qJD(4);
	t734 = cos(qJ(1));
	t787 = qJD(4) * t733;
	t774 = t734 * t787;
	t726 = cos(t728);
	t730 = sin(qJ(4));
	t792 = t734 * t730;
	t781 = t726 * t792;
	t789 = qJD(1) * t734;
	t731 = sin(qJ(1));
	t795 = t731 * t730;
	t695 = -t725 * t774 - t727 * t781 - t733 * t789 + t763 * t795;
	t791 = t734 * t733;
	t712 = t725 * t791 - t795;
	t723 = pkin(7) * qJD(5) - qJD(6);
	t729 = sin(qJ(5));
	t732 = cos(qJ(5));
	t800 = t726 * t734;
	t770 = (t712 * t732 + t729 * t800) * t723 + t695;
	t778 = t725 * t792;
	t794 = t731 * t733;
	t711 = t778 + t794;
	t750 = t725 * t794 + t792;
	t797 = t727 * t733;
	t782 = t726 * t797;
	t696 = t750 * qJD(1) + t711 * qJD(4) - t734 * t782;
	t701 = -t712 * t729 + t732 * t800;
	t790 = qJD(1) * t731;
	t777 = t726 * t790;
	t804 = t725 * t729;
	t780 = t727 * t804;
	t684 = t701 * qJD(5) - t696 * t732 - t729 * t777 - t734 * t780;
	t819 = t711 * t723 + t684;
	t812 = t770 * t719 - t720 * t819;
	t813 = t719 * t819 + t770 * t720;
	t766 = qJD(5) + t797;
	t803 = t725 * t732;
	t765 = qJD(5) * t733 + t727;
	t788 = qJD(4) * t732;
	t740 = t765 * t729 + t730 * t788;
	t821 = t740 * t726;
	t692 = t766 * t803 + t821;
	t807 = t723 * t730;
	t762 = t726 * t807 - t692;
	t793 = t732 * t733;
	t708 = t726 * t793 - t804;
	t823 = t708 * t723;
	t806 = t723 * t732;
	t822 = t730 * (-qJD(4) + t806);
	t741 = -qJD(4) * t729 * t730 + t765 * t732;
	t820 = t741 * t726;
	t818 = (t730 * t790 - t774) * t726 + t727 * t778 + t734 * t823;
	t779 = t725 * t795;
	t817 = (t730 * t789 + t731 * t787) * t726 - t727 * t779 - t731 * t823;
	t799 = t727 * t730;
	t737 = -t725 * t799 + t726 * t787 - t823;
	t816 = t762 * t719 - t737 * t720;
	t815 = t737 * t719 + t762 * t720;
	t764 = qJD(4) * t725 + qJD(1);
	t798 = t727 * t731;
	t697 = t764 * t794 + (t726 * t798 + t763 * t734) * t730;
	t801 = t726 * t731;
	t768 = (-t729 * t801 - t732 * t750) * t723 + t697;
	t736 = qJD(5) * t750 + t725 * t798 - t726 * t789;
	t698 = -t763 * t791 + (t764 * t730 - t782) * t731;
	t786 = qJD(5) * t726;
	t756 = t731 * t786 - t698;
	t686 = t736 * t729 - t756 * t732;
	t709 = t779 - t791;
	t771 = t709 * t723 - t686;
	t674 = t771 * t719 - t768 * t720;
	t814 = t768 * t719 + t771 * t720;
	t805 = t725 * t727;
	t811 = (t734 * t786 - t696) * t729 + (qJD(5) * t712 + t734 * t805 + t777) * t732;
	t809 = t719 * t723;
	t808 = t720 * t723;
	t802 = t726 * t729;
	t796 = t730 * t732;
	t785 = qJD(5) * t729;
	t784 = qJD(5) * t732;
	t754 = t766 * t732;
	t761 = t726 * t754 + (-t740 + t807) * t725;
	t760 = t711 * t806 - t696;
	t759 = t709 * t806 - t698;
	t755 = t729 * t766;
	t745 = qJD(1) * t708;
	t753 = -t723 * t781 + t731 * t745 - (-t725 * t754 - t821) * t734;
	t752 = t762 * t731 + t734 * t745;
	t751 = t733 * t802 + t803;
	t744 = qJD(1) * t751;
	t743 = -t695 * t732 - t711 * t785 - t712 * t723;
	t742 = -t697 * t732 + t709 * t785 + t723 * t750;
	t739 = -t730 * t785 + (-t723 + t788) * t733;
	t738 = (-t725 * t793 - t802) * t723 + t725 * t787 + t726 * t799;
	t691 = t725 * t755 - t820;
	t699 = t729 * t750 - t732 * t801;
	t693 = t741 * t725 + t726 * t755;
	t689 = t734 * t744 + (-t766 * t804 + t820) * t731;
	t687 = t691 * t734 + t731 * t744;
	t685 = t756 * t729 + t736 * t732;
	t681 = t738 * t719 + t761 * t720;
	t680 = t761 * t719 - t738 * t720;
	t679 = t817 * t719 + t752 * t720;
	t678 = t752 * t719 - t817 * t720;
	t677 = t818 * t719 + t753 * t720;
	t676 = t753 * t719 - t818 * t720;
	t1 = [t814, t677, t677, -t760 * t719 + t743 * t720, t813 * pkin(7) + t701 * t809 + t720 * t811, -t813; -t812, t679, t679, t759 * t719 + t742 * t720, -t674 * pkin(7) - t685 * t720 + t699 * t809, t674; 0, t681, t681, (t719 * t733 - t720 * t796) * t805 + (-t719 * t822 + t739 * t720) * t726, t816 * pkin(7) - t691 * t720 - t751 * t809, -t816; t674, t676, t676, t743 * t719 + t760 * t720, t812 * pkin(7) - t701 * t808 + t719 * t811, -t812; t813, t678, t678, t742 * t719 - t759 * t720, t814 * pkin(7) - t685 * t719 - t699 * t808, -t814; 0, t680, t680, (-t719 * t796 - t720 * t733) * t805 + (t739 * t719 + t720 * t822) * t726, -t815 * pkin(7) - t691 * t719 + t751 * t808, t815; t685, t687, t687, -t695 * t729 + t711 * t784, -t684, 0; t811, t689, t689, -t697 * t729 - t709 * t784, -t686, 0; 0, t693, t693, -t730 * t780 + (t729 * t787 + t730 * t784) * t726, t692, 0;];
	JRD_rot = t1;
end