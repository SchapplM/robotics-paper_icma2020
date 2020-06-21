% Calculate kinetic energy for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
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
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = CloosQRC350DE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:24
% EndTime: 2020-06-19 21:37:24
% DurationCPUTime: 0.65s
% Computational Cost: add. (2759->478), mult. (3590->612), div. (0->0), fcn. (3784->12), ass. (0->199)
t1 = V_base(5) * pkin(1);
t2 = V_base(6) - qJD(1);
t3 = sin(qJ(1));
t5 = t2 * t3 * pkin(2);
t6 = cos(qJ(1));
t7 = qJD(2) * t6;
t8 = V_base(5) + t7;
t9 = cos(qJ(2));
t11 = sin(qJ(2));
t13 = rSges(3,1) * t9 - rSges(3,2) * t11;
t15 = t3 * t11;
t17 = t3 * t9;
t20 = -rSges(3,1) * t15 - rSges(3,2) * t17 + rSges(3,3) * t6;
t23 = (t8 * t13 - t2 * t20 + t1 + t5 + V_base(1)) ^ 2;
t24 = V_base(4) * pkin(1);
t26 = t2 * t6 * pkin(2);
t27 = qJD(2) * t3;
t28 = V_base(4) + t27;
t30 = t6 * t11;
t32 = t6 * t9;
t35 = rSges(3,1) * t30 + rSges(3,2) * t32 + rSges(3,3) * t3;
t38 = (-t28 * t13 + t2 * t35 - t24 + t26 + V_base(2)) ^ 2;
t40 = V_base(4) * t3 * pkin(2);
t42 = V_base(5) * t6 * pkin(2);
t46 = (t28 * t20 - t8 * t35 - t40 - t42 + V_base(3)) ^ 2;
t50 = Icges(2,4) * t3;
t51 = Icges(2,1) * t6 + t50;
t53 = Icges(2,4) * t6;
t55 = Icges(2,2) * t3 + t53;
t60 = -Icges(2,1) * t3 + t53;
t63 = Icges(2,2) * t6 - t50;
t69 = Icges(2,5) * t6 + Icges(2,6) * t3;
t83 = -Icges(2,5) * t3 + Icges(2,6) * t6;
t90 = -rSges(2,1) * t3 + rSges(2,2) * t6;
t93 = (V_base(5) * rSges(2,3) - t2 * t90 + t1 + V_base(1)) ^ 2;
t97 = rSges(2,1) * t6 + rSges(2,2) * t3;
t100 = (-V_base(4) * rSges(2,3) + t2 * t97 - t24 + V_base(2)) ^ 2;
t104 = (t90 * V_base(4) - t97 * V_base(5) + V_base(3)) ^ 2;
t130 = (-rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1)) ^ 2;
t134 = (rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2)) ^ 2;
t138 = (-rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3)) ^ 2;
t142 = t8 * t9 * pkin(3);
t145 = pkin(3) * t11 + pkin(2);
t147 = pkin(2) * t3 - t3 * t145;
t148 = t2 * t147;
t149 = qJD(3) * t6;
t150 = V_base(5) + t7 + t149;
t151 = qJ(2) + qJ(3);
t152 = cos(t151);
t154 = sin(t151);
t156 = rSges(4,1) * t152 - rSges(4,2) * t154;
t158 = t3 * t154;
t160 = t3 * t152;
t163 = -rSges(4,1) * t158 - rSges(4,2) * t160 + rSges(4,3) * t6;
t166 = (t150 * t156 - t2 * t163 + t1 + t142 - t148 + t5 + V_base(1)) ^ 2;
t168 = t28 * t9 * pkin(3);
t171 = -pkin(2) * t6 + t6 * t145;
t172 = t2 * t171;
t173 = qJD(3) * t3;
t174 = V_base(4) + t27 + t173;
t176 = t6 * t154;
t178 = t6 * t152;
t181 = rSges(4,1) * t176 + rSges(4,2) * t178 + rSges(4,3) * t3;
t184 = (-t174 * t156 + t2 * t181 - t168 + t172 - t24 + t26 + V_base(2)) ^ 2;
t185 = t28 * t147;
t186 = t8 * t171;
t190 = (-t150 * t181 + t174 * t163 + t185 - t186 - t40 - t42 + V_base(3)) ^ 2;
t195 = Icges(4,4) * t6;
t197 = Icges(4,5) * t3;
t198 = Icges(4,1) * t154 * t6 + t195 * t152 + t197;
t203 = Icges(4,6) * t3;
t204 = Icges(4,2) * t152 * t6 + t195 * t154 + t203;
t206 = Icges(4,5) * t6;
t208 = Icges(4,6) * t6;
t211 = Icges(4,3) * t3 + t208 * t152 + t206 * t154;
t217 = Icges(4,4) * t3;
t219 = -Icges(4,1) * t154 * t3 - t217 * t152 + t206;
t224 = -Icges(4,2) * t152 * t3 - t217 * t154 + t208;
t229 = Icges(4,3) * t6 - t203 * t152 - t197 * t154;
t235 = Icges(4,1) * t152 - Icges(4,4) * t154;
t239 = Icges(4,4) * t152 - Icges(4,2) * t154;
t243 = Icges(4,5) * t152 - Icges(4,6) * t154;
t282 = Icges(3,4) * t6;
t284 = Icges(3,5) * t3;
t285 = Icges(3,1) * t11 * t6 + t282 * t9 + t284;
t290 = Icges(3,6) * t3;
t291 = Icges(3,2) * t6 * t9 + t282 * t11 + t290;
t293 = Icges(3,5) * t6;
t295 = Icges(3,6) * t6;
t298 = Icges(3,3) * t3 + t293 * t11 + t295 * t9;
t304 = Icges(3,4) * t3;
t306 = -Icges(3,1) * t11 * t3 - t304 * t9 + t293;
t311 = -Icges(3,2) * t3 * t9 - t304 * t11 + t295;
t316 = Icges(3,3) * t6 - t284 * t11 - t290 * t9;
t322 = Icges(3,1) * t9 - Icges(3,4) * t11;
t326 = Icges(3,4) * t9 - Icges(3,2) * t11;
t330 = Icges(3,5) * t9 - Icges(3,6) * t11;
t336 = m(3) * (t23 + t38 + t46) + V_base(4) * ((t3 * t55 + t6 * t51) * V_base(4) + (t3 * t63 + t6 * t60) * V_base(5) + t69 * t2) + V_base(5) * ((-t3 * t51 + t6 * t55) * V_base(4) + (-t3 * t60 + t6 * t63) * V_base(5) + t83 * t2) + m(2) * (t93 + t100 + t104) + t2 * (Icges(2,3) * t2 + t69 * V_base(4) + t83 * V_base(5)) + V_base(4) * (Icges(1,1) * V_base(4) + Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6)) + V_base(5) * (Icges(1,4) * V_base(4) + Icges(1,2) * V_base(5) + Icges(1,6) * V_base(6)) + V_base(6) * (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6)) + m(1) * (t130 + t134 + t138) + m(4) * (t166 + t184 + t190) + t174 * ((t176 * t198 + t178 * t204 + t3 * t211) * t174 + (t176 * t219 + t178 * t224 + t3 * t229) * t150 + (t176 * t235 + t178 * t239 + t3 * t243) * t2) + t150 * ((-t158 * t198 - t160 * t204 + t6 * t211) * t174 + (-t158 * t219 - t160 * t224 + t6 * t229) * t150 + (-t158 * t235 - t160 * t239 + t6 * t243) * t2) + t2 * ((t152 * t198 - t154 * t204) * t174 + (t152 * t219 - t154 * t224) * t150 + (t152 * t235 - t154 * t239) * t2) + t8 * ((-t15 * t285 - t17 * t291 + t6 * t298) * t28 + (-t15 * t306 - t17 * t311 + t6 * t316) * t8 + (-t15 * t322 - t17 * t326 + t6 * t330) * t2);
t369 = qJD(4) * t6 * t152;
t370 = sin(qJ(4));
t372 = cos(qJ(4));
t374 = t176 * t370 + t3 * t372;
t375 = qJD(5) * t374;
t376 = V_base(4) + t27 + t173 + t369 + t375;
t379 = t176 * t372 - t3 * t370;
t380 = cos(qJ(5));
t382 = sin(qJ(5));
t384 = t178 * t382 + t379 * t380;
t388 = t178 * t380 - t379 * t382;
t391 = Icges(6,1) * t384 + Icges(6,4) * t388 + Icges(6,5) * t374;
t396 = Icges(6,4) * t384 + Icges(6,2) * t388 + Icges(6,6) * t374;
t401 = Icges(6,5) * t384 + Icges(6,6) * t388 + Icges(6,3) * t374;
t407 = -t158 * t372 - t6 * t370;
t410 = -t160 * t382 + t407 * t380;
t414 = -t160 * t380 - t407 * t382;
t418 = -t158 * t370 + t6 * t372;
t420 = Icges(6,1) * t410 + Icges(6,4) * t414 + Icges(6,5) * t418;
t425 = Icges(6,4) * t410 + Icges(6,2) * t414 + Icges(6,6) * t418;
t430 = Icges(6,5) * t410 + Icges(6,6) * t414 + Icges(6,3) * t418;
t434 = qJD(4) * t3 * t152;
t435 = qJD(5) * t418;
t436 = V_base(5) + t7 + t149 - t434 + t435;
t438 = t152 * t372;
t441 = -t154 * t382 + t438 * t380;
t445 = -t154 * t380 - t438 * t382;
t449 = Icges(6,5) * t152 * t370 + Icges(6,1) * t441 + Icges(6,4) * t445;
t455 = Icges(6,6) * t152 * t370 + Icges(6,4) * t441 + Icges(6,2) * t445;
t461 = Icges(6,3) * t152 * t370 + Icges(6,5) * t441 + Icges(6,6) * t445;
t464 = qJD(4) * t154;
t466 = qJD(5) * t152 * t370;
t467 = V_base(6) - qJD(1) - t464 + t466;
t473 = pkin(4) * t152 - pkin(5) * t154;
t474 = t150 * t473;
t477 = -pkin(4) * t158 - pkin(5) * t160;
t478 = t2 * t477;
t481 = t152 * t370;
t483 = rSges(6,1) * t441 + rSges(6,2) * t445 + rSges(6,3) * t481;
t488 = rSges(6,1) * t410 + rSges(6,2) * t414 + rSges(6,3) * t418;
t491 = (t436 * t483 - t467 * t488 + t1 + t142 - t148 + t474 - t478 + t5 + V_base(1)) ^ 2;
t492 = t174 * t473;
t495 = pkin(4) * t176 + pkin(5) * t178;
t496 = t2 * t495;
t501 = rSges(6,1) * t384 + rSges(6,2) * t388 + rSges(6,3) * t374;
t504 = (-t376 * t483 + t467 * t501 - t168 + t172 - t24 + t26 - t492 + t496 + V_base(2)) ^ 2;
t505 = t174 * t477;
t506 = t150 * t495;
t510 = (t376 * t488 - t436 * t501 + t185 - t186 - t40 - t42 + t505 - t506 + V_base(3)) ^ 2;
t530 = V_base(4) + t27 + t173 + t369;
t535 = Icges(5,5) * t152 * t6 + Icges(5,1) * t379 - Icges(5,4) * t374;
t541 = Icges(5,6) * t152 * t6 + Icges(5,4) * t379 - Icges(5,2) * t374;
t547 = Icges(5,3) * t152 * t6 + Icges(5,5) * t379 - Icges(5,6) * t374;
t555 = -Icges(5,5) * t152 * t3 + Icges(5,1) * t407 - Icges(5,4) * t418;
t561 = -Icges(5,6) * t152 * t3 + Icges(5,4) * t407 - Icges(5,2) * t418;
t567 = -Icges(5,3) * t152 * t3 + Icges(5,5) * t407 - Icges(5,6) * t418;
t570 = V_base(5) + t7 + t149 - t434;
t574 = Icges(5,4) * t152;
t577 = Icges(5,1) * t152 * t372 - Icges(5,5) * t154 - t574 * t370;
t583 = -Icges(5,2) * t152 * t370 - Icges(5,6) * t154 + t574 * t372;
t590 = Icges(5,5) * t152 * t372 - Icges(5,6) * t152 * t370 - Icges(5,3) * t154;
t593 = V_base(6) - qJD(1) - t464;
t634 = rSges(5,1) * t438 - rSges(5,2) * t481 - rSges(5,3) * t154;
t639 = rSges(5,1) * t407 - rSges(5,2) * t418 - rSges(5,3) * t160;
t642 = (t570 * t634 - t593 * t639 + t1 + t142 - t148 + t474 - t478 + t5 + V_base(1)) ^ 2;
t647 = rSges(5,1) * t379 - rSges(5,2) * t374 + rSges(5,3) * t178;
t650 = (-t530 * t634 + t593 * t647 - t168 + t172 - t24 + t26 - t492 + t496 + V_base(2)) ^ 2;
t654 = (t530 * t639 - t570 * t647 + t185 - t186 - t40 - t42 + t505 - t506 + V_base(3)) ^ 2;
t658 = -pkin(7) * qJD(5) + qJD(6);
t660 = t658 * t414 + t149 - t434 + t435 + t7 + V_base(5);
t662 = pkin(7) * qJ(5) - qJ(6);
t663 = cos(t662);
t665 = sin(t662);
t667 = -t410 * t663 - t418 * t665;
t670 = -t374 * t665 - t384 * t663;
t674 = t374 * t663 - t384 * t665;
t677 = Icges(7,1) * t670 + Icges(7,4) * t674 + Icges(7,5) * t388;
t681 = -t410 * t665 + t418 * t663;
t685 = Icges(7,4) * t670 + Icges(7,2) * t674 + Icges(7,6) * t388;
t690 = Icges(7,5) * t670 + Icges(7,6) * t674 + Icges(7,3) * t388;
t694 = t658 * t388 + t173 + t27 + t369 + t375 + V_base(4);
t699 = Icges(7,1) * t667 + Icges(7,4) * t681 + Icges(7,5) * t414;
t704 = Icges(7,4) * t667 + Icges(7,2) * t681 + Icges(7,6) * t414;
t709 = Icges(7,5) * t667 + Icges(7,6) * t681 + Icges(7,3) * t414;
t715 = -t441 * t663 - t481 * t665;
t719 = -t441 * t665 + t481 * t663;
t722 = Icges(7,1) * t715 + Icges(7,4) * t719 + Icges(7,5) * t445;
t727 = Icges(7,4) * t715 + Icges(7,2) * t719 + Icges(7,6) * t445;
t732 = Icges(7,5) * t715 + Icges(7,6) * t719 + Icges(7,3) * t445;
t736 = t658 * t445 - qJD(1) - t464 + t466 + V_base(6);
t781 = rSges(7,1) * t715 + rSges(7,2) * t719 + rSges(7,3) * t445;
t786 = rSges(7,1) * t667 + rSges(7,2) * t681 + rSges(7,3) * t414;
t788 = -pkin(6) * t414 * t467 + pkin(6) * t436 * t445 + t660 * t781 - t736 * t786 + t1 + t142 - t148 + t474 - t478 + t5 + V_base(1);
t789 = t788 ^ 2;
t798 = rSges(7,1) * t670 + rSges(7,2) * t674 + rSges(7,3) * t388;
t800 = -pkin(6) * t376 * t445 + pkin(6) * t388 * t467 - t694 * t781 + t736 * t798 - t168 + t172 - t24 + t26 - t492 + t496 + V_base(2);
t801 = t800 ^ 2;
t808 = pkin(6) * t376 * t414 - pkin(6) * t388 * t436 - t660 * t798 + t694 * t786 + t185 - t186 - t40 - t42 + t505 - t506 + V_base(3);
t809 = t808 ^ 2;
t829 = t2 * ((-t11 * t291 + t9 * t285) * t28 + (-t11 * t311 + t9 * t306) * t8 + (-t11 * t326 + t9 * t322) * t2) + t28 * ((t30 * t285 + t32 * t291 + t3 * t298) * t28 + (t3 * t316 + t30 * t306 + t32 * t311) * t8 + (t3 * t330 + t30 * t322 + t32 * t326) * t2) + t376 * ((t374 * t401 + t384 * t391 + t388 * t396) * t376 + (t374 * t430 + t384 * t420 + t388 * t425) * t436 + (t374 * t461 + t384 * t449 + t388 * t455) * t467) + m(6) * (t491 + t504 + t510) + t436 * ((t410 * t391 + t414 * t396 + t418 * t401) * t376 + (t410 * t420 + t414 * t425 + t418 * t430) * t436 + (t410 * t449 + t414 * t455 + t418 * t461) * t467) + t530 * ((t178 * t547 - t374 * t541 + t379 * t535) * t530 + (t178 * t567 - t374 * t561 + t379 * t555) * t570 + (t178 * t590 - t374 * t583 + t379 * t577) * t593) + t570 * ((-t160 * t547 + t407 * t535 - t418 * t541) * t530 + (-t160 * t567 + t407 * t555 - t418 * t561) * t570 + (-t160 * t590 + t407 * t577 - t418 * t583) * t593) + t593 * ((-t154 * t547 + t438 * t535 - t481 * t541) * t530 + (-t154 * t567 + t438 * t555 - t481 * t561) * t570 + (-t154 * t590 + t438 * t577 - t481 * t583) * t593) + m(5) * (t642 + t650 + t654) + t660 * ((t414 * t690 + t667 * t677 + t681 * t685) * t694 + (t414 * t709 + t667 * t699 + t681 * t704) * t660 + (t414 * t732 + t667 * t722 + t681 * t727) * t736) + t694 * ((t388 * t690 + t670 * t677 + t674 * t685) * t694 + (t388 * t709 + t670 * t699 + t674 * t704) * t660 + (t388 * t732 + t670 * t722 + t674 * t727) * t736) + t467 * ((t441 * t391 + t445 * t396 + t481 * t401) * t376 + (t441 * t420 + t445 * t425 + t481 * t430) * t436 + (t441 * t449 + t445 * t455 + t481 * t461) * t467) + m(7) * (t789 + t801 + t809) + t736 * ((t445 * t690 + t715 * t677 + t719 * t685) * t694 + (t445 * t709 + t715 * t699 + t719 * t704) * t660 + (t445 * t732 + t715 * t722 + t719 * t727) * t736);
t830 = t336 / 0.2e1 + t829 / 0.2e1;
T = t830;
