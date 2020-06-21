% Calculate joint inertia matrix for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350DE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:27
% EndTime: 2020-06-19 21:37:28
% DurationCPUTime: 0.84s
% Computational Cost: add. (2985->490), mult. (5646->660), div. (0->0), fcn. (5457->10), ass. (0->301)
unknown=NaN(21,1);
t2 = pkin(7) * qJ(5) - qJ(6);
t3 = sin(t2);
t4 = cos(qJ(5));
t5 = cos(qJ(4));
t6 = t4 * t5;
t7 = cos(qJ(3));
t8 = cos(qJ(2));
t10 = sin(qJ(3));
t11 = sin(qJ(2));
t13 = t10 * t11 - t7 * t8;
t15 = sin(qJ(5));
t18 = t10 * t8 + t7 * t11;
t20 = t6 * t13 + t15 * t18;
t22 = cos(t2);
t23 = sin(qJ(4));
t26 = t22 * t23 * t13 - t3 * t20;
t30 = -t3 * t23 * t13 - t22 * t20;
t33 = t15 * t5;
t36 = -t33 * t13 + t4 * t18;
t38 = Ifges(7,4) * t30 + Ifges(7,2) * t26 + Ifges(7,6) * t36;
t40 = t4 * t23;
t41 = t11 * pkin(3);
t44 = -t18 * pkin(4) + t13 * pkin(5) - pkin(2) - t41;
t46 = t23 * t13;
t48 = -t46 * pkin(6) - t40 * t44;
t50 = t5 * t44;
t52 = t20 * pkin(6) + t50;
t54 = -t22 * t48 - t3 * t52;
t57 = -t36 * mrSges(7,2) + t26 * mrSges(7,3);
t60 = t54 ^ 2;
t63 = t22 * t52 - t3 * t48;
t64 = t63 ^ 2;
t65 = t15 ^ 2;
t66 = t23 ^ 2;
t68 = t44 ^ 2;
t69 = t65 * t66 * t68;
t75 = t36 * mrSges(7,1) - t30 * mrSges(7,3);
t82 = Ifges(6,5) * t23 * t13 + Ifges(6,1) * t20 + Ifges(6,4) * t36;
t88 = Ifges(6,6) * t23 * t13 + Ifges(6,4) * t20 + Ifges(6,2) * t36;
t90 = t4 ^ 2;
t93 = t5 ^ 2;
t94 = t93 * t68;
t99 = Ifges(5,5) * t5 * t13;
t101 = Ifges(5,6) * t23 * t13;
t102 = Ifges(5,3) * t18;
t117 = -pkin(2) - t41;
t123 = pkin(2) ^ 2;
t125 = Ifges(2,3) + t26 * t38 + 0.2e1 * t54 * t57 + m(7) * (t60 + t64 + t69) + 0.2e1 * t63 * t75 + t20 * t82 + t36 * t88 + m(6) * (t90 * t66 * t68 + t69 + t94) + t18 * (t99 - t101 + t102) + m(5) * (t66 * t68 + t94) + t13 * (Ifges(4,1) * t13 + Ifges(4,4) * t18) + t18 * (Ifges(4,4) * t13 + Ifges(4,2) * t18) + 0.2e1 * t117 * (-t18 * mrSges(4,1) + t13 * mrSges(4,2)) + m(3) * t123;
t126 = t117 ^ 2;
t141 = Ifges(7,5) * t30;
t142 = Ifges(7,6) * t26;
t143 = Ifges(7,3) * t36;
t144 = t141 + t142 + t143;
t149 = Ifges(7,1) * t30 + Ifges(7,4) * t26 + Ifges(7,5) * t36;
t151 = t15 * t23;
t154 = -t26 * mrSges(7,1) + t30 * mrSges(7,2);
t160 = t46 * mrSges(6,1) - t20 * mrSges(6,3);
t166 = -t46 * mrSges(6,2) + t36 * mrSges(6,3);
t172 = -t36 * mrSges(6,1) + t20 * mrSges(6,2);
t175 = t23 * t44;
t178 = -t18 * mrSges(5,2) - t46 * mrSges(5,3);
t181 = t5 * t13;
t184 = t18 * mrSges(5,1) - t181 * mrSges(5,3);
t189 = Ifges(5,4) * t23;
t192 = Ifges(5,1) * t5 * t13 + Ifges(5,5) * t18 - t189 * t13;
t194 = Ifges(5,4) * t5;
t199 = -Ifges(5,2) * t23 * t13 + Ifges(5,6) * t18 + t194 * t13;
t201 = Ifges(6,5) * t20;
t202 = Ifges(6,6) * t36;
t204 = Ifges(6,3) * t23 * t13;
t205 = t201 + t202 + t204;
t207 = m(4) * t126 - t8 * (-Ifges(3,1) * t8 + Ifges(3,4) * t11) + t11 * (-Ifges(3,4) * t8 + Ifges(3,2) * t11) - 0.2e1 * pkin(2) * (-t11 * mrSges(3,1) - t8 * mrSges(3,2)) + t36 * t144 + t30 * t149 + 0.2e1 * t151 * t44 * t154 + 0.2e1 * t151 * t44 * t160 - 0.2e1 * t40 * t44 * t166 + 0.2e1 * t50 * t172 - 0.2e1 * t175 * t178 - 0.2e1 * t50 * t184 + t181 * t192 - t46 * t199 + t46 * t205;
t209 = t10 * pkin(3);
t210 = t209 - pkin(5);
t211 = t6 * t210;
t212 = t7 * pkin(3);
t213 = t212 + pkin(4);
t214 = t15 * t213;
t215 = t5 * pkin(6);
t216 = t211 + t214 - t215;
t218 = t23 * t210;
t219 = t40 * pkin(6);
t220 = t218 - t219;
t222 = -t22 * t216 - t3 * t220;
t226 = -t3 * t216 + t22 * t220;
t230 = -t33 * t210 + t4 * t213;
t231 = t230 * t15;
t232 = t231 * t175;
t241 = t22 * t4 * t23 - t3 * t5;
t244 = t151 * mrSges(7,1) - t241 * mrSges(7,3);
t245 = t63 * t244;
t247 = t241 * t149 / 0.2e1;
t251 = t3 * t4 * t23 + t22 * t5;
t253 = t251 * t38 / 0.2e1;
t257 = Ifges(6,4) * t15;
t260 = -Ifges(6,1) * t4 * t23 + Ifges(6,5) * t5 + t257 * t23;
t262 = t20 * t260 / 0.2e1;
t263 = Ifges(6,4) * t4;
t268 = Ifges(6,2) * t15 * t23 + Ifges(6,6) * t5 - t263 * t23;
t270 = t36 * t268 / 0.2e1;
t272 = t5 * t205 / 0.2e1;
t273 = t211 + t214;
t274 = t273 * t4;
t284 = t46 * mrSges(5,1) + t181 * mrSges(5,2);
t287 = t23 * t192 / 0.2e1;
t289 = t5 * t199 / 0.2e1;
t290 = Ifges(5,5) * t23;
t291 = Ifges(5,6) * t5;
t294 = t18 * (-t290 - t291) / 0.2e1;
t297 = Ifges(4,5) * t13;
t298 = m(7) * (t222 * t54 + t226 * t63 + t232) + t222 * t57 + t230 * t154 + t245 + t247 + t253 + t226 * t75 + t262 + t270 + t272 + m(6) * (-t274 * t175 + t218 * t50 + t232) + t273 * t166 + t230 * t160 + t213 * t284 - t287 - t289 + t294 - Ifges(3,5) * t8 + Ifges(3,6) * t11 + t297;
t299 = Ifges(4,6) * t18;
t300 = Ifges(7,5) * t241;
t301 = Ifges(7,6) * t251;
t303 = Ifges(7,3) * t15 * t23;
t304 = t300 + t301 + t303;
t306 = t36 * t304 / 0.2e1;
t311 = Ifges(7,5) * t15 * t23 + Ifges(7,1) * t241 + Ifges(7,4) * t251;
t313 = t30 * t311 / 0.2e1;
t318 = Ifges(7,6) * t15 * t23 + Ifges(7,4) * t241 + Ifges(7,2) * t251;
t320 = t26 * t318 / 0.2e1;
t323 = -t151 * mrSges(7,2) + t251 * mrSges(7,3);
t324 = t54 * t323;
t327 = -t251 * mrSges(7,1) + t241 * mrSges(7,2);
t329 = t151 * t44 * t327;
t332 = t5 * mrSges(6,1) + t40 * mrSges(6,3);
t334 = t151 * t44 * t332;
t337 = -t5 * mrSges(6,2) + t151 * mrSges(6,3);
t339 = t40 * t44 * t337;
t345 = t151 * t144 / 0.2e1;
t348 = t40 * t82 / 0.2e1;
t350 = t151 * t88 / 0.2e1;
t351 = t5 * t210;
t355 = Ifges(6,5) * t4 * t23;
t357 = Ifges(6,6) * t15 * t23;
t358 = Ifges(6,3) * t5;
t359 = -t355 + t357 + t358;
t361 = t46 * t359 / 0.2e1;
t364 = -t151 * mrSges(6,1) - t40 * mrSges(6,2);
t365 = t50 * t364;
t367 = -Ifges(5,2) * t5 - t189;
t369 = t46 * t367 / 0.2e1;
t371 = -Ifges(5,1) * t23 - t194;
t373 = t181 * t371 / 0.2e1;
t374 = -t212 * t13 * mrSges(4,3) + t209 * t18 * mrSges(4,3) + t218 * t172 + t351 * t178 - t218 * t184 + t299 + t306 + t313 + t320 + t324 + t329 + t334 - t339 + t345 - t348 + t350 + t361 + t365 - t369 + t373;
t376 = t222 * t323;
t378 = t251 * t318;
t379 = t241 * t311;
t380 = t226 * t244;
t382 = t230 * t327;
t384 = t222 ^ 2;
t385 = t226 ^ 2;
t386 = t230 ^ 2;
t390 = t5 * t359;
t391 = t273 * t337;
t393 = t273 ^ 2;
t394 = t210 ^ 2;
t395 = t66 * t394;
t399 = t230 * t332;
t401 = Ifges(3,3) + Ifges(4,3) + 0.2e1 * t376 + t378 + t379 + 0.2e1 * t380 + 0.2e1 * t382 + m(7) * (t384 + t385 + t386) + t390 + 0.2e1 * t391 + m(6) * (t393 + t386 + t395) + 0.2e1 * t399;
t404 = t5 * mrSges(5,1) - t23 * mrSges(5,2);
t405 = t213 * t404;
t407 = t23 * t371;
t408 = t5 * t367;
t410 = t213 ^ 2;
t414 = t10 ^ 2;
t415 = pkin(3) ^ 2;
t417 = t7 ^ 2;
t423 = t93 * t210 * mrSges(5,3);
t425 = t151 * t304;
t426 = t218 * t364;
t428 = t40 * t260;
t429 = t151 * t268;
t430 = t209 * mrSges(4,2);
t432 = t212 * mrSges(4,1);
t435 = t66 * t210 * mrSges(5,3);
t437 = 0.2e1 * t405 - t407 - t408 + m(5) * (t93 * t394 + t395 + t410) + m(4) * (t414 * t415 + t417 * t415) - 0.2e1 * t423 + t425 + 0.2e1 * t426 - t428 + t429 - 0.2e1 * t430 + 0.2e1 * t432 - 0.2e1 * t435;
t439 = t6 * pkin(5);
t440 = t15 * pkin(4);
t441 = -t439 + t440 - t215;
t443 = t23 * pkin(5);
t444 = -t443 - t219;
t446 = -t22 * t441 - t3 * t444;
t450 = t22 * t444 - t3 * t441;
t454 = t4 * pkin(4) + t33 * pkin(5);
t455 = t454 * t15;
t456 = t455 * t175;
t463 = -t439 + t440;
t464 = t463 * t4;
t473 = t245 + t247 + t253 + t262 + t270 + t272 - t287 - t289 + t294 + m(7) * (t446 * t54 + t450 * t63 + t456) + t446 * t57 + t454 * t154 + t450 * t75 + m(6) * (-t464 * t175 - t443 * t50 + t456) + t463 * t166 + t454 * t160 + pkin(4) * t284 + t297;
t475 = t5 * pkin(5);
t478 = -t443 * t172 - t475 * t178 + t443 * t184 + t299 + t306 + t313 + t320 + t324 + t329 + t334 - t339 + t345 - t348 + t350 + t361 + t365 - t369 + t373;
t480 = t463 * t337;
t481 = t454 * t332;
t482 = pkin(4) * t404;
t483 = t93 * pkin(5);
t485 = t66 * pkin(5);
t486 = t485 * t210;
t491 = Ifges(4,3) + t480 + t481 + t482 + m(5) * (pkin(4) * t213 - t483 * t210 - t486) + t376 + t378 + t379 + t380 + t382 + t390 + t391 + t399 + t405 - t407 - t408;
t492 = t454 * t327;
t495 = t454 * t230;
t499 = t446 * t323;
t500 = t450 * t244;
t505 = t443 * t364;
t506 = t483 * mrSges(5,3);
t507 = t485 * mrSges(5,3);
t508 = t492 + m(7) * (t446 * t222 + t450 * t226 + t495) + t499 + t500 + m(6) * (t463 * t273 - t486 + t495) - t423 - t505 + t506 + t507 + t425 + t426 - t428 + t429 - t430 + t432 - t435;
t519 = t446 ^ 2;
t520 = t450 ^ 2;
t521 = t454 ^ 2;
t525 = t463 ^ 2;
t526 = pkin(5) ^ 2;
t527 = t66 * t526;
t532 = pkin(4) ^ 2;
t537 = -t408 + 0.2e1 * t492 + 0.2e1 * t499 + 0.2e1 * t500 + 0.2e1 * t480 + 0.2e1 * t481 + 0.2e1 * t482 + m(7) * (t519 + t520 + t521) + m(6) * (t525 + t521 + t527) + m(5) * (t93 * t526 + t527 + t532) + 0.2e1 * t507;
t540 = Ifges(7,5) * t22 * t15;
t542 = Ifges(7,6) * t3 * t15;
t543 = Ifges(7,3) * t4;
t544 = -t540 - t542 + t543;
t549 = Ifges(7,4) * t3;
t552 = -Ifges(7,1) * t22 * t15 + Ifges(7,5) * t4 - t549 * t15;
t555 = Ifges(7,4) * t22;
t560 = -Ifges(7,2) * t3 * t15 + Ifges(7,6) * t4 - t555 * t15;
t563 = t3 * t15;
t566 = -t4 * mrSges(7,2) - t563 * mrSges(7,3);
t570 = t22 * t15;
t580 = t4 * mrSges(7,1) + t570 * mrSges(7,3);
t583 = Ifges(6,1) * t15 + t263;
t587 = Ifges(6,2) * t4 + t257;
t594 = t36 * t544 / 0.2e1 + t30 * t552 / 0.2e1 + t26 * t560 / 0.2e1 + t54 * t566 + m(7) * (-t563 * pkin(6) * t54 + t570 * pkin(6) * t63) + t4 * t144 / 0.2e1 + t63 * t580 + t20 * t583 / 0.2e1 + t36 * t587 / 0.2e1 + t15 * t82 / 0.2e1 + t4 * t88 / 0.2e1 + t102;
t598 = t44 * mrSges(6,3);
t602 = t563 * mrSges(7,1) - t570 * mrSges(7,2);
t613 = Ifges(6,5) * t15;
t614 = Ifges(6,6) * t4;
t615 = t613 + t614;
t620 = -t4 * mrSges(6,1) + t15 * mrSges(6,2);
t624 = t570 * pkin(6) * t75 - t90 * t23 * t598 + t151 * t44 * t602 - t65 * t23 * t598 - t563 * pkin(6) * t57 - t570 * t149 / 0.2e1 - t563 * t38 / 0.2e1 + t46 * t615 / 0.2e1 + t50 * t620 + t175 * mrSges(5,2) - t50 * mrSges(5,1) + t99 - t101;
t627 = t570 * t311 / 0.2e1;
t629 = t563 * t318 / 0.2e1;
t633 = t151 * t544 / 0.2e1;
t636 = t40 * t583 / 0.2e1;
t638 = t151 * t587 / 0.2e1;
t641 = -t218 * mrSges(5,1) - t351 * mrSges(5,2) - t231 * mrSges(6,3) + t274 * mrSges(6,3) + t218 * t620 - t290 - t291 - t627 - t629 + t633 - t636 + t638;
t644 = t251 * t560 / 0.2e1;
t646 = t241 * t552 / 0.2e1;
t657 = t4 * t304 / 0.2e1;
t659 = t5 * t615 / 0.2e1;
t661 = t15 * t260 / 0.2e1;
t663 = t4 * t268 / 0.2e1;
t665 = t563 * pkin(6) * t323;
t667 = t570 * pkin(6) * t244;
t668 = t222 * t566 + t644 + t646 + t226 * t580 + t230 * t602 + m(7) * (-t563 * pkin(6) * t222 + t570 * pkin(6) * t226) + t657 + t659 + t661 + t663 - t665 + t667;
t675 = t443 * mrSges(5,1) + t475 * mrSges(5,2) - t455 * mrSges(6,3) + t464 * mrSges(6,3) - t443 * t620 - t290 - t291 - t627 - t629 + t633 - t636 + t638;
t686 = t644 + t646 + t657 + t659 + t661 + t663 + t446 * t566 + t450 * t580 + m(7) * (-t563 * pkin(6) * t446 + t570 * pkin(6) * t450) + t454 * t602 - t665 + t667;
t697 = t3 ^ 2;
t699 = pkin(6) ^ 2;
t701 = t22 ^ 2;
t712 = -t22 * mrSges(7,1) - t3 * mrSges(7,2);
t719 = Ifges(7,5) * t3;
t720 = Ifges(7,6) * t22;
t721 = Ifges(7,3) * pkin(7);
t722 = -t719 + t720 - t721;
t727 = -Ifges(7,1) * t3 - Ifges(7,5) * pkin(7) + t555;
t732 = Ifges(7,2) * t22 - Ifges(7,6) * pkin(7) - t549;
t737 = pkin(7) * mrSges(7,2) + t22 * mrSges(7,3);
t739 = t22 * pkin(6);
t741 = t3 * pkin(6);
t751 = -pkin(7) * mrSges(7,1) + t3 * mrSges(7,3);
t758 = t151 * t44 * t712 + t151 * t44 * mrSges(6,1) + t40 * t44 * mrSges(6,2) + t36 * t722 / 0.2e1 + t30 * t727 / 0.2e1 + t26 * t732 / 0.2e1 + t54 * t737 + m(7) * (t739 * t54 + t741 * t63) + t739 * t57 - pkin(7) * t144 / 0.2e1 + t63 * t751 - t3 * t149 / 0.2e1 + t22 * t38 / 0.2e1 + t741 * t75 + t201 + t202 + t204;
t761 = t251 * t732 / 0.2e1;
t763 = t241 * t727 / 0.2e1;
t772 = t3 * t311 / 0.2e1;
t774 = t22 * t318 / 0.2e1;
t776 = pkin(7) * t304 / 0.2e1;
t777 = t739 * t323;
t778 = t741 * t244;
t782 = t151 * t722 / 0.2e1;
t783 = t222 * t737 + t761 + t763 + t226 * t751 + t230 * t712 + m(7) * (t739 * t222 + t741 * t226) - t772 + t774 - t776 + t777 + t778 + t358 - t273 * mrSges(6,2) - t355 + t357 + t230 * mrSges(6,1) + t782;
t794 = t761 + t763 - t772 + t774 - t776 + t777 + t778 + t358 - t355 + t357 + t446 * t737 + t450 * t751 + m(7) * (t739 * t446 + t741 * t450) + t454 * t712 - t463 * mrSges(6,2) + t454 * mrSges(6,1) + t782;
t813 = t570 * pkin(6) * t751 - t563 * pkin(6) * t737 - pkin(7) * t544 / 0.2e1 + t22 * t560 / 0.2e1 - t3 * t552 / 0.2e1 + t4 * t722 / 0.2e1 + t739 * t566 + t741 * t580 + t613 + t614 - t563 * t732 / 0.2e1 - t570 * t727 / 0.2e1;
unknown(1,1) = t125 + t207;
unknown(2,1) = t298 + t374;
unknown(3,1) = t401 + t437;
unknown(4,1) = t473 + t478;
unknown(5,1) = t491 + t508;
unknown(6,1) = t425 - t428 + t429 - 0.2e1 * t505 + 0.2e1 * t506 + Ifges(4,3) + t378 + t379 + t390 - t407 + t537;
unknown(7,1) = t594 + t624;
unknown(8,1) = t641 + t668;
unknown(9,1) = t675 + t686;
unknown(10,1) = 0.2e1 * t570 * pkin(6) * t580 - 0.2e1 * t563 * pkin(6) * t566 + Ifges(5,3) + t4 * t544 - t563 * t560 - t570 * t552 + m(7) * (t697 * t65 * t699 + t701 * t65 * t699) + t4 * t587 + t15 * t583;
unknown(11,1) = t758;
unknown(12,1) = t783;
unknown(13,1) = t794;
unknown(14,1) = t813;
unknown(15,1) = Ifges(6,3) + 0.2e1 * t739 * t737 + 0.2e1 * t741 * t751 - t3 * t727 + m(7) * (t697 * t699 + t701 * t699) + t22 * t732 - pkin(7) * t722;
unknown(16,1) = t63 * mrSges(7,1) - t54 * mrSges(7,2) + t141 + t142 + t143;
unknown(17,1) = t226 * mrSges(7,1) - t222 * mrSges(7,2) + t300 + t301 + t303;
unknown(18,1) = t450 * mrSges(7,1) - t446 * mrSges(7,2) + t300 + t301 + t303;
unknown(19,1) = t570 * pkin(6) * mrSges(7,1) + t563 * pkin(6) * mrSges(7,2) - t540 - t542 + t543;
unknown(20,1) = t741 * mrSges(7,1) - t739 * mrSges(7,2) - t719 + t720 - t721;
unknown(21,1) = Ifges(7,3);
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [unknown(1), unknown(2), unknown(4), unknown(7), unknown(11), unknown(16); unknown(2), unknown(3), unknown(5), unknown(8), unknown(12), unknown(17); unknown(4), unknown(5), unknown(6), unknown(9), unknown(13), unknown(18); unknown(7), unknown(8), unknown(9), unknown(10), unknown(14), unknown(19); unknown(11), unknown(12), unknown(13), unknown(14), unknown(15), unknown(20); unknown(16), unknown(17), unknown(18), unknown(19), unknown(20), unknown(21);];
Mq = res;
