% Calculate kinetic energy for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
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

function T = CloosQRC350DE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:24
% EndTime: 2020-06-19 21:37:25
% DurationCPUTime: 0.52s
% Computational Cost: add. (2517->384), mult. (3481->556), div. (0->0), fcn. (3738->12), ass. (0->177)
t1 = sin(qJ(1));
t3 = qJD(1) * t1 * pkin(2);
t4 = cos(qJ(1));
t5 = qJD(2) * t4;
t6 = cos(qJ(2));
t8 = sin(qJ(2));
t10 = t6 * rSges(3,1) - t8 * rSges(3,2);
t12 = t1 * t8;
t14 = t1 * t6;
t17 = -t12 * rSges(3,1) - t14 * rSges(3,2) + t4 * rSges(3,3);
t20 = (qJD(1) * t17 + t5 * t10 - t3) ^ 2;
t22 = qJD(1) * t4 * pkin(2);
t23 = qJD(2) * t1;
t25 = t4 * t8;
t27 = t4 * t6;
t30 = t25 * rSges(3,1) + t27 * rSges(3,2) + t1 * rSges(3,3);
t33 = (-qJD(1) * t30 - t23 * t10 - t22) ^ 2;
t37 = (t23 * t17 - t5 * t30) ^ 2;
t40 = qJD(1) ^ 2;
t44 = (-t1 * rSges(2,1) + t4 * rSges(2,2)) ^ 2;
t49 = (t4 * rSges(2,1) + t1 * rSges(2,2)) ^ 2;
t53 = qJD(3) * t1;
t54 = t23 + t53;
t55 = qJ(2) + qJ(3);
t56 = sin(t55);
t57 = t4 * t56;
t60 = Icges(4,4) * t4;
t61 = cos(t55);
t63 = Icges(4,5) * t1;
t64 = Icges(4,1) * t4 * t56 + t60 * t61 + t63;
t66 = t4 * t61;
t70 = Icges(4,6) * t1;
t71 = Icges(4,2) * t4 * t61 + t60 * t56 + t70;
t73 = Icges(4,5) * t4;
t75 = Icges(4,6) * t4;
t78 = Icges(4,3) * t1 + t73 * t56 + t75 * t61;
t84 = Icges(4,4) * t1;
t86 = -Icges(4,1) * t1 * t56 - t84 * t61 + t73;
t91 = -Icges(4,2) * t1 * t61 - t84 * t56 + t75;
t96 = Icges(4,3) * t4 - t63 * t56 - t70 * t61;
t99 = qJD(3) * t4;
t100 = t5 + t99;
t104 = Icges(4,1) * t61 - Icges(4,4) * t56;
t108 = Icges(4,4) * t61 - Icges(4,2) * t56;
t112 = Icges(4,5) * t61 - Icges(4,6) * t56;
t118 = t1 * t56;
t120 = t1 * t61;
t153 = Icges(3,4) * t4;
t155 = Icges(3,5) * t1;
t156 = Icges(3,1) * t4 * t8 + t153 * t6 + t155;
t161 = Icges(3,6) * t1;
t162 = Icges(3,2) * t4 * t6 + t153 * t8 + t161;
t169 = Icges(3,4) * t1;
t171 = Icges(3,5) * t4;
t172 = -Icges(3,1) * t1 * t8 - t169 * t6 + t171;
t177 = Icges(3,6) * t4;
t178 = -Icges(3,2) * t1 * t6 - t169 * t8 + t177;
t185 = Icges(3,1) * t6 - Icges(3,4) * t8;
t189 = Icges(3,4) * t6 - Icges(3,2) * t8;
t195 = qJD(4) * t56;
t197 = sin(qJ(4));
t198 = qJD(5) * t61 * t197;
t200 = -pkin(7) * qJD(5) + qJD(6);
t201 = cos(qJ(4));
t202 = t61 * t201;
t203 = sin(qJ(5));
t205 = cos(qJ(5));
t207 = -t202 * t203 - t56 * t205;
t209 = t200 * t207 - qJD(1) - t195 + t198;
t212 = t202 * t205 - t56 * t203;
t214 = pkin(7) * qJ(5) - qJ(6);
t215 = cos(t214);
t217 = t61 * t197;
t218 = sin(t214);
t220 = -t212 * t215 - t217 * t218;
t223 = -t1 * t197 + t57 * t201;
t226 = t66 * t203 + t223 * t205;
t230 = t1 * t201 + t57 * t197;
t232 = -t226 * t215 - t230 * t218;
t236 = t230 * t215 - t226 * t218;
t240 = -t223 * t203 + t66 * t205;
t242 = Icges(7,1) * t232 + Icges(7,4) * t236 + Icges(7,5) * t240;
t246 = -t212 * t218 + t217 * t215;
t250 = Icges(7,4) * t232 + Icges(7,2) * t236 + Icges(7,6) * t240;
t255 = Icges(7,5) * t232 + Icges(7,6) * t236 + Icges(7,3) * t240;
t259 = qJD(4) * t4 * t61;
t260 = qJD(5) * t230;
t262 = t200 * t240 + t23 + t259 + t260 + t53;
t266 = -t118 * t201 - t4 * t197;
t269 = -t120 * t203 + t266 * t205;
t273 = -t118 * t197 + t4 * t201;
t275 = -t269 * t215 - t273 * t218;
t279 = t273 * t215 - t269 * t218;
t283 = -t120 * t205 - t266 * t203;
t285 = Icges(7,1) * t275 + Icges(7,4) * t279 + Icges(7,5) * t283;
t290 = Icges(7,4) * t275 + Icges(7,2) * t279 + Icges(7,6) * t283;
t295 = Icges(7,5) * t275 + Icges(7,6) * t279 + Icges(7,3) * t283;
t299 = qJD(4) * t1 * t61;
t300 = qJD(5) * t273;
t302 = t200 * t283 - t299 + t300 + t5 + t99;
t307 = Icges(7,1) * t220 + Icges(7,4) * t246 + Icges(7,5) * t207;
t312 = Icges(7,4) * t220 + Icges(7,2) * t246 + Icges(7,6) * t207;
t317 = Icges(7,5) * t220 + Icges(7,6) * t246 + Icges(7,3) * t207;
t357 = -qJD(1) - t195 + t198;
t361 = Icges(6,1) * t226 + Icges(6,4) * t240 + Icges(6,5) * t230;
t366 = Icges(6,4) * t226 + Icges(6,2) * t240 + Icges(6,6) * t230;
t371 = Icges(6,5) * t226 + Icges(6,6) * t240 + Icges(6,3) * t230;
t374 = t23 + t53 + t259 + t260;
t379 = Icges(6,1) * t269 + Icges(6,4) * t283 + Icges(6,5) * t273;
t384 = Icges(6,4) * t269 + Icges(6,2) * t283 + Icges(6,6) * t273;
t389 = Icges(6,5) * t269 + Icges(6,6) * t283 + Icges(6,3) * t273;
t392 = t5 + t99 - t299 + t300;
t398 = Icges(6,5) * t61 * t197 + Icges(6,1) * t212 + Icges(6,4) * t207;
t404 = Icges(6,6) * t61 * t197 + Icges(6,4) * t212 + Icges(6,2) * t207;
t410 = Icges(6,3) * t61 * t197 + Icges(6,5) * t212 + Icges(6,6) * t207;
t433 = m(3) * (t20 + t33 + t37) + m(2) * (t40 * t44 + t40 * t49) + t54 * ((t1 * t78 + t57 * t64 + t66 * t71) * t54 + (t1 * t96 + t57 * t86 + t66 * t91) * t100 - (t1 * t112 + t57 * t104 + t66 * t108) * qJD(1)) + t100 * ((-t118 * t64 - t120 * t71 + t4 * t78) * t54 + (-t118 * t86 - t120 * t91 + t4 * t96) * t100 - (-t118 * t104 - t120 * t108 + t4 * t112) * qJD(1)) - qJD(1) * ((-t56 * t71 + t61 * t64) * t54 + (-t56 * t91 + t61 * t86) * t100 - (t61 * t104 - t56 * t108) * qJD(1)) - qJD(1) * ((t6 * t156 - t8 * t162) * qJD(2) * t1 + (t6 * t172 - t8 * t178) * qJD(2) * t4 - (t6 * t185 - t8 * t189) * qJD(1)) + t209 * ((t207 * t255 + t220 * t242 + t246 * t250) * t262 + (t207 * t295 + t220 * t285 + t246 * t290) * t302 + (t207 * t317 + t220 * t307 + t246 * t312) * t209) + t302 * ((t275 * t242 + t279 * t250 + t283 * t255) * t262 + (t275 * t285 + t279 * t290 + t283 * t295) * t302 + (t275 * t307 + t279 * t312 + t283 * t317) * t209) + t262 * ((t232 * t242 + t236 * t250 + t240 * t255) * t262 + (t232 * t285 + t236 * t290 + t240 * t295) * t302 + (t232 * t307 + t236 * t312 + t240 * t317) * t209) + t357 * ((t207 * t366 + t212 * t361 + t217 * t371) * t374 + (t207 * t384 + t212 * t379 + t217 * t389) * t392 + (t207 * t404 + t212 * t398 + t217 * t410) * t357) + t374 * ((t226 * t361 + t230 * t371 + t240 * t366) * t374 + (t226 * t379 + t230 * t389 + t240 * t384) * t392 + (t226 * t398 + t230 * t410 + t240 * t404) * t357);
t451 = t23 + t53 + t259;
t456 = Icges(5,5) * t4 * t61 + Icges(5,1) * t223 - Icges(5,4) * t230;
t462 = Icges(5,6) * t4 * t61 + Icges(5,4) * t223 - Icges(5,2) * t230;
t468 = Icges(5,3) * t4 * t61 + Icges(5,5) * t223 - Icges(5,6) * t230;
t476 = -Icges(5,5) * t1 * t61 + Icges(5,1) * t266 - Icges(5,4) * t273;
t482 = -Icges(5,6) * t1 * t61 + Icges(5,4) * t266 - Icges(5,2) * t273;
t488 = -Icges(5,3) * t1 * t61 + Icges(5,5) * t266 - Icges(5,6) * t273;
t491 = t5 + t99 - t299;
t495 = Icges(5,4) * t61;
t498 = Icges(5,1) * t61 * t201 - Icges(5,5) * t56 - t495 * t197;
t504 = -Icges(5,2) * t61 * t197 - Icges(5,6) * t56 + t495 * t201;
t511 = Icges(5,5) * t61 * t201 - Icges(5,6) * t61 * t197 - Icges(5,3) * t56;
t514 = -qJD(1) - t195;
t552 = t6 * pkin(3);
t553 = t5 * t552;
t556 = t8 * pkin(3) + pkin(2);
t558 = t1 * pkin(2) - t1 * t556;
t559 = qJD(1) * t558;
t562 = t61 * pkin(4) - t56 * pkin(5);
t563 = t100 * t562;
t566 = -t118 * pkin(4) - t120 * pkin(5);
t567 = qJD(1) * t566;
t575 = t220 * rSges(7,1) + t246 * rSges(7,2) + t207 * rSges(7,3);
t580 = t275 * rSges(7,1) + t279 * rSges(7,2) + t283 * rSges(7,3);
t583 = (t392 * t207 * pkin(6) - t357 * t283 * pkin(6) - t209 * t580 + t302 * t575 - t3 + t553 + t559 + t563 + t567) ^ 2;
t584 = t23 * t552;
t587 = -t4 * pkin(2) + t4 * t556;
t588 = qJD(1) * t587;
t589 = t54 * t562;
t592 = t57 * pkin(4) + t66 * pkin(5);
t593 = qJD(1) * t592;
t602 = t232 * rSges(7,1) + t236 * rSges(7,2) + t240 * rSges(7,3);
t605 = (-t374 * t207 * pkin(6) + t357 * t240 * pkin(6) + t209 * t602 - t262 * t575 - t22 - t584 - t588 - t589 - t593) ^ 2;
t606 = t23 * t558;
t607 = t5 * t587;
t608 = t54 * t566;
t609 = t100 * t592;
t617 = (-t392 * t240 * pkin(6) + t374 * t283 * pkin(6) + t262 * t580 - t302 * t602 + t606 - t607 + t608 - t609) ^ 2;
t623 = t212 * rSges(6,1) + t207 * rSges(6,2) + t217 * rSges(6,3);
t628 = t269 * rSges(6,1) + t283 * rSges(6,2) + t273 * rSges(6,3);
t631 = (-t357 * t628 + t392 * t623 - t3 + t553 + t559 + t563 + t567) ^ 2;
t636 = t226 * rSges(6,1) + t240 * rSges(6,2) + t230 * rSges(6,3);
t639 = (t357 * t636 - t374 * t623 - t22 - t584 - t588 - t589 - t593) ^ 2;
t643 = (t374 * t628 - t392 * t636 + t606 - t607 + t608 - t609) ^ 2;
t649 = t202 * rSges(5,1) - t217 * rSges(5,2) - t56 * rSges(5,3);
t654 = t266 * rSges(5,1) - t273 * rSges(5,2) - t120 * rSges(5,3);
t657 = (t491 * t649 - t514 * t654 - t3 + t553 + t559 + t563 + t567) ^ 2;
t662 = t223 * rSges(5,1) - t230 * rSges(5,2) + t66 * rSges(5,3);
t665 = (-t451 * t649 + t514 * t662 - t22 - t584 - t588 - t589 - t593) ^ 2;
t669 = (t451 * t654 - t491 * t662 + t606 - t607 + t608 - t609) ^ 2;
t674 = t61 * rSges(4,1) - t56 * rSges(4,2);
t679 = -t118 * rSges(4,1) - t120 * rSges(4,2) + t4 * rSges(4,3);
t682 = (qJD(1) * t679 + t100 * t674 - t3 + t553 + t559) ^ 2;
t687 = t57 * rSges(4,1) + t66 * rSges(4,2) + t1 * rSges(4,3);
t690 = (-qJD(1) * t687 - t54 * t674 - t22 - t584 - t588) ^ 2;
t694 = (-t100 * t687 + t54 * t679 + t606 - t607) ^ 2;
t703 = Icges(3,3) * t1 + t8 * t171 + t177 * t6;
t713 = Icges(3,3) * t4 - t155 * t8 - t161 * t6;
t722 = Icges(3,5) * t6 - Icges(3,6) * t8;
t747 = t392 * ((t269 * t361 + t273 * t371 + t283 * t366) * t374 + (t269 * t379 + t273 * t389 + t283 * t384) * t392 + (t269 * t398 + t273 * t410 + t283 * t404) * t357) + t451 * ((t223 * t456 - t230 * t462 + t66 * t468) * t451 + (t223 * t476 - t230 * t482 + t66 * t488) * t491 + (t223 * t498 - t230 * t504 + t66 * t511) * t514) + t491 * ((-t120 * t468 + t266 * t456 - t273 * t462) * t451 + (-t120 * t488 + t266 * t476 - t273 * t482) * t491 + (-t120 * t511 + t266 * t498 - t273 * t504) * t514) + t514 * ((t202 * t456 - t217 * t462 - t56 * t468) * t451 + (t202 * t476 - t217 * t482 - t56 * t488) * t491 + (t202 * t498 - t217 * t504 - t56 * t511) * t514) + m(7) * (t583 + t605 + t617) + m(6) * (t631 + t639 + t643) + m(5) * (t657 + t665 + t669) + m(4) * (t682 + t690 + t694) + t40 * Icges(2,3) + t23 * ((t1 * t703 + t25 * t156 + t27 * t162) * qJD(2) * t1 + (t1 * t713 + t25 * t172 + t27 * t178) * qJD(2) * t4 - (t1 * t722 + t25 * t185 + t27 * t189) * qJD(1)) + t5 * ((-t12 * t156 - t14 * t162 + t4 * t703) * qJD(2) * t1 + (-t12 * t172 - t14 * t178 + t4 * t713) * qJD(2) * t4 - (-t12 * t185 - t14 * t189 + t4 * t722) * qJD(1));
t748 = t433 / 0.2e1 + t747 / 0.2e1;
T = t748;
