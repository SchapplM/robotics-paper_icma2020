% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% CloosQRC350DE
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = CloosQRC350DE_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'CloosQRC350DE_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
JaD_transl=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t28 = cos(qJ(1));
	t27 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t28 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; (r_i_i_C(1) * t27 - r_i_i_C(2) * t28) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (17->13), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->11)
	t17 = sin(qJ(1));
	t25 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t24 = qJD(1) * t19;
	t23 = qJD(2) * t17;
	t22 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t21 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18 + pkin(2);
	t20 = (-r_i_i_C(1) * t18 + r_i_i_C(2) * t16) * qJD(2);
	t1 = [t17 * t20 + (-r_i_i_C(3) * t17 - t19 * t21) * qJD(1), (t16 * t25 - t18 * t22) * r_i_i_C(2) + (-t16 * t22 - t18 * t25) * r_i_i_C(1), 0, 0, 0, 0; t19 * t20 + (-r_i_i_C(3) * t19 + t17 * t21) * qJD(1), (t16 * t24 + t18 * t23) * r_i_i_C(2) + (t16 * t23 - t18 * t24) * r_i_i_C(1), 0, 0, 0, 0; 0, t20, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (77->25), mult. (110->41), div. (0->0), fcn. (71->6), ass. (0->26)
	t44 = cos(qJ(2));
	t61 = pkin(3) * t44;
	t41 = qJ(2) + qJ(3);
	t39 = cos(t41);
	t60 = r_i_i_C(1) * t39;
	t38 = sin(t41);
	t59 = r_i_i_C(2) * t38;
	t40 = qJD(2) + qJD(3);
	t58 = t38 * t40;
	t57 = t39 * t40;
	t43 = sin(qJ(1));
	t56 = qJD(1) * t43;
	t45 = cos(qJ(1));
	t55 = qJD(1) * t45;
	t42 = sin(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(2) * t57;
	t52 = qJD(1) * t59;
	t51 = qJD(2) * t61;
	t50 = -r_i_i_C(1) * t57 + r_i_i_C(2) * t58;
	t49 = -r_i_i_C(1) * t38 - r_i_i_C(2) * t39;
	t48 = t42 * pkin(3) + pkin(2) - t49;
	t47 = t45 * t52 - t55 * t60 + (t58 * r_i_i_C(1) + t53) * t43;
	t46 = -t51 + (t59 - t60) * t40;
	t34 = t43 * t52;
	t1 = [t46 * t43 + (-r_i_i_C(3) * t43 - t48 * t45) * qJD(1), t34 + (-t60 - t61) * t56 + (-pkin(3) * t54 + t49 * t40) * t45, -t45 * t53 + t34 + (-t39 * t56 - t45 * t58) * r_i_i_C(1), 0, 0, 0; t46 * t45 + (-r_i_i_C(3) * t45 + t48 * t43) * qJD(1), (t43 * t54 - t44 * t55) * pkin(3) + t47, t47, 0, 0, 0; 0, t50 - t51, t50, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:01
	% EndTime: 2020-06-23 21:15:01
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (266->56), mult. (390->93), div. (0->0), fcn. (295->8), ass. (0->53)
	t300 = pkin(5) + r_i_i_C(3);
	t255 = qJ(2) + qJ(3);
	t252 = sin(t255);
	t306 = t300 * t252;
	t253 = cos(t255);
	t256 = sin(qJ(4));
	t290 = qJD(4) * t256;
	t254 = qJD(2) + qJD(3);
	t259 = cos(qJ(4));
	t293 = t254 * t259;
	t305 = t252 * t293 + t253 * t290;
	t274 = qJD(4) * t252 + qJD(1);
	t303 = t259 * t274;
	t302 = t274 * t256;
	t257 = sin(qJ(2));
	t299 = pkin(4) * t252;
	t301 = t257 * pkin(3) + t300 * t253 + pkin(2) + t299;
	t297 = pkin(3) * qJD(2);
	t296 = t252 * t254;
	t261 = cos(qJ(1));
	t295 = t253 * t261;
	t258 = sin(qJ(1));
	t294 = t254 * t258;
	t292 = qJD(1) * t258;
	t291 = qJD(1) * t261;
	t289 = qJD(4) * t259;
	t288 = r_i_i_C(2) * t253 * t256;
	t287 = t257 * t297;
	t260 = cos(qJ(2));
	t286 = t260 * t297;
	t285 = t256 * t296;
	t283 = t253 * t294;
	t282 = -r_i_i_C(1) * t259 - pkin(4);
	t278 = t253 * t289;
	t276 = r_i_i_C(2) * t285;
	t275 = qJD(1) * t288;
	t273 = qJD(1) * t252 + qJD(4);
	t272 = t258 * t275 + t261 * t276 + t292 * t306;
	t271 = t282 * t254;
	t270 = qJD(1) * t282;
	t269 = t273 * t261;
	t268 = t252 * t271;
	t267 = t261 * t275 + t294 * t299 + t300 * (t252 * t291 + t283) + (r_i_i_C(1) * t305 + r_i_i_C(2) * t278) * t258;
	t266 = qJD(1) * (-pkin(3) * t260 + t282 * t253);
	t265 = -t254 * t295 + t273 * t258;
	t264 = -t300 * t254 + (-r_i_i_C(1) * t256 - r_i_i_C(2) * t259) * qJD(4);
	t263 = t253 * t271 + t254 * t288 + t300 * t296 + (r_i_i_C(1) * t290 + r_i_i_C(2) * t289) * t252;
	t262 = -t286 + (-pkin(4) * t253 + t306) * t254;
	t232 = t259 * t269 + (t253 * t293 - t302) * t258;
	t231 = t258 * t303 + (t269 + t283) * t256;
	t230 = t265 * t259 + t261 * t302;
	t229 = t265 * t256 - t261 * t303;
	t1 = [-t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t262 * t258 - t301 * t291, t258 * t266 + (t264 * t253 + t268 - t287) * t261 + t272, t261 * t268 + (t258 * t270 + t264 * t261) * t253 + t272, t229 * r_i_i_C(1) + t230 * r_i_i_C(2), 0, 0; t230 * r_i_i_C(1) - t229 * r_i_i_C(2) + t262 * t261 + t301 * t292, (-t276 + t287) * t258 + t261 * t266 + t267, -t258 * t276 + t270 * t295 + t267, t231 * r_i_i_C(1) + t232 * r_i_i_C(2), 0, 0; 0, t263 - t286, t263, t305 * r_i_i_C(2) + (-t278 + t285) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:03
	% EndTime: 2020-06-23 21:15:04
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (562->111), mult. (860->189), div. (0->0), fcn. (777->10), ass. (0->81)
	t357 = sin(qJ(5));
	t361 = cos(qJ(5));
	t355 = qJD(2) + qJD(3);
	t362 = cos(qJ(4));
	t386 = qJD(5) * t362 + t355;
	t358 = sin(qJ(4));
	t408 = qJD(4) * t358;
	t367 = t386 * t357 + t361 * t408;
	t435 = t357 * t408 - t386 * t361;
	t356 = qJ(2) + qJ(3);
	t353 = sin(t356);
	t354 = cos(t356);
	t360 = sin(qJ(1));
	t410 = qJD(1) * t360;
	t364 = cos(qJ(1));
	t417 = t355 * t364;
	t434 = t353 * t417 + t354 * t410;
	t418 = t355 * t362;
	t387 = qJD(5) + t418;
	t428 = t387 * t361;
	t433 = t353 * t428 + t367 * t354;
	t431 = t435 * t354;
	t427 = pkin(5) * t354;
	t426 = r_i_i_C(3) * t358;
	t425 = t361 * r_i_i_C(2);
	t424 = pkin(3) * qJD(2);
	t423 = t353 * t355;
	t422 = t353 * t357;
	t421 = t354 * t355;
	t420 = t354 * t362;
	t419 = t355 * t360;
	t416 = t360 * t358;
	t415 = t360 * t362;
	t414 = t364 * t358;
	t413 = t364 * t362;
	t412 = t434 * t357;
	t411 = qJD(1) * t353;
	t409 = qJD(1) * t364;
	t407 = qJD(4) * t362;
	t340 = t353 * t413 - t416;
	t406 = qJD(5) * t340;
	t405 = qJD(5) * t354;
	t404 = qJD(5) * t357;
	t403 = qJD(5) * t361;
	t359 = sin(qJ(2));
	t401 = t359 * t424;
	t400 = pkin(5) * t411;
	t363 = cos(qJ(2));
	t399 = t363 * t424;
	t398 = t353 * t419;
	t396 = t354 * t419;
	t395 = t353 * t416;
	t394 = t354 * t417;
	t393 = -pkin(4) - t426;
	t392 = pkin(5) + t425;
	t388 = t359 * pkin(3) + pkin(4) * t353 + pkin(2);
	t385 = r_i_i_C(3) * t354 * t407;
	t384 = qJD(4) * t353 + qJD(1);
	t383 = qJD(4) + t411;
	t370 = qJD(1) * (t353 * t361 + t357 * t420);
	t371 = qJD(1) * (-t361 * t420 + t422);
	t376 = t387 * t357;
	t382 = (t360 * t370 + (t353 * t376 + t431) * t364) * r_i_i_C(2) + (t360 * t371 - t433 * t364) * r_i_i_C(1) + t364 * t385 + t360 * t400;
	t381 = t393 * t353;
	t380 = t361 * r_i_i_C(1) - t357 * r_i_i_C(2);
	t372 = t353 * t414 + t415;
	t373 = t353 * t415 + t414;
	t334 = t373 * qJD(1) + t372 * qJD(4) - t362 * t394;
	t379 = -t364 * t405 + t334;
	t336 = -t383 * t413 + (-t354 * t418 + t384 * t358) * t360;
	t378 = t360 * t405 - t336;
	t374 = (t364 * t370 + (-t387 * t422 - t431) * t360) * r_i_i_C(2) + (t433 * t360 + t364 * t371) * r_i_i_C(1) + t355 * r_i_i_C(3) * t395 + pkin(4) * t398 + pkin(5) * t396 + t364 * t400;
	t369 = qJD(1) * (-pkin(3) * t363 + t393 * t354);
	t366 = qJD(5) * t373 - t354 * t409 + t398;
	t365 = -pkin(4) * t421 + (-t353 * t435 + t354 * t376) * r_i_i_C(2) + (t367 * t353 - t354 * t428) * r_i_i_C(1) + pkin(5) * t423 + (-t353 * t407 - t358 * t421) * r_i_i_C(3);
	t337 = t395 - t413;
	t335 = t384 * t415 + (t383 * t364 + t396) * t358;
	t333 = -t384 * t413 + (t383 * t360 - t394) * t358;
	t326 = t366 * t357 - t378 * t361;
	t325 = t378 * t357 + t366 * t361;
	t1 = [t326 * r_i_i_C(1) + t325 * r_i_i_C(2) - t335 * r_i_i_C(3) + (-t399 + (-t354 * pkin(4) + t353 * pkin(5)) * t355) * t360 + (-t388 - t427) * t409, t360 * t369 + (-t401 + (t381 - t427) * t355) * t364 + t382, t381 * t417 + (-pkin(5) * t417 + t393 * t410) * t354 + t382, -t334 * r_i_i_C(3) + (-t333 * t357 + t372 * t403) * r_i_i_C(2) + (t333 * t361 + t372 * t404) * r_i_i_C(1), t412 * r_i_i_C(2) + (t379 * r_i_i_C(1) + r_i_i_C(2) * t406) * t357 + ((-t406 - t434) * r_i_i_C(1) + t379 * r_i_i_C(2)) * t361, 0; (t334 * t361 + t340 * t404 + t412) * r_i_i_C(1) + (-t334 * t357 + t340 * t403) * r_i_i_C(2) + t333 * r_i_i_C(3) + (t392 * t354 + t388) * t410 + (-t399 + t392 * t423 + (-t355 * pkin(4) - t380 * qJD(5)) * t354) * t364, (-t385 + t401) * t360 + t364 * t369 + t374, (-pkin(4) * t409 + (-t358 * t409 - t360 * t407) * r_i_i_C(3)) * t354 + t374, t336 * r_i_i_C(3) + (-t335 * t357 - t337 * t403) * r_i_i_C(2) + (t335 * t361 - t337 * t404) * r_i_i_C(1), t325 * r_i_i_C(1) - t326 * r_i_i_C(2), 0; 0, t365 - t399, t365, (-r_i_i_C(3) * t362 + t380 * t358) * t423 + ((t357 * r_i_i_C(1) + t425) * t358 * qJD(5) + (-t380 * t362 - t426) * qJD(4)) * t354, (r_i_i_C(1) * t376 + r_i_i_C(2) * t428) * t353 + (t435 * r_i_i_C(1) + t367 * r_i_i_C(2)) * t354, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:07
	% EndTime: 2020-06-23 21:15:09
	% DurationCPUTime: 2.85s
	% Computational Cost: add. (1852->168), mult. (2762->275), div. (0->0), fcn. (2475->12), ass. (0->124)
	t550 = qJ(2) + qJ(3);
	t547 = sin(t550);
	t552 = sin(qJ(4));
	t556 = cos(qJ(4));
	t623 = qJD(1) * t547;
	t590 = qJD(4) + t623;
	t558 = cos(qJ(1));
	t618 = qJD(4) * t556;
	t601 = t558 * t618;
	t548 = cos(t550);
	t549 = qJD(2) + qJD(3);
	t631 = t549 * t558;
	t607 = t548 * t631;
	t621 = qJD(1) * t558;
	t554 = sin(qJ(1));
	t629 = t554 * t552;
	t507 = -t547 * t601 - t552 * t607 - t556 * t621 + t590 * t629;
	t625 = t558 * t556;
	t526 = t547 * t625 - t629;
	t545 = pkin(7) * qJD(5) - qJD(6);
	t551 = sin(qJ(5));
	t555 = cos(qJ(5));
	t635 = t548 * t558;
	t597 = (t526 * t555 + t551 * t635) * t545 + t507;
	t673 = t597 * r_i_i_C(1);
	t672 = t597 * r_i_i_C(2);
	t650 = r_i_i_C(3) + pkin(6);
	t632 = t549 * t556;
	t593 = qJD(5) + t632;
	t638 = t547 * t555;
	t592 = qJD(5) * t556 + t549;
	t619 = qJD(4) * t552;
	t566 = t592 * t551 + t555 * t619;
	t665 = t566 * t548;
	t504 = t593 * t638 + t665;
	t642 = t545 * t552;
	t612 = t548 * t642;
	t586 = -t504 + t612;
	t546 = pkin(7) * qJ(5) - qJ(6);
	t540 = sin(t546);
	t541 = cos(t546);
	t588 = t541 * r_i_i_C(1) + t540 * r_i_i_C(2);
	t559 = (t588 * t551 - t650 * t555) * qJD(5);
	t641 = t545 * t555;
	t659 = -t540 * r_i_i_C(1) + t541 * r_i_i_C(2);
	t671 = t659 * t641 - t559;
	t608 = t547 * t631;
	t622 = qJD(1) * t554;
	t670 = t548 * t622 + t608;
	t633 = t549 * t554;
	t660 = t547 * t633 - t548 * t621;
	t613 = t650 * t551;
	t669 = t588 * t555 + t613;
	t626 = t558 * t552;
	t628 = t554 * t556;
	t525 = t547 * t626 + t628;
	t576 = t547 * t628 + t626;
	t508 = qJD(1) * t576 + qJD(4) * t525 - t556 * t607;
	t617 = qJD(5) * t548;
	t584 = t558 * t617 - t508;
	t616 = qJD(5) * t551;
	t624 = t670 * t551;
	t490 = -t526 * t616 + t555 * t584 - t624;
	t645 = t525 * t545;
	t599 = t490 + t645;
	t668 = (-t599 * r_i_i_C(2) + t673) * t541 + (t599 * r_i_i_C(1) + t672) * t540;
	t634 = t549 * t552;
	t627 = t555 * t556;
	t639 = t547 * t551;
	t522 = t548 * t627 - t639;
	t666 = t522 * t545;
	t564 = -t547 * t634 + t548 * t618 - t666;
	t667 = (r_i_i_C(1) * t586 - r_i_i_C(2) * t564) * t540 - (r_i_i_C(1) * t564 + r_i_i_C(2) * t586) * t541;
	t567 = -t551 * t619 + t592 * t555;
	t664 = t567 * t548;
	t591 = qJD(4) * t547 + qJD(1);
	t510 = -t590 * t625 + (-t548 * t632 + t591 * t552) * t554;
	t661 = -t554 * t617 + t510;
	t658 = t548 * (t552 * t622 - t601) + t552 * t608 + t558 * t666;
	t609 = t547 * t629;
	t657 = t548 * (t552 * t621 + t554 * t618) - t549 * t609 - t554 * t666;
	t652 = t551 * t584 + t555 * (qJD(5) * t526 + t670);
	t557 = cos(qJ(2));
	t649 = pkin(3) * t557;
	t648 = pkin(4) * t548;
	t644 = t540 * t545;
	t643 = t541 * t545;
	t640 = t547 * t549;
	t637 = t548 * t551;
	t636 = t548 * t554;
	t630 = t552 * t555;
	t553 = sin(qJ(2));
	t620 = qJD(2) * t553;
	t615 = pkin(5) * t623;
	t614 = qJD(2) * t649;
	t610 = t548 * t633;
	t513 = -t526 * t551 + t555 * t635;
	t600 = qJD(5) * t513 - t508 * t555 - t624 + t645;
	t492 = t660 * t551 + t661 * t555 + t576 * t616;
	t523 = t609 - t625;
	t598 = t523 * t545 - t492;
	t509 = t591 * t628 + (t558 * t590 + t610) * t552;
	t595 = (-t551 * t636 - t555 * t576) * t545 + t509;
	t589 = -pkin(4) * t547 - pkin(5) * t548;
	t582 = t593 * t555;
	t585 = t548 * t582 + (-t566 + t642) * t547;
	t581 = t593 * t551;
	t503 = t547 * t581 - t664;
	t577 = t556 * t637 + t638;
	t569 = qJD(1) * t577;
	t570 = qJD(1) * t522;
	t580 = t554 * t570 + (t547 * t582 - t612 + t665) * t558;
	t583 = (t580 * t540 - t658 * t541) * r_i_i_C(2) + (t658 * t540 + t580 * t541) * r_i_i_C(1) + t554 * t615 + t650 * (t503 * t558 + t554 * t569);
	t579 = t554 * t586 + t558 * t570;
	t578 = t553 * pkin(3) + pkin(2) - t589;
	t565 = (-t547 * t627 - t637) * t545 + t547 * t618 + t548 * t634;
	t568 = -t549 * t648 + (t540 * t585 - t541 * t565) * r_i_i_C(2) + (t540 * t565 + t541 * t585) * r_i_i_C(1) + pkin(5) * t640 + t650 * (t547 * t567 + t548 * t581);
	t562 = -t614 + (pkin(5) * t547 - t648) * t549;
	t560 = (t579 * t540 - t657 * t541) * r_i_i_C(2) + (t657 * t540 + t579 * t541) * r_i_i_C(1) + pkin(5) * t610 + t558 * t615 + t650 * (t558 * t569 + (-t593 * t639 + t664) * t554) + t660 * pkin(4);
	t511 = t551 * t576 - t555 * t636;
	t491 = -t661 * t551 + (qJD(5) * t576 + t660) * t555;
	t480 = -t595 * t540 - t598 * t541;
	t479 = t598 * t540 - t595 * t541;
	t1 = [-t480 * r_i_i_C(1) + t479 * r_i_i_C(2) + t650 * t491 + t562 * t554 - t578 * t621, (-t648 - t649) * t622 + (-pkin(3) * t620 + t549 * t589) * t558 + t583, -pkin(4) * t670 - pkin(5) * t607 + t583, (t508 * t540 - t526 * t643) * r_i_i_C(1) + (-t508 * t541 - t526 * t644) * r_i_i_C(2) - t669 * t507 + t671 * t525, (t513 * t644 + t541 * t652) * r_i_i_C(1) + (-t513 * t643 + t540 * t652) * r_i_i_C(2) - t650 * t490 + t668 * pkin(7), -t668; (t600 * r_i_i_C(1) + t672) * t541 + (t600 * r_i_i_C(2) - t673) * t540 + t650 * t652 + t562 * t558 + t578 * t622, (t554 * t620 - t557 * t621) * pkin(3) + t560, t560, (-t510 * t540 + t576 * t643) * r_i_i_C(1) + (t510 * t541 + t576 * t644) * r_i_i_C(2) - t669 * t509 - t671 * t523, (-t491 * t541 + t511 * t644) * r_i_i_C(1) + (-t491 * t540 - t511 * t643) * r_i_i_C(2) - t650 * t492 + ((t595 * r_i_i_C(1) + t598 * r_i_i_C(2)) * t541 + (-t598 * r_i_i_C(1) + t595 * r_i_i_C(2)) * t540) * pkin(7), t479 * r_i_i_C(1) + t480 * r_i_i_C(2); 0, t568 - t614, t568, ((t540 * t556 - t541 * t630) * r_i_i_C(1) + (-t540 * t630 - t541 * t556) * r_i_i_C(2) - t552 * t613) * t640 + ((qJD(4) * t669 - t588 * t545) * t556 + (-t559 + t659 * (-qJD(4) + t641)) * t552) * t548, (-t503 * t541 - t577 * t644) * r_i_i_C(1) + (-t503 * t540 + t577 * t643) * r_i_i_C(2) + t650 * t504 + t667 * pkin(7), -t667;];
	JaD_transl = t1;
end