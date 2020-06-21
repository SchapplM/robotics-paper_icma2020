% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% CloosQRC350OL
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = CloosQRC350OL_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'CloosQRC350OL_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350OL_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
JaD_transl=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t18 - r_i_i_C(2) * t16;
	t21 = -r_i_i_C(1) * t16 - r_i_i_C(2) * t18 - pkin(2);
	t20 = t22 * qJD(2);
	t1 = [-t17 * t20 + (r_i_i_C(3) * t17 + t19 * t21) * qJD(1), (t16 * t26 - t18 * t23) * r_i_i_C(2) + (-t16 * t23 - t18 * t26) * r_i_i_C(1), 0, 0, 0, 0; t22 * t23 + (-r_i_i_C(3) * t19 + t17 * t21) * qJD(1), (-t16 * t25 - t18 * t24) * r_i_i_C(2) + (-t16 * t24 + t18 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t20, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (77->24), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
	t43 = cos(qJ(2));
	t39 = qJD(2) + qJD(3);
	t40 = qJ(2) + qJ(3);
	t37 = sin(t40);
	t57 = r_i_i_C(2) * t37;
	t38 = cos(t40);
	t58 = r_i_i_C(1) * t38;
	t49 = (t57 - t58) * t39;
	t53 = pkin(3) * qJD(2);
	t63 = -t43 * t53 + t49;
	t60 = pkin(3) * t43;
	t59 = r_i_i_C(1) * t37;
	t56 = r_i_i_C(2) * t38;
	t42 = sin(qJ(1));
	t55 = t39 * t42;
	t44 = cos(qJ(1));
	t54 = t39 * t44;
	t52 = qJD(1) * t42;
	t51 = qJD(1) * t44;
	t47 = -t56 - t59;
	t41 = sin(qJ(2));
	t46 = -t41 * pkin(3) - pkin(2) + t47;
	t45 = t47 * t39 - t41 * t53;
	t34 = t51 * t58;
	t33 = t52 * t57;
	t1 = [t63 * t42 + (r_i_i_C(3) * t42 + t46 * t44) * qJD(1), t33 + (-t58 - t60) * t52 + t45 * t44, -t54 * t56 + t33 + (-t37 * t54 - t38 * t52) * r_i_i_C(1), 0, 0, 0; -t63 * t44 + (-r_i_i_C(3) * t44 + t46 * t42) * qJD(1), t34 + (-t57 + t60) * t51 + t45 * t42, -t55 * t59 + t34 + (-t37 * t51 - t38 * t55) * r_i_i_C(2), 0, 0, 0; 0, t63, t49, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:17
	% EndTime: 2020-06-20 08:27:18
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (266->52), mult. (390->85), div. (0->0), fcn. (295->8), ass. (0->48)
	t251 = qJ(2) + qJ(3);
	t249 = cos(t251);
	t255 = cos(qJ(4));
	t294 = r_i_i_C(1) * t255 + pkin(4);
	t297 = t249 * t294;
	t248 = sin(t251);
	t292 = pkin(5) + r_i_i_C(3);
	t277 = t292 * t248;
	t252 = sin(qJ(4));
	t279 = r_i_i_C(2) * t249 * t252;
	t296 = t277 + t279;
	t250 = qJD(2) + qJD(3);
	t256 = cos(qJ(2));
	t288 = pkin(3) * qJD(2);
	t278 = t256 * t288;
	t295 = (-pkin(4) * t249 + t277) * t250 - t278;
	t291 = pkin(3) * t256;
	t287 = t248 * t250;
	t286 = t249 * t250;
	t285 = t250 * t255;
	t257 = cos(qJ(1));
	t284 = t255 * t257;
	t254 = sin(qJ(1));
	t283 = qJD(1) * t254;
	t282 = qJD(1) * t257;
	t281 = qJD(4) * t252;
	t280 = qJD(4) * t255;
	t276 = t252 * t287;
	t271 = r_i_i_C(2) * t276;
	t275 = t254 * t271 + t282 * t297;
	t270 = qJD(4) * t248 - qJD(1);
	t269 = qJD(1) * t248 - qJD(4);
	t268 = t257 * t271 + t296 * t283;
	t267 = t294 * t250;
	t266 = t270 * t252;
	t265 = t248 * t267;
	t253 = sin(qJ(2));
	t263 = qJD(1) * (-t253 * pkin(3) - pkin(4) * t248 - t292 * t249 - pkin(2));
	t262 = t269 * t254 - t257 * t286;
	t261 = -t292 * t250 + (-r_i_i_C(1) * t252 - r_i_i_C(2) * t255) * qJD(4);
	t260 = -t249 * t267 + t250 * t279 + t292 * t287 + (r_i_i_C(1) * t281 + r_i_i_C(2) * t280) * t248;
	t259 = t261 * t249 - t265;
	t258 = -t253 * t288 + t259;
	t234 = -t269 * t284 + (-t249 * t285 + t266) * t254;
	t233 = t270 * t255 * t254 + (t254 * t286 + t269 * t257) * t252;
	t232 = t262 * t255 + t257 * t266;
	t231 = t262 * t252 - t270 * t284;
	t1 = [t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t295 * t254 + t257 * t263, (-t291 - t297) * t283 + t258 * t257 + t268, -t257 * t265 + (t261 * t257 - t283 * t294) * t249 + t268, t231 * r_i_i_C(1) + t232 * r_i_i_C(2), 0, 0; -t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t254 * t263 - t295 * t257, (-t296 + t291) * t282 + t258 * t254 + t275, t259 * t254 - t282 * t296 + t275, -t233 * r_i_i_C(1) + t234 * r_i_i_C(2), 0, 0; 0, t260 - t278, t260, (t248 * t285 + t249 * t281) * r_i_i_C(2) + (-t249 * t280 + t276) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:19
	% EndTime: 2020-06-20 08:27:20
	% DurationCPUTime: 1.24s
	% Computational Cost: add. (562->111), mult. (860->186), div. (0->0), fcn. (777->10), ass. (0->80)
	t367 = sin(qJ(4));
	t426 = r_i_i_C(3) * t367;
	t437 = pkin(4) + t426;
	t366 = sin(qJ(5));
	t370 = cos(qJ(5));
	t364 = qJD(2) + qJD(3);
	t371 = cos(qJ(4));
	t395 = qJD(5) * t371 + t364;
	t413 = qJD(4) * t367;
	t433 = t395 * t366 + t370 * t413;
	t436 = t366 * t413 - t395 * t370;
	t365 = qJ(2) + qJ(3);
	t363 = cos(t365);
	t369 = sin(qJ(1));
	t362 = sin(t365);
	t392 = qJD(1) * t362 - qJD(4);
	t373 = cos(qJ(1));
	t420 = t364 * t373;
	t432 = -t363 * t420 + t392 * t369;
	t372 = cos(qJ(2));
	t430 = pkin(3) * t372;
	t429 = pkin(5) * t362;
	t428 = pkin(5) * t363;
	t427 = r_i_i_C(2) * t370;
	t425 = pkin(3) * qJD(2);
	t424 = t362 * t364;
	t423 = t363 * t364;
	t422 = t363 * t371;
	t421 = t364 * t369;
	t419 = t369 * t367;
	t418 = t369 * t371;
	t417 = t371 * t373;
	t416 = t373 * t367;
	t415 = qJD(1) * t369;
	t414 = qJD(1) * t373;
	t412 = qJD(4) * t371;
	t411 = qJD(5) * t366;
	t410 = qJD(5) * t369;
	t409 = qJD(5) * t370;
	t407 = pkin(5) * t414;
	t406 = t372 * t425;
	t405 = t362 * t421;
	t403 = t362 * t419;
	t402 = t363 * t421;
	t400 = t363 * t414;
	t368 = sin(qJ(2));
	t397 = -pkin(3) * t368 - pkin(4) * t362 - pkin(2);
	t396 = t364 * t371 + qJD(5);
	t394 = r_i_i_C(3) * t363 * t412;
	t393 = -qJD(4) * t362 + qJD(1);
	t383 = t396 * t366;
	t374 = t362 * t383 + t436 * t363;
	t384 = t396 * t370;
	t375 = -t362 * t384 - t433 * t363;
	t380 = t362 * t370 + t366 * t422;
	t381 = -t362 * t366 + t370 * t422;
	t391 = (t374 * t373 + t380 * t415) * r_i_i_C(2) + (t375 * t373 - t381 * t415) * r_i_i_C(1) + t373 * t394 + t415 * t429;
	t390 = t437 * t362;
	t389 = t437 * t369;
	t388 = r_i_i_C(1) * t370 - r_i_i_C(2) * t366;
	t344 = -t432 * t371 + t393 * t416;
	t387 = qJD(5) * t363 * t373 + t344;
	t350 = t362 * t417 + t419;
	t346 = t350 * qJD(1) - qJD(4) * t403 + t371 * t402 - t373 * t412;
	t386 = -t363 * t410 - t346;
	t385 = (t374 * t369 - t380 * t414) * r_i_i_C(2) + (t375 * t369 + t381 * t414) * r_i_i_C(1) + t369 * t394 + t437 * t400;
	t382 = t393 * t371;
	t348 = t362 * t418 - t416;
	t379 = qJD(5) * t348 - t400;
	t378 = -qJD(5) * t350 - t362 * t420 - t363 * t415;
	t377 = -t368 * t425 + (-t390 - t428) * t364;
	t376 = -pkin(4) * t423 + (-t362 * t436 + t363 * t383) * r_i_i_C(2) + (t433 * t362 - t363 * t384) * r_i_i_C(1) + pkin(5) * t424 + (-t362 * t412 - t367 * t423) * r_i_i_C(3);
	t351 = t366 * t405;
	t349 = -t362 * t416 + t418;
	t347 = -t403 - t417;
	t345 = t369 * t382 + (-t392 * t373 - t402) * t367;
	t343 = t432 * t367 + t373 * t382;
	t336 = t378 * t366 + t387 * t370;
	t335 = -t387 * t366 + t378 * t370;
	t1 = [(-t346 * t370 + t348 * t411 + t351) * r_i_i_C(1) + (t346 * t366 + t348 * t409) * r_i_i_C(2) + t345 * r_i_i_C(3) + (-t406 + (pkin(5) + t427) * t424) * t369 + t397 * t414 + ((-t366 * t414 - t369 * t409) * r_i_i_C(1) + (t366 * t410 - t370 * t414) * r_i_i_C(2) - pkin(4) * t421 - t407) * t363, (-t363 * t437 - t430) * t415 + t377 * t373 + t391, -t390 * t420 + (-pkin(5) * t420 - qJD(1) * t389) * t363 + t391, t344 * r_i_i_C(3) + (-t343 * t366 - t349 * t409) * r_i_i_C(2) + (t343 * t370 - t349 * t411) * r_i_i_C(1), r_i_i_C(1) * t335 - r_i_i_C(2) * t336, 0; t336 * r_i_i_C(1) + t335 * r_i_i_C(2) - t343 * r_i_i_C(3) + (t406 + (pkin(4) * t363 - t429) * t364) * t373 + (t397 - t428) * t415, (-t429 + t430) * t414 + t377 * t369 + t385, -pkin(5) * t402 + (-t364 * t389 - t407) * t362 + t385, t346 * r_i_i_C(3) + (-t345 * t366 - t347 * t409) * r_i_i_C(2) + (t345 * t370 - t347 * t411) * r_i_i_C(1), t351 * r_i_i_C(2) + (t386 * r_i_i_C(1) + t379 * r_i_i_C(2)) * t366 + ((-t379 - t405) * r_i_i_C(1) + t386 * r_i_i_C(2)) * t370, 0; 0, t376 - t406, t376, (-r_i_i_C(3) * t371 + t388 * t367) * t424 + ((r_i_i_C(1) * t366 + t427) * t367 * qJD(5) + (-t388 * t371 - t426) * qJD(4)) * t363, (r_i_i_C(1) * t383 + t396 * t427) * t362 + (t436 * r_i_i_C(1) + t433 * r_i_i_C(2)) * t363, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:23
	% EndTime: 2020-06-20 08:27:25
	% DurationCPUTime: 2.18s
	% Computational Cost: add. (1406->168), mult. (2206->277), div. (0->0), fcn. (2177->12), ass. (0->116)
	t622 = pkin(6) + r_i_i_C(3);
	t529 = qJ(2) + qJ(3);
	t526 = sin(t529);
	t536 = cos(qJ(5));
	t528 = qJD(2) + qJD(3);
	t537 = cos(qJ(4));
	t571 = t528 * t537 + qJD(5);
	t559 = t571 * t536;
	t527 = cos(t529);
	t532 = sin(qJ(4));
	t589 = qJD(6) * t532;
	t575 = t527 * t589;
	t531 = sin(qJ(5));
	t570 = qJD(5) * t537 + t528;
	t595 = qJD(4) * t532;
	t550 = t531 * t570 + t536 * t595;
	t629 = t527 * t550;
	t632 = t526 * t559 + t575 + t629;
	t530 = sin(qJ(6));
	t535 = cos(qJ(6));
	t566 = t535 * r_i_i_C(1) - t530 * r_i_i_C(2);
	t585 = t622 * t531;
	t545 = t566 * t536 + t585;
	t539 = cos(qJ(1));
	t599 = t539 * t532;
	t534 = sin(qJ(1));
	t601 = t534 * t537;
	t507 = t526 * t601 - t599;
	t608 = t527 * t534;
	t495 = t507 * t536 + t531 * t608;
	t598 = t539 * t537;
	t602 = t534 * t532;
	t506 = t526 * t602 + t598;
	t631 = -t495 * t530 - t506 * t535;
	t630 = t495 * t535 - t506 * t530;
	t597 = qJD(1) * t534;
	t604 = t528 * t539;
	t553 = -t526 * t604 - t527 * t597;
	t627 = t530 * r_i_i_C(1) + t535 * r_i_i_C(2);
	t600 = t536 * t537;
	t505 = -t526 * t531 + t527 * t600;
	t626 = qJD(6) * t505;
	t538 = cos(qJ(2));
	t616 = pkin(3) * qJD(2);
	t586 = t538 * t616;
	t619 = pkin(5) * t526;
	t620 = pkin(4) * t527;
	t625 = t528 * (-t619 + t620) + t586;
	t624 = -t531 * t595 + t536 * t570;
	t509 = t526 * t598 + t602;
	t593 = qJD(4) * t539;
	t573 = t537 * t593;
	t574 = t534 * t595;
	t605 = t528 * t534;
	t581 = t527 * t605;
	t493 = qJD(1) * t509 - t526 * t574 + t537 * t581 - t573;
	t596 = qJD(1) * t539;
	t578 = t527 * t596;
	t582 = t526 * t605;
	t591 = qJD(5) * t534;
	t623 = t531 * (t527 * t591 + t493) - t536 * (-qJD(5) * t507 + t578 - t582);
	t621 = pkin(3) * t538;
	t594 = qJD(4) * t537;
	t552 = t532 * t596 + t534 * t594;
	t492 = t552 * t526 - t537 * t597 + (t581 - t593) * t532;
	t615 = t492 * t530;
	t614 = t492 * t535;
	t611 = t526 * t528;
	t610 = t526 * t536;
	t609 = t527 * t531;
	t607 = t527 * t539;
	t606 = t528 * t532;
	t603 = t532 * t536;
	t592 = qJD(5) * t531;
	t590 = qJD(6) * t530;
	t588 = qJD(6) * t535;
	t587 = qJD(6) * t536;
	t583 = t527 * t604;
	t580 = t526 * t599;
	t572 = t493 * t536 - t531 * t582;
	t567 = -pkin(4) * t526 - pkin(5) * t527;
	t491 = (-qJD(4) * t526 + qJD(1)) * t599 + (t583 + (-qJD(1) * t526 + qJD(4)) * t534) * t537;
	t565 = qJD(5) * t607 + t491;
	t487 = t571 * t610 + t629;
	t563 = -t487 - t575;
	t562 = -t527 * t559 + (t550 + t589) * t526;
	t558 = t571 * t531;
	t486 = t526 * t558 - t624 * t527;
	t500 = t505 * t539;
	t542 = -t528 * t580 + qJD(6) * t500 + (-t532 * t597 + t573) * t527;
	t555 = t537 * t609 + t610;
	t557 = t505 * t597 + t632 * t539;
	t561 = (-t557 * t530 + t542 * t535) * r_i_i_C(2) + (t542 * t530 + t557 * t535) * r_i_i_C(1) + t597 * t619 + t622 * (t486 * t539 + t555 * t597);
	t543 = t527 * t552 - t532 * t582 + t534 * t626;
	t556 = -qJD(1) * t500 + t632 * t534;
	t560 = (-t556 * t530 + t543 * t535) * r_i_i_C(2) + (t543 * t530 + t556 * t535) * r_i_i_C(1) + pkin(4) * t578 + t622 * (t486 * t534 - t555 * t596);
	t494 = -t507 * t531 + t536 * t608;
	t533 = sin(qJ(2));
	t554 = qJD(1) * (-t533 * pkin(3) - pkin(2) + t567);
	t548 = qJD(6) * (-t526 * t600 - t609) - t526 * t594 - t527 * t606;
	t551 = -t528 * t620 + (t530 * t562 + t535 * t548) * r_i_i_C(2) + (t530 * t548 - t535 * t562) * r_i_i_C(1) + pkin(5) * t611 + t622 * (t624 * t526 + t527 * t558);
	t547 = -t526 * t606 + t527 * t594 + t626;
	t546 = t528 * t567 - t533 * t616;
	t544 = (-t566 * t531 + t536 * t622) * qJD(5);
	t540 = -t627 * t587 + t544;
	t508 = t580 - t601;
	t498 = t509 * t536 + t531 * t607;
	t497 = -t509 * t531 + t536 * t607;
	t490 = qJD(1) * t506 - t526 * t573 - t532 * t583 - t574;
	t475 = -qJD(5) * t494 - t531 * t578 - t572;
	t473 = -t507 * t592 + (t531 * t596 + t536 * t591) * t527 + t572;
	t471 = -t509 * t592 + t553 * t531 + t565 * t536;
	t470 = -t565 * t531 + (-qJD(5) * t509 + t553) * t536;
	t463 = t471 * t535 + t490 * t530 + (-t498 * t530 - t508 * t535) * qJD(6);
	t462 = t471 * t530 - t490 * t535 + (t498 * t535 - t508 * t530) * qJD(6);
	t1 = [(-t475 * t535 - t615) * r_i_i_C(1) + (t475 * t530 - t614) * r_i_i_C(2) + t622 * t623 + (t631 * r_i_i_C(1) - t630 * r_i_i_C(2)) * qJD(6) - t625 * t534 + t539 * t554, (-t620 - t621) * t597 + t546 * t539 + t561, pkin(4) * t553 - pkin(5) * t583 + t561, (t491 * t530 + t509 * t588) * r_i_i_C(1) + (t491 * t535 - t509 * t590) * r_i_i_C(2) - t545 * t490 + t540 * t508, -t622 * t471 + (t470 * t530 + t497 * t588) * r_i_i_C(2) + (-t470 * t535 + t497 * t590) * r_i_i_C(1), t462 * r_i_i_C(1) + t463 * r_i_i_C(2); -t463 * r_i_i_C(1) + t462 * r_i_i_C(2) + t622 * t470 + t534 * t554 + t625 * t539, (-t619 + t621) * t596 + t546 * t534 + t560, -pkin(4) * t582 + (-t526 * t596 - t581) * pkin(5) + t560, (t493 * t530 + t507 * t588) * r_i_i_C(1) + (t493 * t535 - t507 * t590) * r_i_i_C(2) + t545 * t492 + t540 * t506, -t622 * t473 + (t494 * t588 - t530 * t623) * r_i_i_C(2) + (t494 * t590 + t535 * t623) * r_i_i_C(1), (t473 * t530 + t614) * r_i_i_C(1) + (t473 * t535 - t615) * r_i_i_C(2) + (t630 * r_i_i_C(1) + t631 * r_i_i_C(2)) * qJD(6); 0, t551 - t586, t551, ((-t530 * t537 - t535 * t603) * r_i_i_C(1) + (t530 * t603 - t535 * t537) * r_i_i_C(2) - t532 * t585) * t611 + ((qJD(4) * t545 + qJD(6) * t566) * t537 + (t544 + t627 * (-qJD(4) - t587)) * t532) * t527, t622 * t487 + (t486 * t530 - t555 * t588) * r_i_i_C(2) + (-t486 * t535 - t555 * t590) * r_i_i_C(1), (r_i_i_C(1) * t547 + r_i_i_C(2) * t563) * t535 + (r_i_i_C(1) * t563 - r_i_i_C(2) * t547) * t530;];
	JaD_transl = t1;
end