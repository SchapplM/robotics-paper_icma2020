% Calculate inertial parameters regressor of joint inertia matrix time derivative for
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
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = CloosQRC350DE_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:36
% EndTime: 2020-06-23 21:08:03
% DurationCPUTime: 25.20s
% Computational Cost: add. (2330->549), mult. (6323->916), div. (0->0), fcn. (4861->10), ass. (0->398)
t194 = pkin(7) * qJD(5) - qJD(6);
t215 = sin(qJ(5));
t220 = cos(qJ(4));
t212 = t220 ^ 2;
t193 = t212 - 0.1e1 / 0.2e1;
t219 = cos(qJ(5));
t290 = t219 * t193;
t218 = sin(qJ(2));
t221 = cos(qJ(3));
t388 = t218 * t221;
t217 = sin(qJ(3));
t222 = cos(qJ(2));
t389 = t217 * t222;
t114 = t388 + t389;
t298 = t114 * t220 / 0.2e1;
t385 = t222 * t221;
t391 = t217 * t218;
t468 = t385 - t391;
t490 = t194 * (t215 * t298 - t290 * t468);
t491 = 0.2e1 * t215;
t436 = pkin(6) * t468;
t205 = t219 + 0.1e1;
t206 = t219 - 0.1e1;
t216 = sin(qJ(4));
t363 = qJD(5) * t220;
t309 = t215 * t363;
t371 = qJD(4) * t216;
t322 = t206 * t371;
t399 = t206 * t216;
t481 = t309 * t399 + ((t309 + t322) * t216 - t206 * t212 * qJD(4)) * t205;
t489 = t481 * t436;
t195 = pkin(7) * qJ(5) - qJ(6);
t174 = sin(t195);
t168 = t174 ^ 2;
t488 = t168 * t490;
t209 = t216 ^ 2;
t464 = qJD(2) + qJD(3);
t53 = t464 * t114;
t426 = t220 * t53;
t22 = qJD(4) * t468 * (t209 - t212) + t216 * t426;
t487 = 2 * pkin(4);
t486 = -0.2e1 * t205;
t175 = cos(t195);
t408 = t174 * t175;
t485 = 0.2e1 * t408;
t211 = t219 ^ 2;
t197 = t211 + 0.1e1;
t483 = t468 * t197;
t365 = qJD(5) * t217;
t312 = t215 * t365;
t156 = pkin(6) * t312;
t433 = pkin(6) * t219;
t181 = pkin(5) + t433;
t374 = qJD(3) * t217;
t189 = pkin(4) * t374;
t373 = qJD(3) * t221;
t245 = t181 * t373 - t156 + t189;
t190 = pkin(4) * t373;
t367 = qJD(5) * t215;
t188 = pkin(6) * t367;
t279 = t221 * t188;
t61 = -t181 * t374 + t190 - t279;
t480 = t61 * t218 + t245 * t222;
t224 = 0.2e1 * pkin(6);
t213 = t221 ^ 2;
t479 = 0.4e1 * t213;
t449 = 0.2e1 * t221;
t448 = 0.4e1 * t221;
t342 = pkin(5) * t373;
t117 = -t189 - t342;
t118 = -pkin(5) * t374 + t190;
t199 = pkin(5) * t217;
t201 = pkin(4) * t221;
t430 = pkin(3) + t201;
t131 = -t199 + t430;
t200 = pkin(4) * t217;
t437 = pkin(5) * t221;
t466 = -t437 - t200;
t58 = t131 * t222 + t218 * t466;
t28 = qJD(2) * t58 + t117 * t218 + t118 * t222;
t252 = -t131 * t218 + t222 * t466;
t56 = pkin(2) - t252;
t478 = -t28 * t219 + t56 * t367;
t214 = t222 ^ 2;
t444 = t479 - 0.2e1;
t476 = t444 * t214;
t226 = pkin(6) ^ 2;
t364 = qJD(5) * t219;
t304 = t215 * t364;
t266 = t226 * t304;
t127 = t212 * t266;
t445 = t211 - 0.1e1;
t171 = t445 * t226;
t369 = qJD(4) * t220;
t303 = t216 * t369;
t269 = t171 * t303;
t246 = t127 + t269;
t475 = -t171 * t212 - t211 * t226;
t207 = (pkin(4) ^ 2) - pkin(5) ^ 2;
t394 = t215 * t220;
t346 = pkin(6) * t394;
t474 = -0.2e1 * pkin(4) * t346 + t207;
t473 = t217 ^ 2 - t213;
t400 = t205 * t206;
t335 = t212 * t400;
t472 = -t211 - t335;
t186 = pkin(2) * t218 + pkin(3);
t375 = qJD(2) * t222;
t329 = t217 * t375;
t471 = pkin(2) * t329 + t186 * t373;
t151 = -pkin(4) + t346;
t310 = t219 * t363;
t321 = t215 * t371;
t238 = -t310 + t321;
t89 = t238 * pkin(6);
t470 = t151 * t374 + t89 * t221;
t323 = t220 * t373;
t469 = t217 * t371 - t323;
t347 = pkin(5) * t394;
t431 = t219 * pkin(4);
t465 = (t347 + t431) * qJD(5);
t395 = t215 * t219;
t463 = t114 * t395 - t220 * t483;
t343 = pkin(3) * t374;
t153 = t215 * t343;
t203 = pkin(3) * t221;
t185 = t203 + pkin(4);
t324 = t219 * t373;
t281 = pkin(3) * t324;
t202 = pkin(3) * t217;
t182 = t202 - pkin(5);
t307 = t182 * t367;
t461 = -(t281 - t307) * t220 - t185 * t364 + t153;
t366 = qJD(5) * t216;
t313 = t215 * t366;
t370 = qJD(4) * t219;
t460 = t220 * (-t194 + t370) - t313;
t225 = -0.2e1 * pkin(5);
t459 = -0.4e1 * pkin(6);
t178 = t218 * pkin(3) + pkin(2);
t140 = t178 * t217;
t458 = 0.2e1 * t181 * t218 - 0.2e1 * t140;
t172 = t433 + pkin(5) / 0.2e1;
t457 = 0.2e1 * t172;
t187 = pkin(3) * t375;
t293 = 0.2e1 * t187;
t456 = -0.4e1 * t213;
t454 = -0.2e1 * t218;
t453 = -0.2e1 * t219;
t452 = 0.2e1 * t219;
t451 = -0.2e1 * t220;
t450 = -0.2e1 * t221;
t447 = -0.2e1 * t222;
t446 = pkin(6) * t53;
t442 = pkin(4) * t218;
t441 = pkin(5) * t213;
t440 = pkin(5) * t215;
t439 = pkin(5) * t218;
t438 = pkin(5) * t219;
t208 = t215 ^ 2;
t435 = pkin(6) * t208;
t434 = pkin(6) * t211;
t432 = t215 * pkin(4);
t429 = t468 * t53;
t300 = t221 * t375;
t54 = -t222 * t373 + t391 * t464 - t300;
t428 = t215 * t54;
t180 = pkin(5) + 0.2e1 * t433;
t99 = -t180 * t217 + t430;
t427 = t218 * t99;
t425 = t222 * t99;
t234 = t238 * pkin(4);
t341 = pkin(5) * t367;
t423 = (t234 + t341) * t213;
t278 = pkin(6) * t304;
t255 = t212 * t278;
t179 = t445 * pkin(6);
t268 = t179 * t303;
t421 = 0.4e1 * t255 + 0.4e1 * t268;
t104 = t180 * t221 + t200;
t420 = t104 * t218;
t419 = t104 * t222;
t105 = t181 * t221 + t200;
t418 = t105 * t218;
t417 = t105 * t222;
t416 = t468 * t216;
t414 = t114 * t219;
t411 = t151 * t181;
t410 = t151 * t221;
t407 = t174 * t194;
t406 = t175 * t194;
t405 = t179 * t212;
t404 = t182 * t215;
t403 = t182 * t219;
t402 = t194 * t215;
t401 = t194 * t216;
t397 = t213 * t218;
t396 = t215 * t216;
t393 = t215 * t226;
t392 = t216 * t219;
t390 = t217 * t219;
t191 = pkin(3) * t373;
t116 = -t188 - t191;
t387 = t219 * t116;
t386 = t219 * t220;
t384 = 0.2e1 * t246;
t184 = t203 + t487;
t110 = t184 * t321;
t126 = t220 * t153;
t383 = t110 + t126;
t301 = t218 * t367;
t285 = 0.2e1 * t301;
t328 = t219 * t375;
t382 = (t285 - 0.2e1 * t328) * pkin(4);
t280 = pkin(4) * t321;
t286 = -0.2e1 * t310;
t381 = pkin(4) * t286 + 0.2e1 * t280;
t380 = 0.2e1 * t471;
t154 = pkin(6) * t324;
t379 = -t156 + t154;
t378 = (t312 - t324) * pkin(4);
t348 = pkin(4) * t394;
t164 = -0.2e1 * t348;
t377 = t164 - 0.2e1 * t405;
t173 = t199 - pkin(3) / 0.2e1;
t183 = t202 + t225;
t376 = qJD(2) * t218;
t368 = qJD(5) * t114;
t362 = qJD(5) * t221;
t196 = -0.2e1 * t433;
t361 = t471 * t453;
t353 = -0.2e1 * t391;
t112 = pkin(5) * t353 + t178;
t359 = t112 * t452;
t358 = -0.4e1 * t410;
t357 = 0.2e1 * t403;
t356 = -0.8e1 * t397;
t355 = -0.4e1 * t397;
t354 = -0.2e1 * t394;
t352 = 0.2e1 * qJD(5);
t351 = pkin(4) * t390;
t350 = t217 * t201;
t349 = pkin(4) * t386;
t170 = pkin(6) * t390;
t345 = -0.2e1 * t374;
t344 = 0.2e1 * t367;
t133 = -t199 + t186;
t339 = t174 * t406;
t338 = t174 * t402;
t337 = t175 * t402;
t336 = t185 * t386;
t334 = t220 * t400;
t333 = t219 * t393;
t332 = t184 * t386;
t331 = t217 * t385;
t330 = t213 * t375;
t327 = qJD(3) * t396;
t326 = t217 * t373;
t325 = t218 * t373;
t320 = t215 * t369;
t318 = t216 * t370;
t317 = t219 * t369;
t316 = t114 * t364;
t315 = t208 * t366;
t314 = t468 * t367;
t311 = t216 * t364;
t302 = t218 * t371;
t299 = t218 * t375;
t297 = -0.4e1 * t220 * t200;
t296 = pkin(3) * t214 - t186;
t128 = t186 * t217 - pkin(5);
t294 = -0.2e1 * t342;
t267 = t468 * t318;
t292 = -t267 * t451 - t114 * t321 / 0.2e1 - t54 * t394 / 0.2e1 + t298 * t364;
t291 = 0.2e1 * t339;
t106 = t181 * t217 - t201;
t289 = -0.4e1 * t330;
t288 = -0.2e1 * t326;
t287 = t128 * t344;
t284 = -0.2e1 * t299;
t277 = t381 + t421;
t275 = t208 * t339;
t274 = t218 * t331;
t273 = 0.2e1 * t205 + 0.2e1 * t206;
t272 = t468 ^ 2 * t303;
t271 = t185 * t321;
t270 = t209 * t304;
t80 = t468 * t310;
t265 = t219 * t54 + t53 * t394 - t80;
t264 = 0.2e1 * t281;
t145 = t329 * t225;
t263 = t53 * t399 * t486;
t192 = pkin(2) * t375;
t262 = t192 - t342;
t96 = t215 * t221 + t217 * t386;
t97 = t215 * t217 - t221 * t386;
t260 = (-t218 * t97 + t222 * t96) * qJD(2) + ((t220 * t374 + t221 * t371 + t365) * t222 + t218 * (t362 - t469)) * t219 + (-t401 + t215 * (qJD(3) + t363)) * t468;
t259 = (t218 * t294 + t145 + t187) * t219;
t258 = t334 * t446;
t136 = pkin(4) * t222 - t439;
t138 = pkin(5) * t222 + t442;
t257 = t212 * t285;
t98 = pkin(3) - t106;
t254 = -t218 * t98 - t417;
t42 = t222 * t98 - t418;
t251 = t216 * (t194 * t219 - qJD(4));
t250 = pkin(3) * t329 - pkin(5) * t375 + t178 * t373;
t249 = t342 + t379;
t247 = t349 - t440;
t121 = t348 + t438;
t244 = t216 * t53 - t369 * t468;
t243 = t371 * t468 + t426;
t242 = t80 * t491 - t54 * t395 + (-t208 + t211) * t368;
t241 = t474 + t475;
t139 = pkin(6) * t280;
t240 = t139 + t246;
t107 = t466 * qJD(3) * pkin(3);
t94 = t311 + t320;
t236 = pkin(6) * t126 + t107;
t37 = t218 * t96 + t222 * t97;
t235 = -t194 * t37 + t244;
t12 = (t363 * t468 - t54) * t215 + (t243 + t368) * t219;
t48 = t114 * t215 - t386 * t468;
t6 = t48 * t313 + (-t12 * t216 - t369 * t48) * t219;
t232 = 0.2e1 * t301 - 0.2e1 * t328;
t231 = (t183 * t215 - t332) * qJD(5) - t281;
t230 = 0.2e1 * t472 * t226 + 0.2e1 * t474;
t229 = -t281 + t271 + (-t336 + t404) * qJD(5);
t176 = -0.2e1 * t191;
t169 = t175 ^ 2;
t163 = t431 * t454;
t161 = -0.2e1 * t303;
t160 = 0.2e1 * t303;
t159 = -0.2e1 * t304;
t147 = 0.4e1 * t278;
t146 = 0.2e1 * t266;
t143 = 0.2e1 * t156;
t130 = t181 - t202;
t120 = t347 - t431;
t119 = t191 + 0.2e1 * t188;
t111 = t218 * t215 * t297;
t108 = t170 + t173;
t103 = t140 - t439;
t102 = t170 + t199 / 0.2e1 - pkin(3) / 0.4e1;
t100 = t217 * t457 - t201;
t95 = t172 * t221 + t200 / 0.2e1;
t92 = 0.2e1 * t107;
t91 = -0.2e1 * t170 + t133;
t85 = t94 * pkin(7);
t84 = t172 * t394 - t431 / 0.2e1;
t83 = t468 * t394;
t78 = -t173 * t394 + t351;
t77 = -t180 * t218 + t140;
t73 = t170 * t454 + t112;
t72 = t174 * t220 - t175 * t392;
t71 = t174 * t392 + t175 * t220;
t69 = t293 + 0.4e1 * (-t325 - t329) * pkin(5);
t68 = 0.2e1 * t211 * t303 - 0.2e1 * t270;
t67 = 0.2e1 * t208 * t303 + 0.2e1 * t270;
t63 = -t180 * t373 + t143 - t189;
t62 = -t180 * t374 + t190 - 0.2e1 * t279;
t60 = -t211 * t366 - t215 * t317 + t315;
t57 = 0.4e1 * t121 + 0.4e1 * t405 + 0.4e1 * t434;
t55 = -0.8e1 * t102 * t394 + 0.4e1 * t351;
t52 = t56 * t219;
t51 = t111 + t359;
t50 = 0.2e1 * t121 * t224 - 0.2e1 * t207 - 0.2e1 * t475;
t49 = -t83 - t414;
t47 = t106 * t222 + t418;
t43 = t100 * t222 + t420;
t41 = pkin(2) + t419 + t427;
t40 = pkin(2) - t254;
t38 = t95 * t454 + t425;
t33 = -0.2e1 * t114 * t54;
t32 = -t364 * t408 + (t168 - t169) * t402;
t29 = t335 * t391 * t459 + t73 * t452 + t111;
t26 = -t174 * t416 + t175 * t37;
t25 = t174 * t37 + t175 * t416;
t24 = -0.2e1 * t209 * t429 + 0.2e1 * t272;
t21 = t174 * t460 + t175 * t251;
t20 = t174 * t251 - t175 * t460;
t19 = 0.2e1 * pkin(5) * t178 + pkin(6) * t359 + t230 * t391;
t18 = 0.2e1 * t56 * t28;
t17 = 0.2e1 * t72 * t20;
t16 = 0.2e1 * t71 * t21;
t15 = t175 * t21 - t71 * t407;
t14 = t174 * t20 + t72 * t406;
t13 = -t114 * t367 - t321 * t468 - t265;
t10 = -t71 * t337 + (-t21 * t215 - t364 * t71) * t174;
t9 = -t72 * t338 + (t20 * t215 + t364 * t72) * t175;
t8 = -0.2e1 * t49 * t13;
t7 = -t219 * t265 + (-t267 + (-t83 - 0.2e1 * t414) * qJD(5)) * t215;
t5 = -t49 * t311 + (t13 * t216 - t369 * t49) * t215;
t4 = t174 * t235 + t175 * t260;
t3 = t174 * t260 - t175 * t235;
t2 = t216 * t463 * t291 + (t193 * t314 + t290 * t53 + t292) * t485 - 0.2e1 * t488 + t22 + (-t463 * t369 - t216 * (t197 * t243 + t242) + 0.2e1 * t490) * t169;
t1 = 0.2e1 * t488 + (-t463 * t401 - (t219 * t53 + t314) * t193 - t292) * t485 + (t216 * (t211 * t426 + t242 + t426) - 0.2e1 * t490 + (t209 * t483 + t220 * t463) * qJD(4)) * t169 + t6;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, 0, 0, 0.2e1 * t299, 0, 0, 0.2e1 * t192, 0, 0, 0, -0.2e1 * t429, (-t476 + (t444 * t218 + 0.8e1 * t331) * t218) * qJD(2) + (0.8e1 * t274 + t473 * (0.4e1 * t214 - 0.2e1)) * qJD(3), 0, t33, 0, 0, 0.2e1 * (t178 * t385 + t217 * t296) * qJD(3) + 0.2e1 * (pkin(2) * t385 - t178 * t391 + (t214 * t217 + 0.2e1 * t218 * t385) * pkin(3)) * qJD(2), 0.2e1 * (-t178 * t389 + t221 * t296) * qJD(3) + 0.2e1 * (-pkin(2) * t389 - t178 * t388 + (t214 * t221 + t222 * t353) * pkin(3)) * qJD(2), 0, t178 * t293, -0.2e1 * t212 * t429 - 0.2e1 * t272, 0, 0, t24, 0, t33, 0, 0, t69 * t385 + 0.2e1 * t191 * t214 + (pkin(4) * t476 + (-t444 * t442 - 0.4e1 * (t182 + 0.2e1 * t350 + 0.2e1 * t441) * t222 + t112 * t450) * t218) * qJD(2) + ((t112 * t447 + (-0.8e1 * t214 + 0.4e1) * t437) * t217 + (-0.2e1 * t473 * t214 - 0.4e1 * t274 + t473) * t487) * qJD(3) - t380, t18, 0.2e1 * t48 * t12, 0, 0, t8, 0, t24, 0, (t264 - 0.2e1 * t307 - 0.4e1 * t423 + (t173 * t238 - t378) * t448 + 0.4e1 * (-t213 * t347 + (t121 * t450 - t78) * t217) * qJD(3) + t381) * t214 + (t121 * t479 + t78 * t448 + t164 + t357) * t284 + (t120 * t289 + (0.8e1 * t120 * t388 - t51) * t374 + (t69 * t221 + (-0.2e1 * t103 + 0.4e1 * (-t350 - t441) * t218) * t363) * t219 + ((-0.2e1 * t112 * qJD(5) + t297 * t375) * t221 + t250 * t451 + 0.2e1 * t103 * t371 + (-pkin(5) * t371 * t456 + (qJD(5) * t456 + t448 * t469) * pkin(4)) * t218) * t215 + t382) * t222 - (t103 * t354 + t120 * t355 + t51 * t221 + t163) * t376 + 0.2e1 * t423 + 0.4e1 * t121 * t326 + (t133 * t238 - t262 * t394 + t378) * t449 + (-t133 * t394 - t351) * t345 + t287 + t361, 0, t18, 0.2e1 * t26 * t4, 0, 0, 0.2e1 * t25 * t3, 0, t8, 0, 0, (((-0.8e1 * t179 * t220 - 0.4e1 * t432) * t371 + (0.4e1 * t349 + (-0.4e1 * pkin(5) + (-0.8e1 * t212 - 0.8e1) * t433) * t215) * qJD(5)) * t213 + t57 * t288 + ((t342 / 0.2e1 + t379) * t354 + 0.2e1 * t238 * t102 - t378) * t448 - t55 * t374 + t130 * t344 - 0.2e1 * t387 + t277) * t214 + (t130 * t453 + t213 * t57 + t221 * t55 + t377) * t284 + ((-t172 * t321 + (t432 / 0.2e1 + (t172 * t219 - t435) * t220) * qJD(5)) * t356 + (-t73 * t367 + t259 + (t472 * t224 + t164) * t325 + ((t218 * t286 + (-t220 * t375 + t302) * t491) * pkin(4) + (t206 * t257 + t219 * t232 + (t257 + (-0.2e1 * t212 * t375 + 0.4e1 * t220 * t302) * t206) * t205) * pkin(6)) * t217) * t449 - t29 * t374 + (pkin(6) * t232 + t250) * t354 + t382 + (0.16e2 * t217 * t325 - 0.8e1 * t330) * t84 + (t286 + 0.2e1 * t321) * t77) * t222 - (t221 * t29 + t354 * t77 + t356 * t84 + t163) * t376 + (t147 + t277 + 0.2e1 * t341) * t213 + (t377 - 0.2e1 * t434 - 0.2e1 * t438) * t288 + (-t91 * t310 + (-(t143 - 0.2e1 * t154 + t262) * t220 + t91 * t371) * t215 + t378) * t449 + (-t394 * t91 - t351) * t345 - 0.4e1 * t255 - 0.4e1 * t268 - t380 * t219 + t128 * t215 * t352, ((-t240 - t266) * t479 + t50 * t288 + t249 * t358 + 0.2e1 * t139 + t146 + 0.2e1 * pkin(3) * t342 + 0.4e1 * t470 * t108 + (t264 + (t247 * t479 - 0.2e1 * t349 - 0.2e1 * t404) * qJD(5)) * pkin(6) + t384) * t214 + (pkin(6) * t357 + t108 * t358 + t213 * t50 + (-pkin(3) + 0.2e1 * t199) * pkin(3) + t241) * t284 + (t289 * t411 + (pkin(5) * t187 + pkin(6) * t259 - t112 * t188) * t449 - t89 * t458 - 0.2e1 * t151 * (-pkin(6) * t328 + t250) + (t230 * t213 * qJD(3) + t444 * t151 * t188 - t181 * t89 * t456) * t218 + (-t19 * qJD(3) + t230 * t300 + ((0.8e1 * qJD(3) * t411 + 0.4e1 * (pkin(6) * t432 + t226 * t334) * t371) * t221 + (t349 * t459 + (t212 * t273 + 0.4e1 * t219) * t393) * t362) * t218) * t217) * t222 - (t19 * t221 + (t181 * t355 + t458) * t151) * t376 + 0.2e1 * ((-pkin(6) * t247 + t333) * qJD(5) + t240) * t213 + (pkin(5) * t196 + t241) * t288 - 0.2e1 * (t192 - t249) * t410 - 0.2e1 * t127 - 0.2e1 * t269 + pkin(6) * t361 + pkin(6) * t287 + t186 * t294 + 0.2e1 * t470 * (-t170 + t133) + (t145 + t293) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t376, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t54, 0, 0, 0, -pkin(3) * t376, 0, -t22, 0, 0, t22, 0, 0, 0, 0, 0, 0, t6, 0, 0, t5, 0, t22, 0, t58 * t311 + ((qJD(2) * t252 + t117 * t222 - t118 * t218) * t216 + t58 * t369) * t215, 0, 0, t2, 0, 0, t1, 0, t5, 0, 0, t216 * ((t63 * t222 + (-t279 - t172 * t374 + t190 / 0.2e1) * t454 + (t95 * t447 - t427) * qJD(2)) * t215 + t38 * t364 + t322 * t436 * t486) + (qJD(4) * t38 * t215 + (t263 - (-0.2e1 * qJD(4) * t334 + t273 * t313) * t468) * pkin(6)) * t220, (t42 * t320 + t216 * ((qJD(2) * t254 - t480) * t215 + t42 * t364 - t258) - t489) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t343, t176, 0, 0, t160, 0, 0, t161, 0, 0, 0, 0, t176, t92, t68, 0, 0, t67, 0, t161, 0, 0.2e1 * t126 + 0.2e1 * t229, 0, t92, t17, 0, 0, t16, 0, t67, 0, 0, 0.2e1 * t271 + 0.2e1 * t387 + 0.2e1 * t126 + (-t215 * t130 - t336) * t352 + t421, 0.2e1 * pkin(6) * t229 + 0.2e1 * t236 + 0.2e1 * t246 - 0.2e1 * t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t54, 0, 0, 0, 0, 0, -t22, 0, 0, t22, 0, 0, 0, 0, 0, 0, t6, 0, 0, t5, 0, t22, 0, -t94 * (-t136 * t221 + t138 * t217) + (-qJD(2) * t396 - t327) * (t136 * t217 + t138 * t221), 0, 0, t2, 0, 0, t1, 0, t5, 0, 0, -(((t373 * t457 - 0.2e1 * t156 + t189) * t222 + t62 * t218 + (-t100 * t218 + t419) * qJD(2)) * t215 + t43 * t364) * t216 - t43 * t320 + (t220 * t263 - 0.2e1 * t468 * t481) * pkin(6), (-t47 * t320 + (-((-t106 * t218 + t417) * qJD(2) + t480) * t215 - t47 * t364 - t258) * t216 - t489) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t343, -t191, 0, 0, t160, 0, 0, t161, 0, 0, 0, 0, -t191, t107, t68, 0, 0, t67, 0, t161, 0, t231 + t383, 0, t107, t17, 0, 0, t16, 0, t67, 0, 0, -t219 * t119 + (-t332 + t215 * (t196 + t183)) * qJD(5) + t383 + t421, -0.2e1 * t266 + (t231 + t110) * pkin(6) + t236 + t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, 0, 0, t161, 0, 0, 0, 0, 0, 0, t68, 0, 0, t67, 0, t161, 0, 0.2e1 * t234 - 0.2e1 * t341, 0, 0, t17, 0, 0, t16, 0, t67, 0, 0, (t225 - 0.4e1 * t433) * t367 + t277, 0.2e1 * (-t333 + (-t349 - t440) * pkin(6)) * qJD(5) + 0.2e1 * t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, 0, 0, 0, 0, t12 * t215 + t364 * t48, 0, 0, t7, 0, 0, 0, -t56 * t310 + (-t220 * t28 + t371 * t56) * t215, 0, 0, -t26 * t338 + (t215 * t4 + t26 * t364) * t175, 0, 0, t25 * t337 + (t215 * t3 + t25 * t364) * t174, 0, t7, 0, 0, -t41 * t310 + (-(t218 * t63 + t222 * t62 + (-t420 + t425) * qJD(2)) * t220 + t41 * t371 + (0.4e1 * t316 - 0.2e1 * t428) * pkin(6)) * t215, (-t40 * t310 + (-(qJD(2) * t42 - t218 * t245 + t222 * t61) * t220 + t40 * t371 + (0.2e1 * t316 - t428) * pkin(6)) * t215) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, -t60, 0, 0, 0, t182 * t94 + t327 * t203, 0, 0, t9, 0, 0, t10, 0, -t60, 0, 0, t119 * t396 + t94 * (t196 + t182), (-t116 * t396 - t130 * t94) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, -t60, 0, 0, 0, -t94 * pkin(5), 0, 0, t9, 0, 0, t10, 0, -t60, 0, 0, -t180 * t320 + (-t180 * t219 + 0.2e1 * t435) * t366, -pkin(6) * t181 * t94 + t226 * t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t304, 0, 0, t159, 0, 0, 0, 0, 0, 0, 0.2e1 * t169 * t304 - 0.2e1 * t275, 0, 0, 0.2e1 * t168 * t304 + 0.2e1 * t275, 0, t159, 0, 0, t147, t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, 0, t216 * t478 - t56 * t317, 0, 0, t174 * t4 + t26 * t406, 0, 0, -t175 * t3 + t25 * t407, 0, -t13 * pkin(7), 0, 0, -(t52 + 0.2e1 * t436) * t369 - t216 * (-t478 - 0.2e1 * t446), (-(t52 + t436) * t369 - t216 * (-t478 - t446)) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t371, 0, t220 * t307 + t153 + (-pkin(3) * t323 - qJD(5) * t185 + t182 * t371) * t219, 0, 0, t14, 0, 0, t15, 0, -t85, 0, 0, -(t224 - t403) * t371 + t461, -(-(-pkin(6) + t403) * t371 - t461) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t371, 0, -pkin(4) * t364 + (-t309 - t318) * pkin(5), 0, 0, t14, 0, 0, t15, 0, -t85, 0, 0, -(t224 + t438) * t371 - t465, -pkin(6) * (-(-pkin(6) - t438) * t371 + t465); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, t32, 0, pkin(7) * t367, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, 0, 0, -0.2e1 * t339, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t367, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
