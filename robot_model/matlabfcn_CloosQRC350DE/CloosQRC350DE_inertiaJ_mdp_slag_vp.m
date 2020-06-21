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
% MDP [36x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350DE_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1),zeros(36,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [36 1]), ...
  'CloosQRC350DE_inertiaJ_mdp_slag_vp: MDP has to be [36x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:39
% EndTime: 2020-06-19 21:39:40
% DurationCPUTime: 0.61s
% Computational Cost: add. (1693->390), mult. (3282->577), div. (0->0), fcn. (3280->10), ass. (0->226)
unknown=NaN(21,1);
t1 = cos(qJ(4));
t2 = cos(qJ(3));
t3 = cos(qJ(2));
t5 = sin(qJ(3));
t6 = sin(qJ(2));
t8 = -t2 * t3 + t5 * t6;
t9 = t8 ^ 2;
t11 = sin(qJ(4));
t15 = t1 * t8;
t18 = t2 * t6 + t5 * t3;
t22 = t11 * t8;
t26 = t6 * pkin(3);
t29 = -t18 * pkin(4) + t8 * pkin(5) - pkin(2) - t26;
t30 = t1 * t29;
t34 = t11 * t29;
t38 = cos(qJ(5));
t39 = t38 * t1;
t41 = sin(qJ(5));
t43 = t41 * t18 + t39 * t8;
t45 = t8 * MDP(25);
t48 = t41 * t1;
t51 = t38 * t18 - t48 * t8;
t53 = t8 * MDP(26);
t57 = pkin(7) * qJ(5) - qJ(6);
t58 = cos(t57);
t60 = sin(t57);
t63 = -t60 * t11 * t8 - t58 * t43;
t67 = t58 * t11 * t8 - t60 * t43;
t89 = -0.2e1 * t1 * t9 * t11 * MDP(17) + 0.2e1 * t15 * t18 * MDP(18) - 0.2e1 * t22 * t18 * MDP(19) - 0.2e1 * t30 * t18 * MDP(21) + 0.2e1 * t34 * t18 * MDP(22) + 0.2e1 * t43 * t11 * t45 + 0.2e1 * t51 * t11 * t53 + MDP(1) + 0.2e1 * t63 * t67 * MDP(31) + 0.2e1 * t63 * t51 * MDP(32) + 0.2e1 * t67 * t51 * MDP(33) - 0.2e1 * t3 * t6 * MDP(3) + 0.2e1 * pkin(2) * t6 * MDP(7) + 0.2e1 * pkin(2) * t3 * MDP(8) + 0.2e1 * t8 * t18 * MDP(10);
t90 = -pkin(2) - t26;
t97 = t1 ^ 2;
t103 = t11 ^ 2;
t106 = t3 ^ 2;
t109 = t18 ^ 2;
t111 = t43 ^ 2;
t113 = t41 * t103;
t114 = t29 * t8;
t119 = t38 * t103;
t124 = t63 ^ 2;
t126 = t51 ^ 2;
t128 = t38 * t11;
t131 = -t22 * pkin(6) - t128 * t29;
t134 = t43 * pkin(6) + t30;
t136 = -t60 * t131 + t58 * t134;
t138 = t41 * t11;
t145 = -t58 * t131 - t60 * t134;
t151 = -0.2e1 * t90 * t18 * MDP(14) + 0.2e1 * t90 * t8 * MDP(15) + t97 * t9 * MDP(16) + 0.2e1 * t43 * t51 * MDP(24) + t103 * t9 * MDP(27) + t106 * MDP(2) + t9 * MDP(9) + t109 * MDP(20) + t111 * MDP(23) + 0.2e1 * (t113 * t114 - t30 * t51) * MDP(28) + 0.2e1 * (t119 * t114 + t30 * t43) * MDP(29) + t124 * MDP(30) + t126 * MDP(34) + 0.2e1 * (-t138 * t29 * t67 + t136 * t51) * MDP(35) + 0.2e1 * (t138 * t29 * t63 - t145 * t51) * MDP(36);
t154 = t8 * MDP(11);
t155 = t18 * MDP(12);
t159 = (t103 * t8 - t97 * t8) * MDP(17);
t160 = t5 * pkin(3);
t161 = t160 - pkin(5);
t162 = t11 * t161;
t164 = t2 * pkin(3);
t165 = t164 + pkin(4);
t166 = t165 * t11;
t170 = t1 * t161;
t172 = t165 * t1;
t177 = t43 * t41;
t180 = (t177 * t11 - t51 * t128) * MDP(24);
t184 = (t43 * t1 - t119 * t8) * MDP(25);
t188 = (t51 * t1 + t113 * t8) * MDP(26);
t191 = -t48 * t161 + t38 * t165;
t197 = t39 * t161;
t198 = t41 * t165;
t199 = t197 + t198;
t205 = t58 * t38;
t208 = -t60 * t1 + t205 * t11;
t210 = t60 * t38;
t213 = t58 * t1 + t210 * t11;
t216 = (t208 * t67 + t63 * t213) * MDP(31);
t217 = t6 * MDP(5) + t154 + t155 + t159 + (-t162 * t18 + t166 * t8) * MDP(21) + (-t170 * t18 + t172 * t8) * MDP(22) + t180 + t184 + t188 + (t191 * t11 * t8 - t162 * t51) * MDP(28) + (-t199 * t11 * t8 + t162 * t43) * MDP(29) + t216;
t222 = (t63 * t41 * t11 + t208 * t51) * MDP(32);
t227 = (t67 * t41 * t11 + t213 * t51) * MDP(33);
t228 = t1 * pkin(6);
t229 = t197 + t198 - t228;
t231 = t128 * pkin(6);
t232 = t162 - t231;
t234 = -t60 * t229 + t58 * t232;
t237 = t136 * t41 * t11;
t240 = -t138 * t29 * t213;
t245 = -t58 * t229 - t60 * t232;
t248 = t145 * t41 * t11;
t251 = t138 * t29 * t208;
t256 = t11 * t18 * MDP(18);
t258 = t1 * t18 * MDP(19);
t260 = t63 * t208 * MDP(30);
t262 = t15 * t11 * MDP(16);
t263 = t43 * t38;
t265 = t263 * t11 * MDP(23);
t267 = t15 * t11 * MDP(27);
t268 = t41 * t51;
t270 = t268 * t11 * MDP(34);
t271 = t222 + t227 + (-t191 * t67 + t234 * t51 + t237 + t240) * MDP(35) + (t191 * t63 - t245 * t51 - t248 + t251) * MDP(36) - t3 * MDP(4) - t256 - t258 + t260 - t262 - t265 + t267 + t270;
t273 = t164 * MDP(14);
t275 = t160 * MDP(15);
t277 = t103 * MDP(16);
t280 = 0.2e1 * t11 * t1 * MDP(17);
t285 = t38 ^ 2;
t287 = t285 * t103 * MDP(23);
t290 = 0.2e1 * t119 * t41 * MDP(24);
t293 = 0.2e1 * t128 * t1 * MDP(25);
t294 = 0.2e1 * t172 * MDP(21) - 0.2e1 * t166 * MDP(22) + MDP(13) + MDP(6) + 0.2e1 * t273 - 0.2e1 * t275 + t277 + t280 + t287 - t290 - t293;
t297 = 0.2e1 * t138 * t1 * MDP(26);
t298 = t97 * MDP(27);
t299 = t191 * t1;
t300 = t103 * t161;
t301 = t300 * t41;
t304 = t199 * t1;
t305 = t300 * t38;
t308 = t208 ^ 2;
t309 = t308 * MDP(30);
t312 = 0.2e1 * t208 * t213 * MDP(31);
t316 = 0.2e1 * t208 * t41 * t11 * MDP(32);
t320 = 0.2e1 * t213 * t41 * t11 * MDP(33);
t321 = t41 ^ 2;
t323 = t321 * t103 * MDP(34);
t325 = t234 * t41 * t11;
t326 = -t191 * t213;
t330 = t245 * t41 * t11;
t331 = t191 * t208;
t334 = t297 + t298 + 0.2e1 * (t299 - t301) * MDP(28) + 0.2e1 * (-t304 - t305) * MDP(29) + t309 + t312 + t316 + t320 + t323 + 0.2e1 * (t325 + t326) * MDP(35) + 0.2e1 * (-t330 + t331) * MDP(36);
t336 = t11 * pkin(5);
t338 = pkin(4) * t11;
t342 = t1 * pkin(5);
t344 = pkin(4) * t1;
t348 = t154 + t155 - t262 + t159 - t256 - t258 + (t336 * t18 + t338 * t8) * MDP(21) + (t342 * t18 + t344 * t8) * MDP(22) - t265 + t180 + t184;
t351 = t38 * pkin(4) + t48 * pkin(5);
t357 = t39 * pkin(5);
t358 = t41 * pkin(4);
t359 = -t357 + t358;
t365 = -t357 + t358 - t228;
t367 = -t336 - t231;
t369 = -t60 * t365 + t58 * t367;
t376 = -t58 * t365 - t60 * t367;
t381 = t188 + t267 + (t351 * t11 * t8 + t336 * t51) * MDP(28) + (-t359 * t11 * t8 - t336 * t43) * MDP(29) + t260 + t216 + t222 + t227 + t270 + (-t351 * t67 + t369 * t51 + t237 + t240) * MDP(35) + (t351 * t63 - t376 * t51 - t248 + t251) * MDP(36);
t388 = t351 * t1;
t389 = t103 * pkin(5);
t390 = t389 * t41;
t393 = t359 * t1;
t394 = t389 * t38;
t398 = t369 * t41 * t11;
t399 = -t351 * t213;
t403 = t376 * t41 * t11;
t404 = t351 * t208;
t407 = t297 + t298 + (t388 + t299 + t390 - t301) * MDP(28) + (-t393 - t304 + t394 - t305) * MDP(29) + t309 + t312 + t316 + t320 + t323 + (t398 + t325 + t399 + t326) * MDP(35) + (-t403 - t330 + t404 + t331) * MDP(36);
t421 = MDP(13) + t277 + t280 + 0.2e1 * t344 * MDP(21) - 0.2e1 * t338 * MDP(22) + t287 - t290 - t293 + t297 + t298 + 0.2e1 * (t388 + t390) * MDP(28) + 0.2e1 * (-t393 + t394) * MDP(29) + t309 + t312 + t316 + t320 + t323 + 0.2e1 * (t398 + t399) * MDP(35) + 0.2e1 * (-t403 + t404) * MDP(36);
t432 = t38 * MDP(28);
t434 = t41 * MDP(29);
t436 = t63 * t58;
t437 = t41 * MDP(30);
t439 = t58 * t41;
t441 = t63 * t60;
t449 = t60 * t41;
t459 = t321 * t11;
t467 = t29 * t58;
t471 = t15 * MDP(18) - t22 * MDP(19) + t18 * MDP(20) - t30 * MDP(21) + t34 * MDP(22) + t177 * MDP(23) + (t268 + t263) * MDP(24) + t138 * t45 + t128 * t53 - t30 * t432 + t30 * t434 - t436 * t437 + (-t441 * t41 - t439 * t67) * MDP(31) + (t63 * t38 - t439 * t51) * MDP(32) + (t67 * t38 - t449 * t51) * MDP(33) + t51 * t38 * MDP(34) + (t439 * pkin(6) * t51 + t459 * t29 * t60 + t136 * t38) * MDP(35) + (t449 * pkin(6) * t51 - t145 * t38 - t459 * t467) * MDP(36);
t472 = t11 * MDP(18);
t473 = t1 * MDP(19);
t477 = t128 * t41 * MDP(23);
t480 = (-t285 * t11 + t459) * MDP(24);
t481 = t48 * MDP(25);
t482 = t39 * MDP(26);
t485 = t208 * t58;
t486 = t485 * t437;
t488 = t208 * t60;
t491 = (-t439 * t213 - t488 * t41) * MDP(31);
t492 = t58 * t321;
t496 = (-t492 * t11 + t208 * t38) * MDP(32);
t497 = t60 * t321;
t501 = (-t497 * t11 + t213 * t38) * MDP(33);
t503 = t128 * t41 * MDP(34);
t504 = pkin(6) * t11;
t505 = t492 * t504;
t507 = t191 * t60;
t511 = t497 * t504;
t513 = t191 * t58;
t517 = -t472 - t473 - t162 * MDP(21) - t170 * MDP(22) - t477 + t480 + t481 + t482 - t162 * t432 + t162 * t434 - t486 + t491 + t496 + t501 + t503 + (t234 * t38 + t507 * t41 + t505) * MDP(35) + (-t245 * t38 - t513 * t41 + t511) * MDP(36);
t523 = t351 * t60;
t528 = t351 * t58;
t532 = -t472 - t473 + t336 * MDP(21) + t342 * MDP(22) - t477 + t480 + t481 + t482 + t336 * t432 - t336 * t434 - t486 + t491 + t496 + t501 + t503 + (t369 * t38 + t523 * t41 + t505) * MDP(35) + (-t376 * t38 - t528 * t41 + t511) * MDP(36);
t537 = t58 ^ 2;
t550 = pkin(6) * t38;
t579 = t60 * pkin(6);
t585 = t58 * pkin(6);
t591 = t43 * MDP(25) + t51 * MDP(26) + t22 * MDP(27) + t138 * t29 * MDP(28) + t128 * t29 * MDP(29) - t441 * MDP(30) + (-t60 * t67 + t436) * MDP(31) + (-t63 * pkin(7) - t60 * t51) * MDP(32) + (-t67 * pkin(7) + t58 * t51) * MDP(33) - t51 * pkin(7) * MDP(34) + (-t136 * pkin(7) - t138 * t467 + t579 * t51) * MDP(35) + (t145 * pkin(7) - t449 * t34 - t585 * t51) * MDP(36);
t592 = t128 * MDP(25);
t593 = t138 * MDP(26);
t594 = t1 * MDP(27);
t597 = t488 * MDP(30);
t600 = (-t60 * t213 + t485) * MDP(31);
t604 = (-t208 * pkin(7) - t449 * t11) * MDP(32);
t608 = (-t213 * pkin(7) + t439 * t11) * MDP(33);
t609 = pkin(7) * MDP(34);
t610 = t138 * t609;
t611 = t579 * t138;
t615 = t585 * t138;
t619 = -t592 + t593 + t594 + t191 * MDP(28) - t199 * MDP(29) - t597 + t600 + t604 + t608 - t610 + (-t234 * pkin(7) - t513 + t611) * MDP(35) + (t245 * pkin(7) - t507 - t615) * MDP(36);
t628 = -t592 + t593 + t594 + t351 * MDP(28) - t359 * MDP(29) - t597 + t600 + t604 + t608 - t610 + (-t369 * pkin(7) - t528 + t611) * MDP(35) + (t376 * pkin(7) - t523 - t615) * MDP(36);
t633 = t60 ^ 2;
t647 = pkin(6) * pkin(7);
t666 = pkin(7) ^ 2;
t681 = t208 * MDP(32);
t682 = t213 * MDP(33);
t683 = t138 * MDP(34);
unknown(1,1) = t89 + t151;
unknown(2,1) = t217 + t271;
unknown(3,1) = t294 + t334;
unknown(4,1) = t348 + t381;
unknown(5,1) = MDP(13) + t273 - t275 + t277 + t280 + (t344 + t172) * MDP(21) + (-t338 - t166) * MDP(22) + t287 - t290 - t293 + t407;
unknown(6,1) = t421;
unknown(7,1) = t471;
unknown(8,1) = t517;
unknown(9,1) = t532;
unknown(10,1) = 0.2e1 * t41 * t38 * MDP(24) + t537 * t321 * MDP(30) + 0.2e1 * t492 * t60 * MDP(31) - 0.2e1 * t439 * t38 * MDP(32) - 0.2e1 * t449 * t38 * MDP(33) + 0.2e1 * t439 * t550 * MDP(35) + 0.2e1 * t449 * t550 * MDP(36) + t321 * MDP(23) + t285 * MDP(34) + MDP(20);
unknown(11,1) = t591;
unknown(12,1) = t619;
unknown(13,1) = t628;
unknown(14,1) = t41 * MDP(25) + t38 * MDP(26) + t439 * t60 * MDP(30) + (-t537 * t41 + t633 * t41) * MDP(31) + (t439 * pkin(7) - t210) * MDP(32) + (t449 * pkin(7) + t205) * MDP(33) - t38 * pkin(7) * MDP(34) + (t579 * t38 - t439 * t647) * MDP(35) + (-t585 * t38 - t449 * t647) * MDP(36);
unknown(15,1) = -0.2e1 * t60 * t58 * MDP(31) + 0.2e1 * pkin(7) * t60 * MDP(32) - 0.2e1 * pkin(7) * t58 * MDP(33) - 0.2e1 * t579 * pkin(7) * MDP(35) + 0.2e1 * t585 * pkin(7) * MDP(36) + t633 * MDP(30) + t666 * MDP(34) + MDP(27);
unknown(16,1) = t63 * MDP(32) + t67 * MDP(33) + t51 * MDP(34) + t136 * MDP(35) - t145 * MDP(36);
unknown(17,1) = t234 * MDP(35) - t245 * MDP(36) + t681 + t682 + t683;
unknown(18,1) = t369 * MDP(35) - t376 * MDP(36) + t681 + t682 + t683;
unknown(19,1) = t439 * pkin(6) * MDP(35) + t449 * pkin(6) * MDP(36) - t439 * MDP(32) - t449 * MDP(33) + t38 * MDP(34);
unknown(20,1) = -t60 * MDP(32) + t58 * MDP(33) + t579 * MDP(35) - t585 * MDP(36) - t609;
unknown(21,1) = MDP(34);
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [unknown(1), unknown(2), unknown(4), unknown(7), unknown(11), unknown(16); unknown(2), unknown(3), unknown(5), unknown(8), unknown(12), unknown(17); unknown(4), unknown(5), unknown(6), unknown(9), unknown(13), unknown(18); unknown(7), unknown(8), unknown(9), unknown(10), unknown(14), unknown(19); unknown(11), unknown(12), unknown(13), unknown(14), unknown(15), unknown(20); unknown(16), unknown(17), unknown(18), unknown(19), unknown(20), unknown(21);];
Mq = res;
