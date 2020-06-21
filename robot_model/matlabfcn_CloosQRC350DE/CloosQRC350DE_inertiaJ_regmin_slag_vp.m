% Calculate minimal parameter regressor of joint inertia matrix for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% MM_reg [((6+1)*6/2)x36]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = CloosQRC350DE_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:35
% EndTime: 2020-06-19 21:39:39
% DurationCPUTime: 1.60s
% Computational Cost: add. (1475->183), mult. (3011->387), div. (0->0), fcn. (3280->10), ass. (0->920)
unknown=NaN(21,36);
t1 = cos(qJ(2));
t2 = t1 ^ 2;
t3 = sin(qJ(2));
t10 = cos(qJ(3));
t12 = sin(qJ(3));
t14 = -t1 * t10 + t12 * t3;
t15 = t14 ^ 2;
t18 = t1 * t12 + t10 * t3;
t21 = t3 * pkin(3);
t22 = -pkin(2) - t21;
t27 = cos(qJ(4));
t28 = t27 ^ 2;
t31 = sin(qJ(4));
t34 = t27 * t14;
t37 = t31 * t14;
t40 = t18 ^ 2;
t43 = -pkin(4) * t18 + pkin(5) * t14 - pkin(2) - t21;
t44 = t27 * t43;
t47 = t31 * t43;
t50 = cos(qJ(5));
t51 = t50 * t27;
t53 = sin(qJ(5));
t55 = t14 * t51 + t18 * t53;
t56 = t55 ^ 2;
t57 = t53 * t27;
t60 = -t14 * t57 + t18 * t50;
t69 = t31 ^ 2;
t71 = t53 * t69;
t72 = t43 * t14;
t76 = t50 * t69;
t81 = pkin(7) * qJ(5) - qJ(6);
t82 = cos(t81);
t84 = sin(t81);
t87 = -t14 * t31 * t84 - t55 * t82;
t88 = t87 ^ 2;
t92 = t14 * t31 * t82 - t55 * t84;
t99 = t60 ^ 2;
t100 = t50 * t31;
t101 = t100 * t43;
t103 = -pkin(6) * t37 - t101;
t106 = pkin(6) * t55 + t44;
t108 = -t103 * t84 + t106 * t82;
t110 = t53 * t31;
t116 = -t103 * t82 - t106 * t84;
t121 = t34 * t31;
t124 = -t14 * t28 + t14 * t69;
t125 = t31 * t18;
t126 = t27 * t18;
t127 = t12 * pkin(3);
t128 = t127 - pkin(5);
t129 = t31 * t128;
t131 = t10 * pkin(3);
t132 = t131 + pkin(4);
t133 = t132 * t31;
t136 = t27 * t128;
t138 = t132 * t27;
t141 = t55 * t50;
t142 = t141 * t31;
t144 = t55 * t53;
t146 = -t100 * t60 + t144 * t31;
t149 = -t14 * t76 + t27 * t55;
t152 = t14 * t71 + t27 * t60;
t155 = -t128 * t57 + t132 * t50;
t160 = t51 * t128;
t161 = t53 * t132;
t162 = t160 + t161;
t167 = t82 * t50;
t170 = t167 * t31 - t27 * t84;
t171 = t87 * t170;
t173 = t84 * t50;
t176 = t173 * t31 + t27 * t82;
t178 = t170 * t92 + t176 * t87;
t182 = t31 * t53 * t87 + t170 * t60;
t186 = t31 * t53 * t92 + t176 * t60;
t187 = t60 * t53;
t188 = t187 * t31;
t189 = t27 * pkin(6);
t190 = t160 + t161 - t189;
t192 = t100 * pkin(6);
t193 = t129 - t192;
t195 = -t190 * t84 + t193 * t82;
t198 = t108 * t53 * t31;
t201 = -t110 * t43 * t176;
t205 = -t190 * t82 - t193 * t84;
t208 = t116 * t53 * t31;
t211 = t110 * t43 * t170;
t216 = 0.2e1 * t31 * t27;
t219 = t50 ^ 2;
t220 = t219 * t69;
t222 = 0.2e1 * t76 * t53;
t224 = 0.2e1 * t100 * t27;
t226 = 0.2e1 * t110 * t27;
t227 = t155 * t27;
t228 = t69 * t128;
t229 = t228 * t53;
t231 = t162 * t27;
t232 = t228 * t50;
t234 = t170 ^ 2;
t236 = 0.2e1 * t170 * t176;
t239 = 0.2e1 * t170 * t53 * t31;
t242 = 0.2e1 * t176 * t53 * t31;
t243 = t53 ^ 2;
t244 = t243 * t69;
t246 = t195 * t53 * t31;
t247 = -t155 * t176;
t250 = t205 * t53 * t31;
t251 = t155 * t170;
t253 = t31 * pkin(5);
t255 = pkin(4) * t31;
t258 = t27 * pkin(5);
t260 = pkin(4) * t27;
t265 = pkin(4) * t50 + pkin(5) * t57;
t270 = t51 * pkin(5);
t271 = t53 * pkin(4);
t272 = -t270 + t271;
t277 = -t270 + t271 - t189;
t279 = -t253 - t192;
t281 = -t277 * t84 + t279 * t82;
t287 = -t277 * t82 - t279 * t84;
t293 = t265 * t27;
t294 = t69 * pkin(5);
t295 = t294 * t53;
t297 = t272 * t27;
t298 = t294 * t50;
t301 = t281 * t53 * t31;
t302 = -t265 * t176;
t305 = t287 * t53 * t31;
t306 = t265 * t170;
t319 = t87 * t82;
t321 = t82 * t53;
t323 = t87 * t84;
t329 = t84 * t53;
t337 = t243 * t31;
t344 = t43 * t82;
t347 = t100 * t53;
t349 = -t219 * t31 + t337;
t352 = t170 * t82;
t353 = t352 * t53;
t355 = t170 * t84;
t357 = -t176 * t321 - t355 * t53;
t358 = t82 * t243;
t361 = t170 * t50 - t31 * t358;
t362 = t84 * t243;
t365 = t176 * t50 - t31 * t362;
t366 = pkin(6) * t31;
t367 = t358 * t366;
t369 = t155 * t84;
t372 = t362 * t366;
t374 = t155 * t82;
t380 = t265 * t84;
t384 = t265 * t82;
t389 = t82 ^ 2;
t397 = pkin(6) * t50;
t412 = t84 * pkin(6);
t417 = t82 * pkin(6);
t423 = -t176 * t84 + t352;
t426 = -pkin(7) * t170 - t31 * t329;
t429 = -pkin(7) * t176 + t31 * t321;
t430 = t110 * pkin(7);
t431 = t412 * t110;
t434 = t417 * t110;
t442 = t84 ^ 2;
t452 = pkin(6) * pkin(7);
t464 = pkin(7) ^ 2;
unknown(1,1) = 1;
unknown(1,2) = t2;
unknown(1,3) = -(0.2e1 * t1 * t3);
unknown(1,4) = 0;
unknown(1,5) = 0;
unknown(1,6) = 0;
unknown(1,7) = (0.2e1 * pkin(2) * t3);
unknown(1,8) = (0.2e1 * pkin(2) * t1);
unknown(1,9) = t15;
unknown(1,10) = (0.2e1 * t14 * t18);
unknown(1,11) = 0;
unknown(1,12) = 0;
unknown(1,13) = 0;
unknown(1,14) = -(0.2e1 * t22 * t18);
unknown(1,15) = (0.2e1 * t22 * t14);
unknown(1,16) = (t28 * t15);
unknown(1,17) = -(0.2e1 * t27 * t15 * t31);
unknown(1,18) = (0.2e1 * t34 * t18);
unknown(1,19) = -(0.2e1 * t37 * t18);
unknown(1,20) = t40;
unknown(1,21) = -(0.2e1 * t44 * t18);
unknown(1,22) = (0.2e1 * t47 * t18);
unknown(1,23) = t56;
unknown(1,24) = (0.2e1 * t55 * t60);
unknown(1,25) = (0.2e1 * t55 * t31 * t14);
unknown(1,26) = (0.2e1 * t60 * t31 * t14);
unknown(1,27) = (t69 * t15);
unknown(1,28) = (-0.2e1 * t44 * t60 + 0.2e1 * t71 * t72);
unknown(1,29) = (0.2e1 * t44 * t55 + 0.2e1 * t72 * t76);
unknown(1,30) = t88;
unknown(1,31) = (0.2e1 * t87 * t92);
unknown(1,32) = (0.2e1 * t87 * t60);
unknown(1,33) = (0.2e1 * t92 * t60);
unknown(1,34) = t99;
unknown(1,35) = (-0.2e1 * t110 * t43 * t92 + 0.2e1 * t108 * t60);
unknown(1,36) = (0.2e1 * t110 * t43 * t87 - 0.2e1 * t116 * t60);
unknown(2,1) = 0;
unknown(2,2) = 0;
unknown(2,3) = 0;
unknown(2,4) = -t1;
unknown(2,5) = t3;
unknown(2,6) = 0;
unknown(2,7) = 0;
unknown(2,8) = 0;
unknown(2,9) = 0;
unknown(2,10) = 0;
unknown(2,11) = t14;
unknown(2,12) = t18;
unknown(2,13) = 0;
unknown(2,14) = 0;
unknown(2,15) = 0;
unknown(2,16) = -t121;
unknown(2,17) = t124;
unknown(2,18) = -t125;
unknown(2,19) = -t126;
unknown(2,20) = 0;
unknown(2,21) = (-t129 * t18 + t133 * t14);
unknown(2,22) = (-t136 * t18 + t138 * t14);
unknown(2,23) = -t142;
unknown(2,24) = t146;
unknown(2,25) = t149;
unknown(2,26) = t152;
unknown(2,27) = t121;
unknown(2,28) = (t14 * t155 * t31 - t129 * t60);
unknown(2,29) = (-t14 * t162 * t31 + t129 * t55);
unknown(2,30) = t171;
unknown(2,31) = t178;
unknown(2,32) = t182;
unknown(2,33) = t186;
unknown(2,34) = t188;
unknown(2,35) = (-t155 * t92 + t195 * t60 + t198 + t201);
unknown(2,36) = (t155 * t87 - t205 * t60 - t208 + t211);
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 0;
unknown(3,4) = 0;
unknown(3,5) = 0;
unknown(3,6) = 1;
unknown(3,7) = 0;
unknown(3,8) = 0;
unknown(3,9) = 0;
unknown(3,10) = 0;
unknown(3,11) = 0;
unknown(3,12) = 0;
unknown(3,13) = 1;
unknown(3,14) = (0.2e1 * t131);
unknown(3,15) = -(0.2e1 * t127);
unknown(3,16) = t69;
unknown(3,17) = t216;
unknown(3,18) = 0;
unknown(3,19) = 0;
unknown(3,20) = 0;
unknown(3,21) = (0.2e1 * t138);
unknown(3,22) = -(0.2e1 * t133);
unknown(3,23) = t220;
unknown(3,24) = -t222;
unknown(3,25) = -t224;
unknown(3,26) = t226;
unknown(3,27) = t28;
unknown(3,28) = (0.2e1 * t227 - 0.2e1 * t229);
unknown(3,29) = (-0.2e1 * t231 - 0.2e1 * t232);
unknown(3,30) = t234;
unknown(3,31) = t236;
unknown(3,32) = t239;
unknown(3,33) = t242;
unknown(3,34) = t244;
unknown(3,35) = (0.2e1 * t246 + 0.2e1 * t247);
unknown(3,36) = (-0.2e1 * t250 + 0.2e1 * t251);
unknown(4,1) = 0;
unknown(4,2) = 0;
unknown(4,3) = 0;
unknown(4,4) = 0;
unknown(4,5) = 0;
unknown(4,6) = 0;
unknown(4,7) = 0;
unknown(4,8) = 0;
unknown(4,9) = 0;
unknown(4,10) = 0;
unknown(4,11) = t14;
unknown(4,12) = t18;
unknown(4,13) = 0;
unknown(4,14) = 0;
unknown(4,15) = 0;
unknown(4,16) = -t121;
unknown(4,17) = t124;
unknown(4,18) = -t125;
unknown(4,19) = -t126;
unknown(4,20) = 0;
unknown(4,21) = (t14 * t255 + t18 * t253);
unknown(4,22) = (t14 * t260 + t18 * t258);
unknown(4,23) = -t142;
unknown(4,24) = t146;
unknown(4,25) = t149;
unknown(4,26) = t152;
unknown(4,27) = t121;
unknown(4,28) = (t14 * t265 * t31 + t253 * t60);
unknown(4,29) = (-t14 * t272 * t31 - t253 * t55);
unknown(4,30) = t171;
unknown(4,31) = t178;
unknown(4,32) = t182;
unknown(4,33) = t186;
unknown(4,34) = t188;
unknown(4,35) = (-t265 * t92 + t281 * t60 + t198 + t201);
unknown(4,36) = (t265 * t87 - t287 * t60 - t208 + t211);
unknown(5,1) = 0;
unknown(5,2) = 0;
unknown(5,3) = 0;
unknown(5,4) = 0;
unknown(5,5) = 0;
unknown(5,6) = 0;
unknown(5,7) = 0;
unknown(5,8) = 0;
unknown(5,9) = 0;
unknown(5,10) = 0;
unknown(5,11) = 0;
unknown(5,12) = 0;
unknown(5,13) = 1;
unknown(5,14) = t131;
unknown(5,15) = -t127;
unknown(5,16) = t69;
unknown(5,17) = t216;
unknown(5,18) = 0;
unknown(5,19) = 0;
unknown(5,20) = 0;
unknown(5,21) = (t260 + t138);
unknown(5,22) = (-t255 - t133);
unknown(5,23) = t220;
unknown(5,24) = -t222;
unknown(5,25) = -t224;
unknown(5,26) = t226;
unknown(5,27) = t28;
unknown(5,28) = (t293 + t227 + t295 - t229);
unknown(5,29) = (-t297 - t231 + t298 - t232);
unknown(5,30) = t234;
unknown(5,31) = t236;
unknown(5,32) = t239;
unknown(5,33) = t242;
unknown(5,34) = t244;
unknown(5,35) = (t301 + t246 + t302 + t247);
unknown(5,36) = (-t305 - t250 + t306 + t251);
unknown(6,1) = 0;
unknown(6,2) = 0;
unknown(6,3) = 0;
unknown(6,4) = 0;
unknown(6,5) = 0;
unknown(6,6) = 0;
unknown(6,7) = 0;
unknown(6,8) = 0;
unknown(6,9) = 0;
unknown(6,10) = 0;
unknown(6,11) = 0;
unknown(6,12) = 0;
unknown(6,13) = 1;
unknown(6,14) = 0;
unknown(6,15) = 0;
unknown(6,16) = t69;
unknown(6,17) = t216;
unknown(6,18) = 0;
unknown(6,19) = 0;
unknown(6,20) = 0;
unknown(6,21) = (0.2e1 * t260);
unknown(6,22) = -(0.2e1 * t255);
unknown(6,23) = t220;
unknown(6,24) = -t222;
unknown(6,25) = -t224;
unknown(6,26) = t226;
unknown(6,27) = t28;
unknown(6,28) = (0.2e1 * t293 + 0.2e1 * t295);
unknown(6,29) = (-0.2e1 * t297 + 0.2e1 * t298);
unknown(6,30) = t234;
unknown(6,31) = t236;
unknown(6,32) = t239;
unknown(6,33) = t242;
unknown(6,34) = t244;
unknown(6,35) = (0.2e1 * t301 + 0.2e1 * t302);
unknown(6,36) = (-0.2e1 * t305 + 0.2e1 * t306);
unknown(7,1) = 0;
unknown(7,2) = 0;
unknown(7,3) = 0;
unknown(7,4) = 0;
unknown(7,5) = 0;
unknown(7,6) = 0;
unknown(7,7) = 0;
unknown(7,8) = 0;
unknown(7,9) = 0;
unknown(7,10) = 0;
unknown(7,11) = 0;
unknown(7,12) = 0;
unknown(7,13) = 0;
unknown(7,14) = 0;
unknown(7,15) = 0;
unknown(7,16) = 0;
unknown(7,17) = 0;
unknown(7,18) = t34;
unknown(7,19) = -t37;
unknown(7,20) = t18;
unknown(7,21) = -t44;
unknown(7,22) = t47;
unknown(7,23) = t144;
unknown(7,24) = (t187 + t141);
unknown(7,25) = (t110 * t14);
unknown(7,26) = (t100 * t14);
unknown(7,27) = 0;
unknown(7,28) = -(t44 * t50);
unknown(7,29) = (t44 * t53);
unknown(7,30) = -(t319 * t53);
unknown(7,31) = (-t321 * t92 - t323 * t53);
unknown(7,32) = (-t321 * t60 + t50 * t87);
unknown(7,33) = (-t329 * t60 + t50 * t92);
unknown(7,34) = (t60 * t50);
unknown(7,35) = (pkin(6) * t321 * t60 + t337 * t43 * t84 + t108 * t50);
unknown(7,36) = (pkin(6) * t329 * t60 - t116 * t50 - t337 * t344);
unknown(8,1) = 0;
unknown(8,2) = 0;
unknown(8,3) = 0;
unknown(8,4) = 0;
unknown(8,5) = 0;
unknown(8,6) = 0;
unknown(8,7) = 0;
unknown(8,8) = 0;
unknown(8,9) = 0;
unknown(8,10) = 0;
unknown(8,11) = 0;
unknown(8,12) = 0;
unknown(8,13) = 0;
unknown(8,14) = 0;
unknown(8,15) = 0;
unknown(8,16) = 0;
unknown(8,17) = 0;
unknown(8,18) = -t31;
unknown(8,19) = -t27;
unknown(8,20) = 0;
unknown(8,21) = -t129;
unknown(8,22) = -t136;
unknown(8,23) = -t347;
unknown(8,24) = t349;
unknown(8,25) = t57;
unknown(8,26) = t51;
unknown(8,27) = 0;
unknown(8,28) = -(t129 * t50);
unknown(8,29) = (t129 * t53);
unknown(8,30) = -t353;
unknown(8,31) = t357;
unknown(8,32) = t361;
unknown(8,33) = t365;
unknown(8,34) = t347;
unknown(8,35) = (t195 * t50 + t369 * t53 + t367);
unknown(8,36) = (-t205 * t50 - t374 * t53 + t372);
unknown(9,1) = 0;
unknown(9,2) = 0;
unknown(9,3) = 0;
unknown(9,4) = 0;
unknown(9,5) = 0;
unknown(9,6) = 0;
unknown(9,7) = 0;
unknown(9,8) = 0;
unknown(9,9) = 0;
unknown(9,10) = 0;
unknown(9,11) = 0;
unknown(9,12) = 0;
unknown(9,13) = 0;
unknown(9,14) = 0;
unknown(9,15) = 0;
unknown(9,16) = 0;
unknown(9,17) = 0;
unknown(9,18) = -t31;
unknown(9,19) = -t27;
unknown(9,20) = 0;
unknown(9,21) = t253;
unknown(9,22) = t258;
unknown(9,23) = -t347;
unknown(9,24) = t349;
unknown(9,25) = t57;
unknown(9,26) = t51;
unknown(9,27) = 0;
unknown(9,28) = (t253 * t50);
unknown(9,29) = -(t253 * t53);
unknown(9,30) = -t353;
unknown(9,31) = t357;
unknown(9,32) = t361;
unknown(9,33) = t365;
unknown(9,34) = t347;
unknown(9,35) = (t281 * t50 + t380 * t53 + t367);
unknown(9,36) = (-t287 * t50 - t384 * t53 + t372);
unknown(10,1) = 0;
unknown(10,2) = 0;
unknown(10,3) = 0;
unknown(10,4) = 0;
unknown(10,5) = 0;
unknown(10,6) = 0;
unknown(10,7) = 0;
unknown(10,8) = 0;
unknown(10,9) = 0;
unknown(10,10) = 0;
unknown(10,11) = 0;
unknown(10,12) = 0;
unknown(10,13) = 0;
unknown(10,14) = 0;
unknown(10,15) = 0;
unknown(10,16) = 0;
unknown(10,17) = 0;
unknown(10,18) = 0;
unknown(10,19) = 0;
unknown(10,20) = 1;
unknown(10,21) = 0;
unknown(10,22) = 0;
unknown(10,23) = t243;
unknown(10,24) = (0.2e1 * t53 * t50);
unknown(10,25) = 0;
unknown(10,26) = 0;
unknown(10,27) = 0;
unknown(10,28) = 0;
unknown(10,29) = 0;
unknown(10,30) = (t389 * t243);
unknown(10,31) = (0.2e1 * t358 * t84);
unknown(10,32) = -(0.2e1 * t321 * t50);
unknown(10,33) = -(0.2e1 * t329 * t50);
unknown(10,34) = t219;
unknown(10,35) = (0.2e1 * t321 * t397);
unknown(10,36) = (0.2e1 * t329 * t397);
unknown(11,1) = 0;
unknown(11,2) = 0;
unknown(11,3) = 0;
unknown(11,4) = 0;
unknown(11,5) = 0;
unknown(11,6) = 0;
unknown(11,7) = 0;
unknown(11,8) = 0;
unknown(11,9) = 0;
unknown(11,10) = 0;
unknown(11,11) = 0;
unknown(11,12) = 0;
unknown(11,13) = 0;
unknown(11,14) = 0;
unknown(11,15) = 0;
unknown(11,16) = 0;
unknown(11,17) = 0;
unknown(11,18) = 0;
unknown(11,19) = 0;
unknown(11,20) = 0;
unknown(11,21) = 0;
unknown(11,22) = 0;
unknown(11,23) = 0;
unknown(11,24) = 0;
unknown(11,25) = t55;
unknown(11,26) = t60;
unknown(11,27) = t37;
unknown(11,28) = (t110 * t43);
unknown(11,29) = t101;
unknown(11,30) = -t323;
unknown(11,31) = (-t84 * t92 + t319);
unknown(11,32) = (-pkin(7) * t87 - t60 * t84);
unknown(11,33) = (-pkin(7) * t92 + t60 * t82);
unknown(11,34) = -(t60 * pkin(7));
unknown(11,35) = (-pkin(7) * t108 - t110 * t344 + t412 * t60);
unknown(11,36) = (pkin(7) * t116 - t329 * t47 - t417 * t60);
unknown(12,1) = 0;
unknown(12,2) = 0;
unknown(12,3) = 0;
unknown(12,4) = 0;
unknown(12,5) = 0;
unknown(12,6) = 0;
unknown(12,7) = 0;
unknown(12,8) = 0;
unknown(12,9) = 0;
unknown(12,10) = 0;
unknown(12,11) = 0;
unknown(12,12) = 0;
unknown(12,13) = 0;
unknown(12,14) = 0;
unknown(12,15) = 0;
unknown(12,16) = 0;
unknown(12,17) = 0;
unknown(12,18) = 0;
unknown(12,19) = 0;
unknown(12,20) = 0;
unknown(12,21) = 0;
unknown(12,22) = 0;
unknown(12,23) = 0;
unknown(12,24) = 0;
unknown(12,25) = -t100;
unknown(12,26) = t110;
unknown(12,27) = t27;
unknown(12,28) = t155;
unknown(12,29) = -t162;
unknown(12,30) = -t355;
unknown(12,31) = t423;
unknown(12,32) = t426;
unknown(12,33) = t429;
unknown(12,34) = -t430;
unknown(12,35) = (-pkin(7) * t195 - t374 + t431);
unknown(12,36) = (pkin(7) * t205 - t369 - t434);
unknown(13,1) = 0;
unknown(13,2) = 0;
unknown(13,3) = 0;
unknown(13,4) = 0;
unknown(13,5) = 0;
unknown(13,6) = 0;
unknown(13,7) = 0;
unknown(13,8) = 0;
unknown(13,9) = 0;
unknown(13,10) = 0;
unknown(13,11) = 0;
unknown(13,12) = 0;
unknown(13,13) = 0;
unknown(13,14) = 0;
unknown(13,15) = 0;
unknown(13,16) = 0;
unknown(13,17) = 0;
unknown(13,18) = 0;
unknown(13,19) = 0;
unknown(13,20) = 0;
unknown(13,21) = 0;
unknown(13,22) = 0;
unknown(13,23) = 0;
unknown(13,24) = 0;
unknown(13,25) = -t100;
unknown(13,26) = t110;
unknown(13,27) = t27;
unknown(13,28) = t265;
unknown(13,29) = -t272;
unknown(13,30) = -t355;
unknown(13,31) = t423;
unknown(13,32) = t426;
unknown(13,33) = t429;
unknown(13,34) = -t430;
unknown(13,35) = (-pkin(7) * t281 - t384 + t431);
unknown(13,36) = (pkin(7) * t287 - t380 - t434);
unknown(14,1) = 0;
unknown(14,2) = 0;
unknown(14,3) = 0;
unknown(14,4) = 0;
unknown(14,5) = 0;
unknown(14,6) = 0;
unknown(14,7) = 0;
unknown(14,8) = 0;
unknown(14,9) = 0;
unknown(14,10) = 0;
unknown(14,11) = 0;
unknown(14,12) = 0;
unknown(14,13) = 0;
unknown(14,14) = 0;
unknown(14,15) = 0;
unknown(14,16) = 0;
unknown(14,17) = 0;
unknown(14,18) = 0;
unknown(14,19) = 0;
unknown(14,20) = 0;
unknown(14,21) = 0;
unknown(14,22) = 0;
unknown(14,23) = 0;
unknown(14,24) = 0;
unknown(14,25) = t53;
unknown(14,26) = t50;
unknown(14,27) = 0;
unknown(14,28) = 0;
unknown(14,29) = 0;
unknown(14,30) = (t321 * t84);
unknown(14,31) = (-t389 * t53 + t442 * t53);
unknown(14,32) = (pkin(7) * t321 - t173);
unknown(14,33) = (pkin(7) * t329 + t167);
unknown(14,34) = -(t50 * pkin(7));
unknown(14,35) = (-t321 * t452 + t412 * t50);
unknown(14,36) = (-t329 * t452 - t417 * t50);
unknown(15,1) = 0;
unknown(15,2) = 0;
unknown(15,3) = 0;
unknown(15,4) = 0;
unknown(15,5) = 0;
unknown(15,6) = 0;
unknown(15,7) = 0;
unknown(15,8) = 0;
unknown(15,9) = 0;
unknown(15,10) = 0;
unknown(15,11) = 0;
unknown(15,12) = 0;
unknown(15,13) = 0;
unknown(15,14) = 0;
unknown(15,15) = 0;
unknown(15,16) = 0;
unknown(15,17) = 0;
unknown(15,18) = 0;
unknown(15,19) = 0;
unknown(15,20) = 0;
unknown(15,21) = 0;
unknown(15,22) = 0;
unknown(15,23) = 0;
unknown(15,24) = 0;
unknown(15,25) = 0;
unknown(15,26) = 0;
unknown(15,27) = 1;
unknown(15,28) = 0;
unknown(15,29) = 0;
unknown(15,30) = t442;
unknown(15,31) = -(0.2e1 * t84 * t82);
unknown(15,32) = (0.2e1 * pkin(7) * t84);
unknown(15,33) = -(0.2e1 * pkin(7) * t82);
unknown(15,34) = t464;
unknown(15,35) = -(0.2e1 * t412 * pkin(7));
unknown(15,36) = (0.2e1 * t417 * pkin(7));
unknown(16,1) = 0;
unknown(16,2) = 0;
unknown(16,3) = 0;
unknown(16,4) = 0;
unknown(16,5) = 0;
unknown(16,6) = 0;
unknown(16,7) = 0;
unknown(16,8) = 0;
unknown(16,9) = 0;
unknown(16,10) = 0;
unknown(16,11) = 0;
unknown(16,12) = 0;
unknown(16,13) = 0;
unknown(16,14) = 0;
unknown(16,15) = 0;
unknown(16,16) = 0;
unknown(16,17) = 0;
unknown(16,18) = 0;
unknown(16,19) = 0;
unknown(16,20) = 0;
unknown(16,21) = 0;
unknown(16,22) = 0;
unknown(16,23) = 0;
unknown(16,24) = 0;
unknown(16,25) = 0;
unknown(16,26) = 0;
unknown(16,27) = 0;
unknown(16,28) = 0;
unknown(16,29) = 0;
unknown(16,30) = 0;
unknown(16,31) = 0;
unknown(16,32) = t87;
unknown(16,33) = t92;
unknown(16,34) = t60;
unknown(16,35) = t108;
unknown(16,36) = -t116;
unknown(17,1) = 0;
unknown(17,2) = 0;
unknown(17,3) = 0;
unknown(17,4) = 0;
unknown(17,5) = 0;
unknown(17,6) = 0;
unknown(17,7) = 0;
unknown(17,8) = 0;
unknown(17,9) = 0;
unknown(17,10) = 0;
unknown(17,11) = 0;
unknown(17,12) = 0;
unknown(17,13) = 0;
unknown(17,14) = 0;
unknown(17,15) = 0;
unknown(17,16) = 0;
unknown(17,17) = 0;
unknown(17,18) = 0;
unknown(17,19) = 0;
unknown(17,20) = 0;
unknown(17,21) = 0;
unknown(17,22) = 0;
unknown(17,23) = 0;
unknown(17,24) = 0;
unknown(17,25) = 0;
unknown(17,26) = 0;
unknown(17,27) = 0;
unknown(17,28) = 0;
unknown(17,29) = 0;
unknown(17,30) = 0;
unknown(17,31) = 0;
unknown(17,32) = t170;
unknown(17,33) = t176;
unknown(17,34) = t110;
unknown(17,35) = t195;
unknown(17,36) = -t205;
unknown(18,1) = 0;
unknown(18,2) = 0;
unknown(18,3) = 0;
unknown(18,4) = 0;
unknown(18,5) = 0;
unknown(18,6) = 0;
unknown(18,7) = 0;
unknown(18,8) = 0;
unknown(18,9) = 0;
unknown(18,10) = 0;
unknown(18,11) = 0;
unknown(18,12) = 0;
unknown(18,13) = 0;
unknown(18,14) = 0;
unknown(18,15) = 0;
unknown(18,16) = 0;
unknown(18,17) = 0;
unknown(18,18) = 0;
unknown(18,19) = 0;
unknown(18,20) = 0;
unknown(18,21) = 0;
unknown(18,22) = 0;
unknown(18,23) = 0;
unknown(18,24) = 0;
unknown(18,25) = 0;
unknown(18,26) = 0;
unknown(18,27) = 0;
unknown(18,28) = 0;
unknown(18,29) = 0;
unknown(18,30) = 0;
unknown(18,31) = 0;
unknown(18,32) = t170;
unknown(18,33) = t176;
unknown(18,34) = t110;
unknown(18,35) = t281;
unknown(18,36) = -t287;
unknown(19,1) = 0;
unknown(19,2) = 0;
unknown(19,3) = 0;
unknown(19,4) = 0;
unknown(19,5) = 0;
unknown(19,6) = 0;
unknown(19,7) = 0;
unknown(19,8) = 0;
unknown(19,9) = 0;
unknown(19,10) = 0;
unknown(19,11) = 0;
unknown(19,12) = 0;
unknown(19,13) = 0;
unknown(19,14) = 0;
unknown(19,15) = 0;
unknown(19,16) = 0;
unknown(19,17) = 0;
unknown(19,18) = 0;
unknown(19,19) = 0;
unknown(19,20) = 0;
unknown(19,21) = 0;
unknown(19,22) = 0;
unknown(19,23) = 0;
unknown(19,24) = 0;
unknown(19,25) = 0;
unknown(19,26) = 0;
unknown(19,27) = 0;
unknown(19,28) = 0;
unknown(19,29) = 0;
unknown(19,30) = 0;
unknown(19,31) = 0;
unknown(19,32) = -t321;
unknown(19,33) = -t329;
unknown(19,34) = t50;
unknown(19,35) = (t321 * pkin(6));
unknown(19,36) = (t329 * pkin(6));
unknown(20,1) = 0;
unknown(20,2) = 0;
unknown(20,3) = 0;
unknown(20,4) = 0;
unknown(20,5) = 0;
unknown(20,6) = 0;
unknown(20,7) = 0;
unknown(20,8) = 0;
unknown(20,9) = 0;
unknown(20,10) = 0;
unknown(20,11) = 0;
unknown(20,12) = 0;
unknown(20,13) = 0;
unknown(20,14) = 0;
unknown(20,15) = 0;
unknown(20,16) = 0;
unknown(20,17) = 0;
unknown(20,18) = 0;
unknown(20,19) = 0;
unknown(20,20) = 0;
unknown(20,21) = 0;
unknown(20,22) = 0;
unknown(20,23) = 0;
unknown(20,24) = 0;
unknown(20,25) = 0;
unknown(20,26) = 0;
unknown(20,27) = 0;
unknown(20,28) = 0;
unknown(20,29) = 0;
unknown(20,30) = 0;
unknown(20,31) = 0;
unknown(20,32) = -t84;
unknown(20,33) = t82;
unknown(20,34) = -pkin(7);
unknown(20,35) = t412;
unknown(20,36) = -t417;
unknown(21,1) = 0;
unknown(21,2) = 0;
unknown(21,3) = 0;
unknown(21,4) = 0;
unknown(21,5) = 0;
unknown(21,6) = 0;
unknown(21,7) = 0;
unknown(21,8) = 0;
unknown(21,9) = 0;
unknown(21,10) = 0;
unknown(21,11) = 0;
unknown(21,12) = 0;
unknown(21,13) = 0;
unknown(21,14) = 0;
unknown(21,15) = 0;
unknown(21,16) = 0;
unknown(21,17) = 0;
unknown(21,18) = 0;
unknown(21,19) = 0;
unknown(21,20) = 0;
unknown(21,21) = 0;
unknown(21,22) = 0;
unknown(21,23) = 0;
unknown(21,24) = 0;
unknown(21,25) = 0;
unknown(21,26) = 0;
unknown(21,27) = 0;
unknown(21,28) = 0;
unknown(21,29) = 0;
unknown(21,30) = 0;
unknown(21,31) = 0;
unknown(21,32) = 0;
unknown(21,33) = 0;
unknown(21,34) = 1;
unknown(21,35) = 0;
unknown(21,36) = 0;
MM_reg = unknown;
