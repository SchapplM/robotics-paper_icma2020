% Calculate minimal parameter regressor of coriolis joint torque vector for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% tauc_reg [6x36]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = CloosQRC350OL_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:18:58
% EndTime: 2020-06-20 08:19:24
% DurationCPUTime: 16.14s
% Computational Cost: add. (9885->527), mult. (23262->760), div. (0->0), fcn. (18491->10), ass. (0->257)
t176 = sin(qJ(4));
t175 = sin(qJ(5));
t284 = qJD(5) * t175;
t260 = t176 * t284;
t182 = cos(qJ(2));
t339 = cos(qJ(3));
t266 = t339 * t182;
t242 = qJD(1) * t266;
t177 = sin(qJ(3));
t178 = sin(qJ(2));
t291 = qJD(1) * t178;
t265 = t177 * t291;
t133 = -t242 + t265;
t267 = t339 * t178;
t218 = -t177 * t182 - t267;
t134 = t218 * qJD(1);
t180 = cos(qJ(5));
t181 = cos(qJ(4));
t299 = t180 * t181;
t98 = t133 * t175 + t134 * t299;
t360 = t260 - t98;
t350 = t134 + qJD(4);
t357 = qJD(4) * t350;
t363 = t176 * t357;
t170 = qJD(2) + qJD(3);
t224 = t133 * t176 - t170 * t181;
t351 = qJD(5) - t224;
t119 = -t133 * t181 - t176 * t170;
t85 = t119 * t175 - t180 * t350;
t362 = t351 * t85;
t87 = t180 * t119 + t175 * t350;
t361 = t351 * t87;
t359 = t176 * t350;
t317 = t351 * t175;
t358 = t351 * t180;
t174 = sin(qJ(6));
t179 = cos(qJ(6));
t301 = t179 * t180;
t137 = t174 * t181 + t176 * t301;
t303 = t176 * t179;
t337 = t137 * qJD(6) + t350 * t303 + (t299 * qJD(4) - t360) * t174;
t283 = qJD(5) * t179;
t286 = qJD(4) * t180;
t300 = t179 * t181;
t307 = t174 * t176;
t336 = t134 * t307 - t179 * t98 - (qJD(6) + t286) * t300 - (-t175 * t283 + (-qJD(6) * t180 - qJD(4)) * t174) * t176;
t295 = qJD(6) - t85;
t54 = t174 * t87 + t179 * t351;
t356 = t295 * t54;
t55 = -t174 * t351 + t179 * t87;
t355 = t295 * t55;
t113 = t170 * t218;
t104 = t113 * qJD(1);
t102 = t181 * t104;
t66 = qJD(4) * t224 + t102;
t354 = qJD(5) * t350 + t66;
t285 = qJD(4) * t181;
t264 = t175 * t285;
t305 = t175 * t181;
t97 = -t180 * t133 + t134 * t305;
t353 = t97 + t264;
t281 = qJD(5) * t181;
t352 = t175 * t281 + t176 * t286;
t349 = t360 * pkin(6);
t348 = t224 * t350;
t347 = t119 * t350;
t288 = qJD(3) * t177;
t274 = pkin(3) * t288;
t244 = qJD(2) * t274;
t254 = t339 * qJD(2);
t245 = pkin(3) * t254;
t144 = t170 * pkin(4) + t245;
t278 = t144 * qJD(4);
t346 = -t176 * t244 + t181 * t278;
t248 = t178 * t170;
t233 = t177 * t248;
t293 = t170 * t242;
t206 = -qJD(1) * t233 + t293;
t345 = -t176 * t206 + t181 * t357;
t344 = t181 * t206 + t363;
t207 = -t170 * t265 + t293;
t343 = -t181 * t207 - t363;
t289 = qJD(2) * t177;
t275 = pkin(3) * t289;
t282 = qJD(5) * t180;
t107 = -pkin(4) * t133 + t134 * pkin(5);
t92 = -t176 * t107 + t181 * t245;
t342 = -pkin(4) * t282 - t352 * pkin(5) - t175 * t275 + t180 * t92;
t341 = -t176 * t278 - t181 * t244;
t139 = -t177 * t178 + t266;
t287 = qJD(4) * t176;
t297 = t181 * t113;
t216 = -t139 * t287 + t297;
t225 = qJD(3) * t245;
t335 = pkin(3) * qJD(2);
t273 = t182 * t335;
t57 = t104 * pkin(5) + t293 * pkin(4) + (-pkin(4) * t233 + t273) * qJD(1);
t202 = -t176 * t57 + t181 * t225;
t143 = -pkin(5) * t170 + t275;
t149 = qJD(1) * pkin(2) + pkin(3) * t291;
t95 = -pkin(4) * t134 - pkin(5) * t133 + t149;
t74 = t143 * t176 + t181 * t95;
t36 = -qJD(4) * t74 + t202;
t75 = t143 * t181 - t176 * t95;
t23 = -t175 * (qJD(5) * t75 + t244) + t144 * t282 + t180 * t36;
t340 = pkin(6) * t295;
t338 = pkin(6) * t180;
t164 = t178 * pkin(3) + pkin(2);
t211 = -t176 * t104 + t133 * t285 + t170 * t287;
t279 = qJD(6) * t179;
t280 = qJD(6) * t174;
t33 = -t119 * t284 - t175 * t206 + t354 * t180;
t15 = -t174 * t211 - t179 * t33 + t279 * t351 + t87 * t280;
t334 = t15 * t174;
t34 = t119 * t282 + t354 * t175 + t180 * t206;
t333 = t174 * t34;
t332 = t174 * t295;
t331 = t175 * t33;
t210 = t176 * t225 - t95 * t287;
t37 = (qJD(4) * t143 + t57) * t181 + t210;
t330 = t175 * t37;
t329 = t175 * t211;
t328 = t175 * t74;
t327 = t176 * t66;
t326 = t179 * t34;
t325 = t180 * t34;
t324 = t180 * t211;
t323 = t180 * t74;
t322 = t181 * t211;
t321 = t37 * t180;
t320 = t175 * t92 + t180 * t275 - pkin(4) * t284 + (-t175 * t287 + t180 * t281) * pkin(5);
t318 = qJD(4) * t87;
t315 = t224 * t180;
t314 = t119 * t181;
t313 = t133 * t134;
t312 = t133 * t149;
t311 = t134 * t176;
t310 = t139 * t176;
t309 = t144 * t134;
t308 = t149 * t134;
t306 = t175 * t176;
t304 = t176 * t113;
t302 = t176 * t180;
t108 = -pkin(4) * t218 + pkin(5) * t139 + t164;
t298 = t181 * t108;
t184 = qJD(1) ^ 2;
t296 = t182 * t184;
t162 = pkin(3) * t177 - pkin(5);
t163 = t339 * pkin(3) + pkin(4);
t294 = t162 * t299 + t175 * t163;
t292 = t178 ^ 2 - t182 ^ 2;
t290 = qJD(1) * t182;
t276 = qJD(1) * qJD(2);
t253 = t339 * qJD(3);
t114 = (t254 + t253) * t182 + (-t289 - t288) * t178;
t64 = pkin(4) * t114 + pkin(5) * t113 + t273;
t272 = t64 * t306;
t271 = 0.2e1 * t276;
t269 = -pkin(5) - t338;
t268 = t339 * t170;
t261 = t351 * t284;
t258 = t176 * t282;
t59 = t144 * t175 + t180 * t75;
t24 = -qJD(5) * t59 - t175 * t36 - t180 * t244;
t255 = t24 * t181 - t74 * t97;
t252 = t162 - t338;
t251 = t295 ^ 2;
t250 = t179 * t295;
t247 = -0.2e1 * t178 * t276;
t246 = pkin(3) * t253;
t58 = t144 * t180 - t175 * t75;
t243 = t24 * t137 - t336 * t58;
t129 = t252 * t176;
t167 = pkin(6) * t287;
t227 = t181 * t246;
t215 = t163 * t282 - t175 * t274 + t180 * t227;
t99 = pkin(3) * t290 + t107;
t241 = qJD(6) * t129 + (-pkin(6) * t134 - t180 * t99) * t176 - t167 + t352 * t162 - t215;
t239 = -t180 * t285 - t98;
t238 = t97 + t258;
t42 = (-t139 * t281 - t114) * t175 + (qJD(5) * t218 + t216) * t180;
t237 = qJD(6) * t310 - t42;
t120 = -pkin(6) * t181 + t294;
t236 = qJD(6) * t120 + t176 * t246 - t181 * t99 + t252 * t285 + t349;
t141 = t269 * t176;
t235 = -pkin(6) * t311 + qJD(6) * t141 - t167 + t342;
t168 = t175 * pkin(4);
t131 = t168 + (-pkin(5) * t180 - pkin(6)) * t181;
t91 = t181 * t107 + t176 * t245;
t234 = qJD(6) * t131 + t269 * t285 + t349 - t91;
t232 = -pkin(6) * t139 - t108 * t180;
t17 = pkin(6) * t211 + t23;
t18 = t33 * pkin(6) + t37;
t231 = -t179 * t17 + t174 * t18;
t106 = t139 * t299 + t175 * t218;
t61 = pkin(6) * t106 + t298;
t79 = t232 * t176;
t230 = t174 * t79 + t179 * t61;
t229 = t174 * t61 - t179 * t79;
t228 = t360 * t74 + (t287 + t311) * t59;
t222 = -t279 * t295 + t333;
t221 = -t280 * t295 - t326;
t220 = t282 * t351 - t329;
t219 = t261 + t324;
t217 = t139 * t285 + t304;
t214 = pkin(6) * t34;
t212 = pkin(6) * t87 + t74;
t209 = -t282 * t74 - t350 * t58;
t208 = t238 + t264;
t205 = t174 * t212;
t136 = t174 * t302 - t300;
t48 = -pkin(6) * t351 + t59;
t30 = t174 * t48 + t179 * t212;
t6 = -qJD(6) * t205 + t174 * t17 + t179 * t18 + t48 * t279;
t204 = t24 * t136 + t6 * t306 + t337 * t58 + (t258 + t353) * t30;
t16 = t174 * t33 - t351 * t280 + (qJD(6) * t87 - t211) * t179;
t203 = qJD(6) * t106 + t217;
t196 = t75 * t133 - t181 * t309 - t346;
t193 = t181 * t350;
t192 = t350 * t133;
t185 = t74 * t133 - t176 * t309 + t341;
t183 = qJD(2) ^ 2;
t140 = pkin(4) * t180 + pkin(5) * t305;
t123 = -t162 * t305 + t163 * t180;
t105 = t139 * t305 - t180 * t218;
t89 = t133 ^ 2 - t134 ^ 2;
t82 = -t133 * t170 - t206;
t81 = -t134 * t170 + t104;
t80 = (-t162 * t281 - t274) * t180 + (-qJD(5) * t163 + t162 * t287 - t227) * t175;
t78 = -t106 * t179 + t139 * t307;
t77 = t106 * t174 + t139 * t303;
t72 = t119 * t174 - t224 * t301;
t71 = t119 * t179 + t174 * t315;
t51 = pkin(6) * t315 + t75;
t50 = -pkin(6) * t119 - t323;
t43 = qJD(5) * t106 + t180 * t114 + t216 * t175;
t41 = -t119 * t133 - t134 * t193 - t345;
t40 = -t133 * t224 + t134 * t359 + t344;
t39 = -t314 * t350 - t327;
t38 = -t351 * t359 - t322;
t32 = t232 * t285 + (-pkin(6) * t113 + t108 * t284 - t180 * t64) * t176;
t31 = t179 * t48 - t205;
t27 = pkin(6) * t42 - t108 * t287 + t181 * t64;
t26 = -t174 * t237 + t179 * t203;
t25 = t174 * t203 + t179 * t237;
t19 = (-t66 - t348) * t181 + (-t211 + t347) * t176;
t12 = -t33 * t302 + (t239 + t260) * t87;
t11 = t208 * t295 - t34 * t306;
t10 = -t181 * t34 + t353 * t351 + (t350 * t85 + t220) * t176;
t9 = t181 * t33 + t239 * t351 + (-t350 * t87 + t219) * t176;
t8 = t137 * t15 + t336 * t55;
t7 = t85 * t98 + t87 * t97 + (t175 * t87 + t180 * t85) * t285 + (t331 + t325 + (-t175 * t85 + t180 * t87) * qJD(5)) * t176;
t5 = qJD(6) * t30 + t231;
t3 = t136 * t34 + t16 * t306 + t208 * t54 - t295 * t337;
t2 = -t137 * t34 + t15 * t306 - t208 * t55 - t295 * t336;
t1 = -t136 * t15 + t137 * t16 - t336 * t54 + t337 * t55;
t4 = [0, t182 * t247, t292 * t271, -t183 * t178, -t183 * t182, 0, pkin(2) * t182 * t271, pkin(2) * t247, t104 * t139 - t113 * t133, t104 * t218 + t113 * t134 + t133 * t114 - t139 * t206, t113 * t170, -t114 * t170, 0, t149 * t114 - 0.2e1 * t134 * t273 + t164 * t207, t104 * t164 + t113 * t149 + (qJD(1) * t139 - t133) * t273, t139 * t181 * t66 + t119 * t216, (-t119 * t176 + t181 * t224) * t113 + (-t327 + t322 + (-t176 * t224 - t314) * qJD(4)) * t139, t113 * t193 - t119 * t114 - t344 * t139 + t218 * t66, -t113 * t359 - t114 * t224 - t345 * t139 + t211 * t218, -t293 * t218 - qJD(4) * t114 + (t114 * t267 + (t182 * t114 + t218 * t248) * t177) * qJD(1), t344 * t108 + t74 * t114 + t346 * t139 + t144 * t304 - t193 * t64 - t218 * t37, -t108 * t176 * t207 + t75 * t114 + t341 * t139 + t144 * t297 - t218 * t36 + t298 * t357 + t359 * t64, t106 * t33 + t42 * t87, -t105 * t33 - t106 * t34 - t42 * t85 - t43 * t87, t87 * t304 - t106 * t211 + t351 * t42 + (t176 * t33 + t285 * t87) * t139, -t85 * t304 + t105 * t211 - t351 * t43 + (-t176 * t34 - t285 * t85) * t139, -t211 * t310 + t217 * t351, t105 * t37 + t43 * t74 + (t108 * t34 + t64 * t85 + (t108 * t317 + t139 * t58) * qJD(4)) * t181 + (t64 * t317 + t113 * t58 + t139 * t24 + (-qJD(4) * t85 + t220) * t108) * t176, t106 * t37 + t42 * t74 + (t108 * t33 + t64 * t87 + (t108 * t358 - t139 * t59) * qJD(4)) * t181 + (t64 * t358 - t113 * t59 - t139 * t23 + (-t219 - t318) * t108) * t176, t15 * t78 - t25 * t55, t15 * t77 + t16 * t78 + t25 * t54 - t26 * t55, -t105 * t15 + t25 * t295 - t34 * t78 + t43 * t55, -t105 * t16 + t26 * t295 - t34 * t77 - t43 * t54, t105 * t34 - t295 * t43, (-qJD(6) * t229 + t174 * t32 + t179 * t27) * t295 - t230 * t34 - t6 * t105 - t30 * t43 - t54 * t272 - t24 * t77 - t58 * t26 + (-t54 * t258 + (-t16 * t176 - t285 * t54) * t175) * t108, -(qJD(6) * t230 + t174 * t27 - t179 * t32) * t295 + t229 * t34 + t5 * t105 - t31 * t43 - t55 * t272 + t24 * t78 + t58 * t25 + (-t55 * t258 + (t15 * t176 - t285 * t55) * t175) * t108; 0, t178 * t296, -t292 * t184, 0, 0, 0, -pkin(2) * t296, t184 * pkin(2) * t178, t313, t89, t81, t82, 0, t312 + (t134 * t290 + (-qJD(2) - t170) * t288) * pkin(3), -t308 + (t133 * t290 + (-t254 - t268) * qJD(3)) * pkin(3), t39, t19, t41, t40, -t192, -t345 * t162 - t163 * t211 + t193 * t99 + t224 * t274 - t246 * t359 + t185, -t119 * t274 - t162 * t343 + t163 * t66 - t227 * t350 - t359 * t99 + t196, t12, t7, t9, t10, t38, t80 * t351 - t123 * t211 + (-t85 * t99 + (t162 * t85 - t328) * qJD(4)) * t181 + (t85 * t246 + t162 * t34 + (-t351 * t99 - t37) * t175 + t209) * t176 + t255, -t215 * t351 + t294 * t211 + (t162 * t261 - t99 * t87 - t23 + (t162 * t87 - t323) * qJD(4)) * t181 + (t87 * t246 + t162 * t33 + (-t37 + (qJD(4) * t162 - t99) * t351) * t180) * t176 + t228, t8, t1, t2, t3, t11, -(t120 * t174 + t129 * t179) * t34 - t123 * t16 - (t174 * t241 - t179 * t236) * t295 + (t306 * t99 - t80) * t54 + t204, (-t120 * t179 + t129 * t174) * t34 - t80 * t55 + t123 * t15 - (t174 * t236 + t179 * t241) * t295 + t238 * t31 + (t31 * t285 + (t55 * t99 - t5) * t176) * t175 + t243; 0, 0, 0, 0, 0, 0, 0, 0, t313, t89, t81, t82, 0, t312 + (-qJD(3) + t170) * t275, -t308 + (-t253 + t268) * t335, t39, t19, t41, t40, -t192, -pkin(4) * t211 + t345 * pkin(5) - t224 * t275 + t350 * t91 + t185, pkin(4) * t66 + pkin(5) * t343 + t119 * t275 + t350 * t92 + t196, t12, t7, t9, t10, t38, -t140 * t211 - t85 * t91 + t320 * t351 + (-pkin(5) * t85 - t328) * t285 + (-pkin(5) * t34 + t209 - t330) * t176 + t255, t168 * t211 - t91 * t87 + (-pkin(5) * t33 - t321) * t176 + t342 * t351 + (-t74 * t286 - t23 + (-t318 - t324) * pkin(5)) * t181 + t228, t8, t1, t2, t3, t11, -(t131 * t174 + t141 * t179) * t34 - t140 * t16 - (t174 * t235 - t179 * t234) * t295 - t320 * t54 + t204, (-t131 * t179 + t141 * t174) * t34 - t5 * t306 + t140 * t15 - (t174 * t234 + t179 * t235) * t295 - t320 * t55 + t208 * t31 + t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119 * t224, t119 ^ 2 - t224 ^ 2, t133 * t287 - t170 * t285 + t102 - t348, t211 + t347, -t206, -t144 * t119 - t143 * t285 - t181 * t57 + t350 * t75 - t210, t143 * t287 - t144 * t224 + t285 * t95 - t350 * t74 - t202, t358 * t87 + t331, (t33 - t362) * t180 + (-t34 - t361) * t175, -t119 * t87 + t351 * t358 - t329, t119 * t85 - t317 * t351 - t324, -t351 * t119, -t119 * t58 - t75 * t85 - t321, t119 * t59 - t75 * t87 + t330, -t15 * t175 * t179 + (-t175 * t280 + t179 * t282 + t72) * t55, -t54 * t72 + t55 * t71 + (-t174 * t55 - t179 * t54) * t282 + (t334 - t16 * t179 + (t174 * t54 - t179 * t55) * qJD(6)) * t175, -t72 * t295 + (-t283 * t295 + t15) * t180 + (t351 * t55 - t221) * t175, -t71 * t295 + (qJD(5) * t332 + t16) * t180 + (-t351 * t54 - t222) * t175, -t295 * t317 - t325, -(t174 * t50 + t179 * t51) * t295 + t58 * t71 + (t6 + (-t174 * t58 + t179 * t340) * qJD(5)) * t180 + (pkin(6) * t221 - t174 * t24 - t279 * t58 - t30 * t351 + t74 * t54) * t175, (t174 * t51 - t179 * t50) * t295 - t58 * t72 + (-t5 + (-pkin(6) * t332 - t179 * t58) * qJD(5)) * t180 + (pkin(6) * t222 - t24 * t179 + t280 * t58 - t31 * t351 + t74 * t55) * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87 * t85, -t85 ^ 2 + t87 ^ 2, t33 + t362, -t34 + t361, -t211, t351 * t59 - t74 * t87 + t24, t351 * t58 + t74 * t85 - t23, -t250 * t55 + t334, (t15 + t356) * t179 + (t16 + t355) * t174, t250 * t295 - t55 * t87 - t333, -t174 * t251 + t54 * t87 - t326, t295 * t87, t30 * t87 - t59 * t54 + (-pkin(6) * t251 - t24) * t179 + t214 * t174, t31 * t87 - t59 * t55 + (t295 * t340 + t24) * t174 + t214 * t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t54, -t54 ^ 2 + t55 ^ 2, t15 - t356, t16 - t355, -t34, -t295 * t31 + t58 * t55 + t6, -t58 * t54 - t231 + (-qJD(6) + t295) * t30;];
tauc_reg = t4;