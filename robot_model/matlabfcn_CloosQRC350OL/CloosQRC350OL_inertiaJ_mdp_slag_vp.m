% Calculate joint inertia matrix for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% MDP [36x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350OL_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(36,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [36 1]), ...
  'CloosQRC350OL_inertiaJ_mdp_slag_vp: MDP has to be [36x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:19:04
% EndTime: 2020-06-20 08:19:10
% DurationCPUTime: 3.10s
% Computational Cost: add. (1246->266), mult. (2738->386), div. (0->0), fcn. (3076->10), ass. (0->134)
t266 = sin(qJ(3));
t267 = sin(qJ(2));
t271 = cos(qJ(3));
t272 = cos(qJ(2));
t226 = -t266 * t272 - t271 * t267;
t227 = -t266 * t267 + t271 * t272;
t264 = sin(qJ(5));
t269 = cos(qJ(5));
t270 = cos(qJ(4));
t336 = t269 * t270;
t194 = t226 * t264 + t227 * t336;
t268 = cos(qJ(6));
t263 = sin(qJ(6));
t265 = sin(qJ(4));
t343 = t263 * t265;
t179 = -t194 * t268 + t227 * t343;
t320 = t179 * MDP(32);
t338 = t265 * t268;
t178 = t194 * t263 + t227 * t338;
t321 = t178 * MDP(33);
t363 = t320 + t321;
t323 = MDP(36) * t263;
t286 = MDP(35) * t268 - t323;
t359 = t286 * pkin(6) - t268 * MDP(32) + t263 * MDP(33);
t362 = MDP(24) + t359;
t329 = MDP(26) * t265;
t330 = MDP(25) * t265;
t361 = t264 * t329 - t269 * t330;
t337 = t265 * t269;
t218 = t263 * t337 - t268 * t270;
t219 = t263 * t270 + t268 * t337;
t342 = t264 * t265;
t239 = MDP(34) * t342;
t307 = t219 * MDP(32) - t218 * MDP(33) + t239;
t360 = 0.2e1 * t227;
t285 = MDP(35) * t263 + MDP(36) * t268;
t276 = -t285 * pkin(6) + t263 * MDP(32) + t268 * MDP(33);
t358 = MDP(26) + t276;
t357 = 2 * pkin(2);
t356 = -0.2e1 * t218;
t355 = 2 * MDP(28);
t354 = 2 * MDP(29);
t353 = 0.2e1 * MDP(35);
t352 = 0.2e1 * MDP(36);
t351 = pkin(6) * t269;
t251 = t267 * pkin(3) + pkin(2);
t250 = pkin(3) * t271 + pkin(4);
t350 = pkin(4) + t250;
t349 = t194 * t270;
t198 = -pkin(4) * t226 + pkin(5) * t227 + t251;
t348 = t198 * t270;
t347 = t250 * t264;
t259 = t265 ^ 2;
t346 = t259 * t264;
t345 = t259 * t269;
t344 = t263 * t264;
t341 = t264 * t268;
t340 = t264 * t269;
t339 = t264 * t270;
t180 = (-pkin(6) * t227 - t198 * t269) * t265;
t293 = pkin(6) * t194 + t348;
t169 = t263 * t180 + t268 * t293;
t311 = t198 * t342;
t335 = t169 * t342 + t218 * t311;
t249 = pkin(3) * t266 - pkin(5);
t204 = t347 + (t249 * t269 - pkin(6)) * t270;
t294 = (t249 - t351) * t265;
t184 = t263 * t204 + t268 * t294;
t208 = -t249 * t339 + t250 * t269;
t334 = t184 * t342 + t208 * t218;
t255 = t264 * pkin(4);
t216 = t255 + (-pkin(5) * t269 - pkin(6)) * t270;
t297 = (-pkin(5) - t351) * t265;
t195 = t263 * t216 + t268 * t297;
t228 = pkin(4) * t269 + pkin(5) * t339;
t333 = t195 * t342 + t228 * t218;
t332 = pkin(5) * t346 + t228 * t270;
t331 = MDP(24) * t265;
t328 = MDP(28) * t270;
t327 = MDP(30) * t179;
t326 = MDP(30) * t263;
t325 = MDP(30) * t268;
t192 = -t269 * t226 + t227 * t339;
t324 = MDP(34) * t192;
t322 = (t180 * t268 - t263 * t293) * MDP(36);
t319 = t194 * MDP(23);
t318 = t194 * MDP(24);
t317 = t226 * MDP(19);
t316 = t226 * MDP(20);
t313 = t270 * MDP(21);
t254 = t270 * MDP(27);
t312 = pkin(6) * t342;
t310 = t227 * t345;
t309 = t249 * t345;
t258 = t264 ^ 2;
t308 = t258 * t338;
t306 = MDP(31) * t263 * t268;
t305 = MDP(32) * t342;
t304 = MDP(33) * t342;
t303 = t227 * t330;
t302 = t227 * t329;
t301 = t265 * t270 * MDP(17);
t299 = t263 * t312;
t298 = t268 * t312;
t296 = pkin(4) * t227 + pkin(5) * t226;
t295 = t208 * t270 - t249 * t346;
t292 = -t226 * t249 + t227 * t250;
t291 = MDP(18) * t270 - MDP(19) * t265;
t290 = -t265 * MDP(22) + t313;
t287 = -t169 * MDP(35) - t322;
t284 = MDP(28) * t269 - MDP(29) * t264 + MDP(21);
t283 = (MDP(14) * t271 - MDP(15) * t266) * pkin(3);
t281 = 0.2e1 * t290;
t280 = -pkin(6) * t258 * t323 - MDP(23) * t340 - MDP(18);
t279 = (-t218 * t263 + t219 * t268) * MDP(31) + t219 * t326 + t263 * t305 + t268 * t304 + t254 + t361;
t261 = t269 ^ 2;
t278 = (-t218 * t269 + t258 * t343) * MDP(33) + (t219 * t269 - t308) * MDP(32) + (t258 - t261) * t331 + t269 * t239 + MDP(25) * t339 + MDP(26) * t336 + (-t219 * t325 + (t218 * t268 + t219 * t263) * MDP(31)) * t264;
t262 = t270 ^ 2;
t277 = (t178 * t219 - t179 * t218) * MDP(31) + t194 * t264 * t331 + t219 * t327 + (-t310 + t349) * MDP(25) + t226 * MDP(12) + t363 * t342 + (-t270 * MDP(26) + t269 * t331 - t307) * t192 + (t346 * MDP(26) + (t259 - t262) * MDP(17) + t265 * t254 + MDP(11)) * t227;
t275 = t262 * MDP(27) + t304 * t356 + MDP(13) + 0.2e1 * t301 + 0.2e1 * t361 * t270 + (MDP(30) * t219 + MDP(31) * t356 + 0.2e1 * t305) * t219 + (MDP(23) * t261 - 0.2e1 * MDP(24) * t340 + MDP(34) * t258 + MDP(16)) * t259;
t274 = -t287 - t324 + t363;
t273 = -t270 * t227 * MDP(16) - t226 * MDP(18) + t264 * t322 - t269 * t319;
t260 = t268 ^ 2;
t257 = t263 ^ 2;
t245 = pkin(5) * t345;
t236 = pkin(6) * t308;
t229 = pkin(5) * t336 - t255;
t209 = t249 * t336 + t347;
t203 = t228 * t219;
t196 = t216 * t268 - t263 * t297;
t191 = t208 * t219;
t185 = t204 * t268 - t263 * t294;
t183 = t219 * t311;
t1 = [t267 * MDP(7) * t357 + t251 * MDP(15) * t360 + MDP(1) + (0.2e1 * t303 + t319) * t194 + (0.2e1 * t178 * MDP(31) + t327) * t179 + (MDP(2) * t272 - 0.2e1 * t267 * MDP(3) + (MDP(8) * t357)) * t272 + (t262 * MDP(16) + t259 * MDP(27) + MDP(9) - 0.2e1 * t301) * t227 ^ 2 + (-0.2e1 * t302 - 0.2e1 * t318 - 0.2e1 * t320 - 0.2e1 * t321 + t324) * t192 + 0.2e1 * t287 * t192 + 0.2e1 * ((t310 + t349) * MDP(29) + t192 * t328 + (t227 * t259 * MDP(28) + (-MDP(35) * t178 + MDP(36) * t179) * t265) * t264) * t198 + (-0.2e1 * t251 * MDP(14) + t316 - 0.2e1 * t290 * t198 + (MDP(10) + t291) * t360) * t226; (t179 * t208 - t185 * t192 + t183) * MDP(36) + (-t178 * t208 - t184 * t192 + t335) * MDP(35) - t267 * MDP(5) + (t292 * MDP(21) + (t192 * t249 + t208 * t227) * MDP(28) + (t194 * t249 - t209 * t227) * MDP(29) + t273) * t265 + t272 * MDP(4) + (t292 * MDP(22) - t317) * t270 + t277; t295 * t355 + (-t209 * t270 - t309) * t354 + t334 * t353 + (t185 * t342 + t191) * t352 + t275 + MDP(6) + 0.2e1 * t283 + t250 * t281; (-t178 * t228 - t192 * t195 + t335) * MDP(35) + (t296 * MDP(22) - t317) * t270 + (t296 * MDP(21) + (-pkin(5) * t192 + t227 * t228) * MDP(28) + (-pkin(5) * t194 + t227 * t229) * MDP(29) + t273) * t265 + (t179 * t228 - t192 * t196 + t183) * MDP(36) + t277; (t333 + t334) * MDP(35) + (t295 + t332) * MDP(28) + t283 + (t191 + t203) * MDP(36) + (t245 - t309) * MDP(29) + (t350 * MDP(21) + (-t209 + t229) * MDP(29)) * t270 + (-t350 * MDP(22) + (t185 + t196) * MDP(36) * t264) * t265 + t275; t332 * t355 + (t229 * t270 + t245) * t354 + t333 * t353 + (t196 * t342 + t203) * t352 + t275 + pkin(4) * t281; t316 + t291 * t227 + (-t313 + (-t285 * t258 + MDP(22)) * t265) * t198 + (-t198 * t328 + t274 + t302 + t318) * t269 + (t319 + t303 + MDP(29) * t348 - t179 * t325 + (-t178 * t268 + t179 * t263) * MDP(31) - t362 * t192) * t264; (t184 * t269 - t208 * t344 + t236) * MDP(35) + (t185 * t269 - t208 * t341) * MDP(36) + (-MDP(22) * t249 - MDP(19)) * t270 + (-t284 * t249 + t280) * t265 + t278; (t195 * t269 - t228 * t344 + t236) * MDP(35) + (t196 * t269 - t228 * t341) * MDP(36) + (pkin(5) * MDP(22) - MDP(19)) * t270 + (t284 * pkin(5) + t280) * t265 + t278; MDP(34) * t261 + MDP(20) + (MDP(30) * t260 + MDP(23) - 0.2e1 * t306) * t258 + 0.2e1 * t362 * t340; t194 * MDP(25) + t179 * t326 + (t178 * t263 + t179 * t268) * MDP(31) + (t227 * MDP(27) + (MDP(29) * t269 + (MDP(28) - t286) * t264) * t198) * t265 - t358 * t192; t208 * MDP(28) - t209 * MDP(29) + (-t208 * t268 - t299) * MDP(35) + (t208 * t263 - t298) * MDP(36) + t279; t228 * MDP(28) + t229 * MDP(29) + (-t228 * t268 - t299) * MDP(35) + (t228 * t263 - t298) * MDP(36) + t279; (MDP(25) - t263 * t325 + (t257 - t260) * MDP(31)) * t264 + t358 * t269; MDP(30) * t257 + MDP(27) + 0.2e1 * t306; t274; MDP(35) * t184 + MDP(36) * t185 + t307; MDP(35) * t195 + MDP(36) * t196 + t307; MDP(34) * t269 + t359 * t264; t276; MDP(34);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
