% Calculate vector of inverse dynamics joint torques for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = CloosQRC350OL_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:23
% EndTime: 2020-06-23 22:03:28
% DurationCPUTime: 3.24s
% Computational Cost: add. (2141->239), mult. (4742->366), div. (0->0), fcn. (3785->12), ass. (0->116)
t236 = sin(qJ(3));
t241 = cos(qJ(3));
t242 = cos(qJ(2));
t285 = qJD(1) * t242;
t237 = sin(qJ(2));
t286 = qJD(1) * t237;
t201 = t236 * t286 - t241 * t285;
t204 = t236 * t242 + t237 * t241;
t202 = t204 * qJD(1);
t234 = sin(qJ(5));
t239 = cos(qJ(5));
t240 = cos(qJ(4));
t288 = t239 * t240;
t184 = t201 * t234 - t202 * t288;
t235 = sin(qJ(4));
t277 = qJD(5) * t235;
t301 = t234 * t277 - t184;
t273 = qJDD(1) * t237;
t274 = qJD(1) * qJD(2);
t300 = t242 * t274 + t273;
t229 = qJD(2) + qJD(3);
t193 = -t201 * t240 - t229 * t235;
t199 = qJD(4) - t202;
t180 = -t193 * t234 + t239 * t199;
t299 = pkin(3) * t236;
t298 = pkin(3) * t241;
t230 = qJ(2) + qJ(3);
t226 = cos(t230);
t217 = g(3) * t226;
t221 = t237 * pkin(3) + pkin(2);
t297 = pkin(3) * qJD(2);
t189 = t229 * t204;
t272 = qJDD(1) * t242;
t178 = -t189 * qJD(1) - t236 * t273 + t241 * t272;
t263 = t237 * t274;
t283 = qJD(3) * t236;
t265 = t237 * t283;
t251 = -qJD(1) * t265 + (-t263 + t272) * t236;
t257 = t229 * t242;
t179 = (qJD(1) * t257 + t273) * t241 + t251;
t200 = qJDD(1) * pkin(2) + t300 * pkin(3);
t163 = pkin(4) * t179 + pkin(5) * t178 + t200;
t208 = qJD(1) * pkin(2) + pkin(3) * t286;
t183 = pkin(4) * t202 - pkin(5) * t201 + t208;
t284 = qJD(2) * t236;
t270 = pkin(3) * t284;
t206 = -pkin(5) * t229 + t270;
t175 = -t183 * t235 + t206 * t240;
t228 = qJDD(2) + qJDD(3);
t271 = qJDD(2) * t236;
t282 = qJD(3) * t241;
t197 = -pkin(5) * t228 + (qJD(2) * t282 + t271) * pkin(3);
t157 = t175 * qJD(4) + t163 * t240 + t197 * t235;
t296 = t157 * t235;
t192 = t201 * t235 - t229 * t240;
t191 = qJD(5) - t192;
t295 = t191 * t239;
t205 = -t236 * t237 + t241 * t242;
t293 = t205 * t235;
t292 = t205 * t240;
t291 = t234 * t240;
t166 = t192 * qJD(4) + t178 * t240 - t228 * t235;
t244 = (-t229 * t285 - t273) * t241 - t251;
t176 = qJDD(4) + t244;
t159 = t180 * qJD(5) + t166 * t239 + t176 * t234;
t290 = t235 * t159;
t238 = cos(qJ(6));
t289 = t238 * t239;
t287 = t242 * qJD(1) ^ 2;
t281 = qJD(4) * t235;
t280 = qJD(4) * t239;
t279 = qJD(4) * t240;
t278 = qJD(5) * t234;
t276 = qJD(5) * t240;
t174 = t183 * t240 + t206 * t235;
t275 = t174 * qJD(4);
t177 = qJD(6) + t180;
t269 = t241 * t297;
t268 = t242 * t297;
t267 = pkin(3) * t283;
t266 = pkin(3) * t282;
t187 = -pkin(4) * t201 - t202 * pkin(5);
t261 = t202 * t208 + t217;
t219 = -pkin(5) + t299;
t260 = -pkin(3) * t285 + qJD(4) * t219 - t187;
t259 = qJD(2) * (-qJD(3) + t229);
t258 = qJD(3) * (-qJD(2) - t229);
t223 = qJDD(2) * t298;
t225 = sin(t230);
t256 = g(3) * t225 + t201 * t208 + t223;
t190 = -t237 * t284 + t241 * t257 - t265;
t255 = -t205 * t276 - t190;
t254 = t239 * (-t163 * t235 + t197 * t240 - t275) + t234 * (pkin(4) * t228 - qJD(2) * t267 + t223);
t207 = pkin(4) * t229 + t269;
t253 = t175 * t234 - t207 * t239;
t181 = t193 * t239 + t199 * t234;
t233 = sin(qJ(6));
t252 = t181 * t233 + t191 * t238;
t250 = -t189 * t235 + t205 * t279;
t249 = -t189 * t240 - t205 * t281;
t165 = t193 * qJD(4) + t178 * t235 + t228 * t240 + qJDD(5);
t248 = -qJD(4) * t181 + t165 * t239 - t191 * t278;
t247 = qJD(5) * t204 - t249;
t155 = -t253 * qJD(5) + t254;
t169 = t175 * t239 + t207 * t234;
t246 = -g(3) * (t225 * t291 - t226 * t239) - t155 * t240 + t301 * t174 + (-t235 * t202 + t281) * t169;
t154 = t252 * qJD(6) - t159 * t238 + t165 * t233;
t158 = -t181 * qJD(5) - t166 * t234 + t176 * t239 + qJDD(6);
t167 = -t181 * t238 + t191 * t233;
t245 = (t154 * (t233 * t240 + t235 * t289) + ((t184 + (qJD(6) + t280) * t240) * t238 + (-t238 * t278 + (-qJD(6) * t239 - t199) * t233) * t235) * t167) * MDP(18) + (t158 * t234 * t235 + ((-t201 + t277) * t239 + t199 * t291) * t177) * MDP(19) + (-t239 * t290 + (-t239 * t279 + t301) * t181) * MDP(15) + (-t199 * t235 * t191 + t240 * t165) * MDP(16) + (-t199 * t240 * t193 - t166 * t235) * MDP(13) + (t202 * t229 + t178) * MDP(8) + (-t201 * t229 + t244) * MDP(9) + (t201 ^ 2 - t202 ^ 2) * MDP(7) + t228 * MDP(10) + (-t199 * MDP(14) - MDP(6) * t202) * t201;
t220 = pkin(4) + t298;
t188 = pkin(4) * t204 + pkin(5) * t205 + t221;
t186 = -t204 * t234 + t205 * t288;
t172 = pkin(4) * t190 - pkin(5) * t189 + t268;
t162 = t255 * t234 - t247 * t239;
t1 = [qJDD(1) * MDP(1) + (-0.2e1 * t263 + t272) * t242 * MDP(2) + (-qJD(2) ^ 2 * t237 + t242 * qJDD(2)) * MDP(3) + 0.2e1 * t300 * MDP(5) * pkin(2) + (t178 * t205 + t189 * t201) * MDP(6) + (-t178 * t204 - t179 * t205 + t189 * t202 + t190 * t201) * MDP(7) + (-t189 * t229 + t205 * t228) * MDP(8) + (-t190 * t229 - t204 * t228) * MDP(9) + (t179 * t221 + t190 * t208 + t200 * t204 + t202 * t268) * MDP(11) + (t178 * t221 - t189 * t208 + t200 * t205 - t201 * t268) * MDP(12) + (t166 * t292 + t249 * t193) * MDP(13) + (-t176 * t204 - t190 * t199) * MDP(14) + (t159 * t186 + t162 * t181) * MDP(15) + (t165 * t293 + t250 * t191) * MDP(16) + (t157 * t186 + t162 * t174 + (t159 * t188 + t172 * t181 + (-t169 * t205 + t188 * t295) * qJD(4)) * t240 + (-t155 * t205 + t169 * t189 + t172 * t295 + t248 * t188) * t235) * MDP(17) + ((-t154 * t186 + t167 * (qJD(6) * t293 - t162)) * t238 + (t154 * t293 + t167 * (qJD(6) * t186 + t250)) * t233) * MDP(18) + ((-t158 * t204 + t177 * t255) * t239 + (-t158 * t292 + t177 * t247) * t234) * MDP(19); MDP(3) * t272 + t245 + qJDD(2) * MDP(4) + (t219 * t290 + (-(-t219 * t276 - t267) * t191 - t220 * t165) * t234 + (t235 * t266 + t260 * t240) * t181 + (-t296 + (-t219 * t165 - t275) * t240 + (-qJD(5) * t220 + t260 * t235 - t240 * t266) * t191) * t239 + t246) * MDP(17) + ((-t202 * t285 + t228 * t241 + t236 * t258) * pkin(3) + t256) * MDP(11) + ((t201 * t285 + (-qJDD(2) - t228) * t236 + t241 * t258) * pkin(3) + t261) * MDP(12) + t237 * MDP(2) * t287 + (-pkin(2) * t287 + g(3) * t237) * MDP(5); t245 + (t259 * t299 + t256) * MDP(11) + (-(t187 * t240 + t235 * t269) * t181 + (-pkin(4) * t165 - t191 * t270) * t234 + (-t240 * t275 - t296 + (-pkin(4) * qJD(5) - t187 * t235 + t240 * t269) * t191) * t239 + ((-t191 * t280 - t159) * t235 + t248 * t240) * pkin(5) + t246) * MDP(17) + ((t241 * t259 - t271) * pkin(3) + t261) * MDP(12); -t193 * t192 * MDP(13) + t176 * MDP(14) + (t159 * t234 + t181 * t295) * MDP(15) - t191 * t193 * MDP(16) + (t169 * t193 - t175 * t181 + (-t235 * t217 + t157) * t234) * MDP(17) + (-t154 * t238 * t234 + (-t191 * t289 + (qJD(6) * t234 - t193) * t233) * t167) * MDP(18) + (-t191 * t234 * t177 + t158 * t239) * MDP(19); t165 * MDP(16) + (-t174 * t180 - g(3) * (t225 * t234 - t226 * t288) - t254 + (-t191 + qJD(5)) * t253) * MDP(17) + (t177 * t238 * t167 + t154 * t233) * MDP(18) + (-MDP(15) * t180 + MDP(19) * t177) * t181; -t167 * t252 * MDP(18) + t158 * MDP(19);];
tau = t1;
