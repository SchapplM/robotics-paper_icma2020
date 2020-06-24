% Calculate vector of inverse dynamics joint torques for
% CloosQRC350DE
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = CloosQRC350DE_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:35
% EndTime: 2020-06-23 21:07:42
% DurationCPUTime: 4.07s
% Computational Cost: add. (2295->240), mult. (4916->366), div. (0->0), fcn. (3892->10), ass. (0->110)
t237 = cos(qJ(3));
t233 = sin(qJ(3));
t234 = sin(qJ(2));
t283 = t233 * t234;
t258 = qJD(1) * t283;
t238 = cos(qJ(2));
t278 = qJD(1) * t238;
t199 = t237 * t278 - t258;
t232 = sin(qJ(4));
t298 = -qJD(5) * t232 + t199;
t282 = t234 * t237;
t204 = t233 * t238 + t282;
t200 = qJD(1) * t204;
t231 = sin(qJ(5));
t235 = cos(qJ(5));
t236 = cos(qJ(4));
t281 = t235 * t236;
t250 = t200 * t281 + t231 * t298;
t269 = qJD(1) * qJD(2);
t256 = t238 * t269;
t297 = qJDD(1) * t234 + t256;
t268 = qJD(1) * qJD(3);
t296 = t238 * t268 + t297;
t230 = qJD(2) + qJD(3);
t192 = t199 * t232 - t230 * t236;
t191 = qJD(5) - t192;
t295 = t191 * t231;
t193 = -t199 * t236 - t230 * t232;
t198 = qJD(4) + t200;
t179 = -t193 * t231 + t235 * t198;
t227 = pkin(7) * qJD(5) - qJD(6);
t175 = t179 - t227;
t294 = pkin(3) * t233;
t293 = pkin(3) * t237;
t292 = pkin(3) * qJD(2);
t257 = t234 * t269;
t266 = qJDD(1) * t238;
t177 = t268 * t282 + (t257 - t266) * t237 + t296 * t233;
t178 = -t230 * t258 + t233 * t266 + t237 * t296;
t223 = -pkin(3) * t234 - pkin(2);
t162 = -pkin(3) * t256 - pkin(4) * t178 + pkin(5) * t177 + qJDD(1) * t223;
t279 = qJD(1) * t223;
t182 = -pkin(4) * t200 - pkin(5) * t199 + t279;
t264 = t233 * t292;
t205 = -pkin(5) * t230 + t264;
t174 = -t182 * t232 + t205 * t236;
t229 = qJDD(2) + qJDD(3);
t265 = qJDD(2) * t233;
t276 = qJD(3) * t237;
t196 = -pkin(5) * t229 + (qJD(2) * t276 + t265) * pkin(3);
t156 = qJD(4) * t174 + t162 * t236 + t196 * t232;
t291 = t156 * t232;
t180 = t193 * t235 + t198 * t231;
t228 = pkin(7) * qJ(5) - qJ(6);
t219 = sin(t228);
t290 = t180 * t219;
t289 = t191 * t235;
t203 = -t237 * t238 + t283;
t287 = t203 * t232;
t286 = t203 * t236;
t285 = t231 * t236;
t284 = t232 * t235;
t280 = t238 * qJD(1) ^ 2;
t277 = qJD(2) * t238;
t275 = qJD(4) * t232;
t274 = qJD(4) * t235;
t273 = qJD(4) * t236;
t271 = qJD(5) * t236;
t173 = t182 * t236 + t205 * t232;
t270 = t173 * qJD(4);
t263 = t237 * t292;
t262 = pkin(3) * t277;
t261 = qJD(3) * t294;
t260 = pkin(3) * t276;
t187 = -pkin(4) * t199 + pkin(5) * t200;
t224 = -pkin(5) + t294;
t254 = pkin(3) * t278 + qJD(4) * t224 - t187;
t253 = qJD(2) * (-qJD(3) + t230);
t252 = qJD(3) * (-qJD(2) - t230);
t251 = -t203 * g(3) - t200 * t279;
t190 = t230 * t283 - t237 * t277 - t238 * t276;
t249 = -t203 * t271 - t190;
t226 = qJDD(2) * t293;
t248 = t235 * (-t162 * t232 + t196 * t236 - t270) + t231 * (pkin(4) * t229 - qJD(2) * t261 + t226);
t206 = pkin(4) * t230 + t263;
t247 = t174 * t231 - t206 * t235;
t246 = t204 * g(3) + t199 * t279 + t226;
t189 = t230 * t204;
t245 = t189 * t232 + t203 * t273;
t244 = t189 * t236 - t203 * t275;
t165 = qJD(4) * t193 + t177 * t232 + t229 * t236 + qJDD(5);
t243 = -qJD(4) * t180 - qJD(5) * t295 + t165 * t235;
t242 = -qJD(5) * t204 - t244;
t154 = -qJD(5) * t247 + t248;
t168 = t174 * t235 + t206 * t231;
t241 = -t154 * t236 + (-t203 * t235 - t204 * t285) * g(3) - t250 * t173 + (t200 * t232 + t275) * t168;
t166 = qJD(4) * t192 + t177 * t236 - t229 * t232;
t176 = qJDD(4) + t178;
t158 = qJD(5) * t179 + t166 * t235 + t176 * t231;
t220 = cos(t228);
t153 = (-t191 * t227 - t158) * t220 + (t180 * t227 - t165) * t219;
t157 = -pkin(7) * qJDD(5) - qJD(5) * t180 - t166 * t231 + t176 * t235 + qJDD(6);
t164 = -t180 * t220 - t191 * t219;
t240 = (t153 * (-t219 * t236 + t220 * t284) + ((-t227 * t235 + t198) * t232 * t219 + ((-t227 + t274) * t236 + t250) * t220) * t164) * MDP(18) + (t157 * t231 * t232 + (t198 * t285 - t235 * t298) * t175) * MDP(19) + (-t158 * t284 + (-t235 * t273 - t250) * t180) * MDP(15) + (-t191 * t198 * t232 + t236 * t165) * MDP(16) + (-t193 * t198 * t236 - t166 * t232) * MDP(13) + (-t200 * t230 + t177) * MDP(8) + (-t199 * t230 + t178) * MDP(9) + (t199 ^ 2 - t200 ^ 2) * MDP(7) + t229 * MDP(10) + (-MDP(14) * t198 + MDP(6) * t200) * t199;
t225 = pkin(4) + t293;
t188 = -pkin(4) * t204 + pkin(5) * t203 + t223;
t186 = t203 * t281 + t204 * t231;
t171 = pkin(4) * t190 + pkin(5) * t189 - t262;
t161 = t231 * t249 - t235 * t242;
t1 = [qJDD(1) * MDP(1) + (-0.2e1 * t257 + t266) * t238 * MDP(2) + (qJD(2) ^ 2 * t234 - qJDD(2) * t238) * MDP(3) + 0.2e1 * t297 * MDP(5) * pkin(2) + (t177 * t203 - t189 * t199) * MDP(6) + (t177 * t204 + t178 * t203 + t189 * t200 + t190 * t199) * MDP(7) + (t189 * t230 + t203 * t229) * MDP(8) + (-t190 * t230 + t204 * t229) * MDP(9) + (0.2e1 * t200 * t262 + (qJD(1) * t190 - qJDD(1) * t204 - t178) * t223) * MDP(11) + ((-qJD(1) * t203 + t199) * t262 + (qJD(1) * t189 + qJDD(1) * t203 + t177) * t223) * MDP(12) + (t166 * t286 + t193 * t244) * MDP(13) + (t176 * t204 - t190 * t198) * MDP(14) + (t158 * t186 + t161 * t180) * MDP(15) + (t165 * t287 + t191 * t245) * MDP(16) + (t156 * t186 + t161 * t173 + (t158 * t188 + t171 * t180 + (-t168 * t203 + t188 * t289) * qJD(4)) * t236 + (-t154 * t203 - t168 * t189 + t171 * t289 + t188 * t243) * t232) * MDP(17) + ((-t153 * t186 + t164 * (-t227 * t287 - t161)) * t220 + (-t153 * t287 + t164 * (t186 * t227 - t245)) * t219) * MDP(18) + ((t157 * t204 + t175 * t249) * t235 + (-t157 * t286 + t175 * t242) * t231) * MDP(19); t240 + (t232 * t224 * t158 + (-(-t224 * t271 - t261) * t191 - t225 * t165) * t231 + (t232 * t260 + t236 * t254) * t180 + (-t291 + (-t165 * t224 - t270) * t236 + (-qJD(5) * t225 + t232 * t254 - t236 * t260) * t191) * t235 + t241) * MDP(17) - MDP(3) * t266 + ((-t200 * t278 + t229 * t237 + t233 * t252) * pkin(3) + t246) * MDP(11) + t234 * MDP(2) * t280 + (-pkin(2) * t280 + g(3) * t234) * MDP(5) + ((-t199 * t278 + (-qJDD(2) - t229) * t233 + t237 * t252) * pkin(3) + t251) * MDP(12) + qJDD(2) * MDP(4); (-(t187 * t236 + t232 * t263) * t180 + (-pkin(4) * t165 - t191 * t264) * t231 + (-t236 * t270 - t291 + (-pkin(4) * qJD(5) - t187 * t232 + t236 * t263) * t191) * t235 + ((-t191 * t274 - t158) * t232 + t243 * t236) * pkin(5) + t241) * MDP(17) + t240 + (t253 * t294 + t246) * MDP(11) + ((t237 * t253 - t265) * pkin(3) + t251) * MDP(12); -t193 * t192 * MDP(13) + t176 * MDP(14) + (t158 * t231 + t180 * t289) * MDP(15) - t191 * t193 * MDP(16) + (t168 * t193 - t174 * t180 + (g(3) * t287 + t156) * t231) * MDP(17) + (-t153 * t220 * t231 + (-t220 * t289 + (t227 * t231 + t193) * t219) * t164) * MDP(18) + (t157 * t235 - t175 * t295) * MDP(19); -t180 * t179 * MDP(15) + t165 * MDP(16) + (-t186 * g(3) - t173 * t179 - t248 + (qJD(5) - t191) * t247) * MDP(17) + (-t153 * t219 + (-pkin(7) * t290 + (pkin(7) * t191 + t175) * t220) * t164) * MDP(18) + (-pkin(7) * t157 + t175 * t180) * MDP(19); -t164 * (t191 * t220 - t290) * MDP(18) + t157 * MDP(19);];
tau = t1;
