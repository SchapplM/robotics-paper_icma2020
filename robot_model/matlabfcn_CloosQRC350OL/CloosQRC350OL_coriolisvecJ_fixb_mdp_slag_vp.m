% Calculate minimal parameter regressor of Coriolis joint torque vector for
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = CloosQRC350OL_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:22
% EndTime: 2020-06-23 22:03:26
% DurationCPUTime: 2.84s
% Computational Cost: add. (1786->199), mult. (4290->327), div. (0->0), fcn. (3386->10), ass. (0->104)
t200 = sin(qJ(3));
t261 = sin(qJ(2));
t227 = qJD(1) * t261;
t221 = t200 * t227;
t204 = cos(qJ(3));
t205 = cos(qJ(2));
t244 = qJD(1) * t205;
t228 = t204 * t244;
t177 = t221 - t228;
t214 = -t200 * t205 - t204 * t261;
t178 = t214 * qJD(1);
t198 = sin(qJ(5));
t202 = cos(qJ(5));
t203 = cos(qJ(4));
t248 = t202 * t203;
t162 = t177 * t198 + t178 * t248;
t199 = sin(qJ(4));
t237 = qJD(5) * t199;
t263 = t198 * t237 - t162;
t262 = pkin(2) * MDP(5);
t195 = qJD(2) + qJD(3);
t173 = -t177 * t203 - t195 * t199;
t176 = qJD(4) + t178;
t158 = -t173 * t198 + t202 * t176;
t192 = t261 * pkin(3) + pkin(2);
t260 = pkin(3) * qJD(2);
t259 = pkin(3) * qJD(3);
t169 = t195 * t214;
t164 = t169 * qJD(1);
t247 = t204 * t205;
t217 = t195 * t247;
t226 = qJD(2) * t261;
t220 = qJD(1) * t226;
t245 = -qJD(3) * t221 - t200 * t220;
t165 = qJD(1) * t217 + t245;
t243 = qJD(2) * t205;
t225 = qJD(1) * t243;
t147 = pkin(3) * t225 + pkin(4) * t165 + pkin(5) * t164;
t184 = qJD(1) * pkin(2) + pkin(3) * t227;
t161 = -pkin(4) * t178 - pkin(5) * t177 + t184;
t235 = t200 * t260;
t182 = -pkin(5) * t195 + t235;
t154 = -t161 * t199 + t182 * t203;
t231 = qJD(2) * t259;
t222 = t204 * t231;
t140 = t154 * qJD(4) + t147 * t203 + t199 * t222;
t258 = t140 * t199;
t172 = t177 * t199 - t195 * t203;
t171 = qJD(5) - t172;
t257 = t171 * t202;
t255 = t176 * t198;
t254 = t177 * t184;
t253 = t178 * t184;
t230 = t200 * t261;
t181 = -t230 + t247;
t252 = t181 * t199;
t251 = t181 * t203;
t151 = t172 * qJD(4) + t164 * t203;
t137 = t158 * qJD(5) + t151 * t202 - t165 * t198;
t250 = t199 * t137;
t201 = cos(qJ(6));
t249 = t201 * t202;
t153 = t161 * t203 + t182 * t199;
t242 = qJD(4) * t153;
t241 = qJD(4) * t199;
t240 = qJD(4) * t202;
t239 = qJD(4) * t203;
t238 = qJD(5) * t198;
t236 = qJD(5) * t203;
t157 = qJD(6) + t158;
t234 = t204 * t260;
t233 = pkin(3) * t243;
t232 = t204 * t259;
t167 = -pkin(4) * t177 + t178 * pkin(5);
t190 = pkin(3) * t200 - pkin(5);
t224 = -pkin(3) * t244 + qJD(4) * t190 - t167;
t223 = qJD(3) * (-qJD(2) - t195);
t219 = (-qJD(3) + t195) * t260;
t170 = -qJD(3) * t230 - t200 * t226 + t217;
t218 = -t181 * t236 - t170;
t183 = pkin(4) * t195 + t234;
t216 = t154 * t198 - t183 * t202;
t159 = t173 * t202 + t255;
t197 = sin(qJ(6));
t215 = t159 * t197 + t171 * t201;
t213 = t169 * t199 + t181 * t239;
t212 = t169 * t203 - t181 * t241;
t211 = -t198 * t200 * t231 + t202 * (-t147 * t199 + t203 * t222 - t242);
t136 = -t216 * qJD(5) + t211;
t148 = t154 * t202 + t183 * t198;
t210 = -t136 * t203 + t263 * t153 + (t178 * t199 + t241) * t148;
t152 = t173 * qJD(4) + t164 * t199;
t209 = -qJD(4) * t159 + t152 * t202 - t171 * t238;
t208 = -qJD(5) * t214 - t212;
t135 = t215 * qJD(6) - t137 * t201 + t152 * t197;
t138 = -t159 * qJD(5) - t151 * t198 - t165 * t202;
t145 = -t159 * t201 + t171 * t197;
t207 = (t135 * (t197 * t203 + t199 * t249) + ((t162 + (qJD(6) + t240) * t203) * t201 + (-t201 * t238 + (-qJD(6) * t202 - t176) * t197) * t199) * t145) * MDP(18) + (t138 * t198 * t199 + ((-t177 + t237) * t202 + t203 * t255) * t157) * MDP(19) + (-t202 * t250 + (-t202 * t239 + t263) * t159) * MDP(15) + (-t176 * t199 * t171 + t152 * t203) * MDP(16) + (-t176 * t203 * t173 - t151 * t199) * MDP(13) + (-t178 * t195 + t164) * MDP(8) + (-t245 + (-t177 - t228) * t195) * MDP(9) + (t177 ^ 2 - t178 ^ 2) * MDP(7) + (-t176 * MDP(14) + MDP(6) * t178) * t177;
t191 = pkin(3) * t204 + pkin(4);
t168 = -pkin(4) * t214 + pkin(5) * t181 + t192;
t166 = t181 * t248 + t198 * t214;
t149 = pkin(4) * t170 + pkin(5) * t169 + t233;
t143 = t218 * t198 - t208 * t202;
t1 = [-0.2e1 * t205 * MDP(2) * t220 - qJD(2) ^ 2 * t261 * MDP(3) + 0.2e1 * t225 * t262 + (t164 * t181 - t169 * t177) * MDP(6) + (t164 * t214 - t165 * t181 + t169 * t178 + t170 * t177) * MDP(7) + (t165 * t192 + t170 * t184 - 0.2e1 * t178 * t233) * MDP(11) + (t164 * t192 + t169 * t184 + (qJD(1) * t181 - t177) * t233) * MDP(12) + (t151 * t251 + t212 * t173) * MDP(13) + (-t165 * t214 - t170 * t176) * MDP(14) + (t137 * t166 + t143 * t159) * MDP(15) + (t152 * t252 + t213 * t171) * MDP(16) + (t140 * t166 + t143 * t153 + (t137 * t168 + t149 * t159 + (-t148 * t181 + t168 * t257) * qJD(4)) * t203 + (-t136 * t181 - t148 * t169 + t149 * t257 + t209 * t168) * t199) * MDP(17) + ((-t135 * t166 + t145 * (qJD(6) * t252 - t143)) * t201 + (t135 * t252 + t145 * (qJD(6) * t166 + t213)) * t197) * MDP(18) + ((t138 * t214 + t157 * t218) * t202 + (-t138 * t251 + t157 * t208) * t198) * MDP(19) + (t169 * MDP(8) - t170 * MDP(9)) * t195; (t190 * t250 + (-(-t190 * t236 - t200 * t259) * t171 - t191 * t152) * t198 + (t199 * t232 + t224 * t203) * t159 + (-t258 + (-t152 * t190 - t242) * t203 + (-qJD(5) * t191 + t224 * t199 - t203 * t232) * t171) * t202 + t210) * MDP(17) + (t254 + (t178 * t244 + t200 * t223) * pkin(3)) * MDP(11) + (-t253 + (t177 * t244 + t204 * t223) * pkin(3)) * MDP(12) + t207 + (t261 * MDP(2) - t262) * t205 * qJD(1) ^ 2; (t200 * t219 + t254) * MDP(11) + (t204 * t219 - t253) * MDP(12) + (-(t167 * t203 + t199 * t234) * t159 + (-pkin(4) * t152 - t171 * t235) * t198 + (-t153 * t239 - t258 + (-pkin(4) * qJD(5) - t167 * t199 + t203 * t234) * t171) * t202 + ((-t171 * t240 - t137) * t199 + t209 * t203) * pkin(5) + t210) * MDP(17) + t207; -t173 * t172 * MDP(13) - t165 * MDP(14) + (t137 * t198 + t159 * t257) * MDP(15) - t171 * t173 * MDP(16) + (t140 * t198 + t148 * t173 - t154 * t159) * MDP(17) + (-t135 * t201 * t198 + (-t171 * t249 + (qJD(6) * t198 - t173) * t197) * t145) * MDP(18) + (-t171 * t198 * t157 + t138 * t202) * MDP(19); t152 * MDP(16) + (-t153 * t158 - t211 + (qJD(5) - t171) * t216) * MDP(17) + (t157 * t201 * t145 + t135 * t197) * MDP(18) + (-MDP(15) * t158 + MDP(19) * t157) * t159; -t145 * t215 * MDP(18) + t138 * MDP(19);];
tauc = t1;
