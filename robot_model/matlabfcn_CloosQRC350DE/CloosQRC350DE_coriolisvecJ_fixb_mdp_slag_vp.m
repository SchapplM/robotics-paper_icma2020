% Calculate minimal parameter regressor of Coriolis joint torque vector for
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = CloosQRC350DE_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350DE_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:35
% EndTime: 2020-06-23 21:07:41
% DurationCPUTime: 3.61s
% Computational Cost: add. (1924->201), mult. (4431->328), div. (0->0), fcn. (3443->10), ass. (0->99)
t196 = sin(qJ(3));
t200 = cos(qJ(3));
t201 = cos(qJ(2));
t235 = qJD(1) * t201;
t197 = sin(qJ(2));
t236 = qJD(1) * t197;
t173 = -t196 * t236 + t200 * t235;
t195 = sin(qJ(4));
t254 = -qJD(5) * t195 + t173;
t176 = t196 * t201 + t197 * t200;
t174 = qJD(1) * t176;
t194 = sin(qJ(5));
t198 = cos(qJ(5));
t199 = cos(qJ(4));
t239 = t198 * t199;
t212 = t174 * t239 + t254 * t194;
t193 = qJD(2) + qJD(3);
t253 = pkin(2) * MDP(5);
t168 = t173 * t195 - t193 * t199;
t167 = qJD(5) - t168;
t252 = t167 * t194;
t169 = -t173 * t199 - t193 * t195;
t172 = qJD(4) + t174;
t154 = -t169 * t194 + t198 * t172;
t191 = pkin(7) * qJD(5) - qJD(6);
t153 = t154 - t191;
t240 = t196 * t197;
t175 = -t200 * t201 + t240;
t251 = qJD(1) * t175;
t250 = pkin(3) * qJD(2);
t249 = pkin(3) * qJD(3);
t234 = qJD(2) * t201;
t217 = qJD(1) * t234;
t233 = qJD(3) * t201;
t160 = t193 * t200 * t236 + (qJD(1) * t233 + t217) * t196;
t161 = t193 * t251;
t143 = -pkin(3) * t217 + pkin(4) * t161 + pkin(5) * t160;
t188 = -pkin(3) * t197 - pkin(2);
t237 = qJD(1) * t188;
t157 = -pkin(4) * t174 - pkin(5) * t173 + t237;
t226 = t196 * t250;
t177 = -pkin(5) * t193 + t226;
t150 = -t157 * t195 + t177 * t199;
t222 = qJD(2) * t249;
t214 = t200 * t222;
t136 = t150 * qJD(4) + t143 * t199 + t195 * t214;
t248 = t136 * t195;
t244 = t172 * t194;
t155 = t169 * t198 + t244;
t192 = pkin(7) * qJ(5) - qJ(6);
t186 = sin(t192);
t247 = t155 * t186;
t246 = t167 * t198;
t243 = t175 * t195;
t242 = t175 * t199;
t241 = t195 * t198;
t149 = t157 * t199 + t177 * t195;
t232 = qJD(4) * t149;
t231 = qJD(4) * t195;
t230 = qJD(4) * t198;
t229 = qJD(4) * t199;
t227 = qJD(5) * t199;
t225 = t200 * t250;
t224 = pkin(3) * t234;
t223 = t200 * t249;
t221 = t173 * t237;
t220 = t174 * t237;
t163 = -pkin(4) * t173 + t174 * pkin(5);
t189 = pkin(3) * t196 - pkin(5);
t216 = pkin(3) * t235 + qJD(4) * t189 - t163;
t215 = qJD(3) * (-qJD(2) - t193);
t213 = (-qJD(3) + t193) * t250;
t166 = t193 * t240 + (-t233 - t234) * t200;
t211 = -t175 * t227 - t166;
t178 = pkin(4) * t193 + t225;
t210 = t150 * t194 - t178 * t198;
t165 = t193 * t176;
t209 = t165 * t195 + t175 * t229;
t208 = t165 * t199 - t175 * t231;
t207 = -t194 * t196 * t222 + t198 * (-t143 * t195 + t199 * t214 - t232);
t132 = -t210 * qJD(5) + t207;
t144 = t150 * t198 + t178 * t194;
t206 = -t132 * t199 - t212 * t149 + (t174 * t195 + t231) * t144;
t148 = t169 * qJD(4) + t160 * t195;
t205 = -qJD(4) * t155 - qJD(5) * t252 + t148 * t198;
t204 = -qJD(5) * t176 - t208;
t147 = t168 * qJD(4) + t160 * t199;
t133 = t154 * qJD(5) + t147 * t198 - t161 * t194;
t187 = cos(t192);
t131 = (-t167 * t191 - t133) * t187 + (t155 * t191 - t148) * t186;
t134 = -t155 * qJD(5) - t147 * t194 - t161 * t198;
t141 = -t155 * t187 - t167 * t186;
t203 = (t131 * (-t186 * t199 + t187 * t241) + ((-t191 * t198 + t172) * t195 * t186 + ((-t191 + t230) * t199 + t212) * t187) * t141) * MDP(18) + (t134 * t194 * t195 + (-t254 * t198 + t199 * t244) * t153) * MDP(19) + (-t133 * t241 + (-t198 * t229 - t212) * t155) * MDP(15) + (-t172 * t195 * t167 + t148 * t199) * MDP(16) + (-t172 * t199 * t169 - t147 * t195) * MDP(13) + (-t174 * t193 + t160) * MDP(8) + (-t173 * t193 - t161) * MDP(9) + (t173 ^ 2 - t174 ^ 2) * MDP(7) + (-t172 * MDP(14) + MDP(6) * t174) * t173;
t190 = pkin(3) * t200 + pkin(4);
t164 = -pkin(4) * t176 + pkin(5) * t175 + t188;
t162 = t175 * t239 + t176 * t194;
t145 = pkin(4) * t166 + pkin(5) * t165 - t224;
t139 = t211 * t194 - t204 * t198;
t1 = [0.2e1 * t217 * t253 + (t160 * t175 - t165 * t173) * MDP(6) + (t160 * t176 - t161 * t175 + t165 * t174 + t166 * t173) * MDP(7) + ((qJD(1) * t166 + t161) * t188 + 0.2e1 * t174 * t224) * MDP(11) + ((qJD(1) * t165 + t160) * t188 + (t173 - t251) * t224) * MDP(12) + (t147 * t242 + t208 * t169) * MDP(13) + (-t161 * t176 - t166 * t172) * MDP(14) + (t133 * t162 + t139 * t155) * MDP(15) + (t148 * t243 + t209 * t167) * MDP(16) + (t136 * t162 + t139 * t149 + (t133 * t164 + t145 * t155 + (-t144 * t175 + t164 * t246) * qJD(4)) * t199 + (-t132 * t175 - t144 * t165 + t145 * t246 + t205 * t164) * t195) * MDP(17) + ((-t131 * t162 + t141 * (-t191 * t243 - t139)) * t187 + (-t131 * t243 + t141 * (t162 * t191 - t209)) * t186) * MDP(18) + ((t134 * t176 + t153 * t211) * t198 + (-t134 * t242 + t153 * t204) * t194) * MDP(19) + (qJD(2) ^ 2 * MDP(3) - 0.2e1 * MDP(2) * t217) * t197 + (t165 * MDP(8) - t166 * MDP(9)) * t193; (t195 * t189 * t133 + (-(-t189 * t227 - t196 * t249) * t167 - t190 * t148) * t194 + (t195 * t223 + t216 * t199) * t155 + (-t248 + (-t148 * t189 - t232) * t199 + (-qJD(5) * t190 + t216 * t195 - t199 * t223) * t167) * t198 + t206) * MDP(17) + (-t220 + (-t173 * t235 + t200 * t215) * pkin(3)) * MDP(12) + (t221 + (-t174 * t235 + t196 * t215) * pkin(3)) * MDP(11) + t203 + (t197 * MDP(2) - t253) * t201 * qJD(1) ^ 2; (t196 * t213 + t221) * MDP(11) + (t200 * t213 - t220) * MDP(12) + (-(t163 * t199 + t195 * t225) * t155 + (-pkin(4) * t148 - t167 * t226) * t194 + (-t149 * t229 - t248 + (-pkin(4) * qJD(5) - t163 * t195 + t199 * t225) * t167) * t198 + ((-t167 * t230 - t133) * t195 + t205 * t199) * pkin(5) + t206) * MDP(17) + t203; -t169 * t168 * MDP(13) - t161 * MDP(14) + (t133 * t194 + t155 * t246) * MDP(15) - t167 * t169 * MDP(16) + (t136 * t194 + t144 * t169 - t150 * t155) * MDP(17) + (-t131 * t187 * t194 + (-t187 * t246 + (t191 * t194 + t169) * t186) * t141) * MDP(18) + (t134 * t198 - t153 * t252) * MDP(19); -t155 * t154 * MDP(15) + t148 * MDP(16) + (-t149 * t154 - t207 + (qJD(5) - t167) * t210) * MDP(17) + (-t131 * t186 + (-pkin(7) * t247 + (pkin(7) * t167 + t153) * t187) * t141) * MDP(18) + (-pkin(7) * t134 + t153 * t155) * MDP(19); -t141 * (t167 * t187 - t247) * MDP(18) + t134 * MDP(19);];
tauc = t1;
