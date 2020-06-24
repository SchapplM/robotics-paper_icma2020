% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = CloosQRC350OL_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:22
% EndTime: 2020-06-23 22:03:33
% DurationCPUTime: 8.70s
% Computational Cost: add. (28529->290), mult. (59877->466), div. (0->0), fcn. (44325->10), ass. (0->206)
t165 = qJDD(2) + qJDD(3);
t173 = sin(qJ(2));
t177 = cos(qJ(3));
t172 = sin(qJ(3));
t178 = cos(qJ(2));
t213 = t178 * t172;
t151 = (-t173 * t177 - t213) * qJD(1);
t212 = t178 * t177;
t234 = -t173 * t172 + t212;
t153 = t234 * qJD(1);
t221 = t153 * t151;
t130 = t165 + t221;
t237 = t130 * t177;
t162 = t178 * qJDD(1);
t208 = qJD(1) * qJD(2);
t197 = t173 * t208;
t159 = t162 - t197;
t196 = t178 * t208;
t186 = qJDD(1) * t173 + t196;
t126 = t151 * qJD(3) + t177 * t159 - t172 * t186;
t166 = qJD(2) + qJD(3);
t220 = t166 * t151;
t236 = t126 + t220;
t235 = t186 + t196;
t171 = sin(qJ(4));
t176 = cos(qJ(4));
t138 = t153 * t176 - t166 * t171;
t148 = qJD(4) + t151;
t170 = sin(qJ(5));
t175 = cos(qJ(5));
t122 = -t138 * t170 + t175 * t148;
t123 = t138 * t175 + t148 * t170;
t102 = t122 * t123;
t188 = -t126 * t171 - t165 * t176;
t107 = -qJD(4) * t138 + t188;
t106 = qJDD(5) - t107;
t87 = t102 - t106;
t149 = t151 ^ 2;
t150 = t153 ^ 2;
t164 = t166 ^ 2;
t180 = qJD(1) ^ 2;
t233 = pkin(2) * t180;
t163 = t173 * g(3);
t211 = t178 * t180;
t160 = pkin(2) * t211 - t163;
t201 = t173 * t211;
t194 = qJDD(2) - t201;
t141 = t194 * pkin(3) - t160;
t179 = qJD(2) ^ 2;
t146 = -t173 * t233 - g(3) * t178 + (-t173 ^ 2 * t180 - t179) * pkin(3);
t128 = t172 * t141 + t177 * t146;
t132 = -pkin(4) * t151 + pkin(5) * t153;
t101 = -pkin(4) * t164 - pkin(5) * t165 + t132 * t151 + t128;
t195 = t159 * t172 + t177 * t186;
t125 = -qJD(3) * t153 - t195;
t140 = qJDD(1) * pkin(2) + t235 * pkin(3);
t89 = t236 * pkin(5) + (t153 * t166 - t125) * pkin(4) + t140;
t73 = t101 * t171 + t176 * t89;
t232 = t171 * t73;
t231 = t176 * t73;
t169 = sin(qJ(6));
t174 = cos(qJ(6));
t137 = -t153 * t171 - t166 * t176;
t189 = -t126 * t176 + t165 * t171;
t108 = qJD(4) * t137 - t189;
t124 = qJDD(4) + t125;
t190 = -t108 * t175 - t124 * t170;
t86 = qJD(5) * t122 - t190;
t230 = t174 * t106 + t169 * t86;
t136 = qJD(5) - t137;
t109 = t123 * t169 + t136 * t174;
t110 = -t123 * t174 + t136 * t169;
t229 = t109 * t110;
t120 = qJD(6) + t122;
t228 = t120 * t169;
t227 = t120 * t174;
t226 = t136 * t170;
t225 = t136 * t175;
t224 = t137 * t138;
t223 = t148 * t171;
t222 = t148 * t176;
t219 = t166 * t172;
t218 = t166 * t177;
t217 = t171 * t175;
t131 = -t221 + t165;
t216 = t172 * t131;
t214 = t175 * t176;
t111 = -t123 ^ 2 - t136 ^ 2;
t210 = -qJD(4) + t148;
t209 = -qJD(6) + t120;
t207 = t170 * t229;
t206 = t175 * t229;
t205 = t171 * t102;
t204 = t176 * t102;
t203 = t172 * t224;
t202 = t177 * t224;
t200 = -pkin(4) * t177 - pkin(3);
t127 = -t177 * t141 + t172 * t146;
t100 = pkin(4) * t165 - pkin(5) * t164 - t132 * t153 - t127;
t74 = t101 * t176 - t171 * t89;
t58 = t175 * t100 - t170 * t74;
t59 = t100 * t170 + t175 * t74;
t36 = -t170 * t58 + t175 * t59;
t29 = t176 * t36 + t232;
t35 = t170 * t59 + t175 * t58;
t199 = pkin(4) * t35 - pkin(5) * t29;
t45 = t176 * t74 + t232;
t198 = pkin(4) * t100 - pkin(5) * t45;
t43 = pkin(6) * t87 + t59;
t48 = (t122 * t136 + t86) * pkin(6) + t73;
t26 = t169 * t43 + t174 * t48;
t27 = t169 * t48 - t174 * t43;
t14 = t169 * t27 + t174 * t26;
t15 = t169 * t26 - t174 * t27;
t53 = t209 * t110 + t230;
t191 = -t106 * t169 + t174 * t86;
t54 = t209 * t109 + t191;
t193 = t169 * t54 - t174 * t53;
t192 = -t171 * t74 + t231;
t187 = t127 * t177 - t128 * t172;
t85 = -qJD(5) * t123 - t108 * t170 + t175 * t124;
t113 = (-qJD(3) + t166) * t153 - t195;
t118 = -t137 ^ 2 - t138 ^ 2;
t96 = t210 * t138 + t188;
t97 = t210 * t137 + t189;
t80 = -t171 * t97 + t176 * t96;
t185 = pkin(4) * t118 - pkin(5) * t80 - t45;
t88 = -t109 ^ 2 - t110 ^ 2;
t25 = -t170 * t88 + t175 * t193;
t31 = t169 * t53 + t174 * t54;
t17 = t171 * t31 + t176 * t25;
t24 = t170 * t193 + t175 * t88;
t8 = -pkin(6) * t193 - t15;
t9 = pkin(6) * t31 + t14;
t184 = pkin(4) * t24 - pkin(5) * t17 + t176 * t8 - t9 * t217;
t77 = -t111 * t170 + t175 * t87;
t78 = (-qJD(5) - t136) * t122 + t190;
t47 = -t171 * t78 + t176 * t77;
t76 = t111 * t175 + t170 * t87;
t183 = pkin(4) * t76 - pkin(5) * t47 - t176 * t59 - t73 * t217;
t57 = t111 * pkin(6) + t58;
t11 = t15 * t170 + t175 * t57;
t12 = t15 * t175 - t170 * t57;
t6 = t12 * t176 + t14 * t171;
t182 = -pkin(5) * t6 + pkin(4) * t11 + (-t14 * t217 - t176 * t15) * pkin(6);
t181 = t178 * (pkin(4) * t172 + pkin(5) * t177) - t173 * (pkin(5) * t172 + t200) + pkin(2);
t168 = t178 ^ 2;
t143 = -t150 + t164;
t142 = t149 - t164;
t133 = t150 - t149;
t129 = -t149 - t150;
t116 = -t220 + t126;
t112 = (qJD(3) + t166) * t153 + t195;
t104 = (t137 * t176 + t138 * t171) * t148;
t103 = (-t137 * t171 + t138 * t176) * t148;
t95 = t108 * t176 - t138 * t223;
t94 = -t108 * t171 - t138 * t222;
t93 = -t107 * t171 - t137 * t222;
t92 = -t107 * t176 + t137 * t223;
t91 = (t122 * t175 + t123 * t170) * t136;
t90 = (t122 * t170 - t123 * t175) * t136;
t84 = qJDD(6) + t85;
t82 = t106 * t171 + t176 * t91;
t81 = t106 * t176 - t171 * t91;
t79 = -t171 * t96 - t176 * t97;
t72 = (-t109 * t174 - t110 * t169) * t120;
t71 = (t109 * t169 - t110 * t174) * t120;
t70 = -t123 * t226 + t175 * t86;
t69 = t123 * t225 + t170 * t86;
t68 = -t122 * t225 - t170 * t85;
t67 = -t122 * t226 + t175 * t85;
t66 = t170 * t73;
t65 = qJD(6) * t109 - t191;
t64 = -qJD(6) * t110 + t230;
t63 = t176 * t70 - t205;
t62 = t176 * t68 + t205;
t61 = -t171 * t70 - t204;
t60 = -t171 * t68 + t204;
t56 = -t170 * t84 + t175 * t72;
t55 = t170 * t72 + t175 * t84;
t52 = t110 * t228 - t174 * t65;
t51 = t110 * t227 + t169 * t65;
t50 = t109 * t227 + t169 * t64;
t49 = t109 * t228 - t174 * t64;
t46 = -t171 * t77 - t176 * t78;
t42 = t175 * t52 + t207;
t41 = t175 * t50 - t207;
t40 = t170 * t52 - t206;
t39 = t170 * t50 + t206;
t38 = t171 * t71 + t176 * t56;
t37 = -t171 * t56 + t176 * t71;
t32 = -pkin(4) * t46 + t66;
t22 = t171 * t51 + t176 * t42;
t21 = -t171 * t49 + t176 * t41;
t20 = -t171 * t42 + t176 * t51;
t19 = -t171 * t41 - t176 * t49;
t18 = pkin(5) * t46 - t171 * t59 + t73 * t214;
t16 = -t171 * t25 + t176 * t31;
t13 = t170 * pkin(6) * t14;
t7 = t170 * t9;
t5 = -t12 * t171 + t14 * t176;
t4 = -pkin(4) * t16 + t7;
t3 = -pkin(4) * t5 + t13;
t2 = pkin(5) * t16 + t171 * t8 + t9 * t214;
t1 = pkin(5) * t5 + (t14 * t214 - t15 * t171) * pkin(6);
t10 = [0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, (t159 - t197) * t178, 0, t178 * t194 - t173 * (-t168 * t180 + t179), t235 * t173, 0, 0, 0.2e1 * t186 * pkin(2), 0, (t160 + t163) * t178 - t168 * t233, pkin(2) ^ 2 * qJDD(1), t178 * (t126 * t177 - t153 * t219) - t173 * (t126 * t172 + t153 * t218), t178 * (-t112 * t177 - t172 * t236) - t173 * (-t112 * t172 + t177 * t236), t178 * (-t143 * t172 + t237) - t173 * (t130 * t172 + t143 * t177), t178 * (-t125 * t172 - t151 * t218) - t173 * (t125 * t177 - t151 * t219), t178 * (t142 * t177 - t216) - t173 * (t131 * t177 + t142 * t172), (t178 * (t151 * t177 + t153 * t172) - t173 * (t151 * t172 - t153 * t177)) * t166, t140 * t213 - t173 * (-pkin(3) * t112 - t140 * t177) + pkin(2) * t112, t140 * t212 - t173 * (-pkin(3) * t236 + t140 * t172) + pkin(2) * t236, t178 * t187 - t173 * (-pkin(3) * t129 + t127 * t172 + t128 * t177) + pkin(2) * t129, (t173 * pkin(3) + pkin(2)) * t140, t178 * (t177 * t95 + t203) - t173 * (t172 * t95 - t202), 0, 0, t178 * (t177 * t93 - t203) - t173 * (t172 * t93 + t202), 0, t178 * (t104 * t177 - t124 * t172) - t173 * (t104 * t172 + t124 * t177), 0, 0, t234 * (pkin(5) * t79 + t192) + (pkin(4) * t213 - t173 * t200 + pkin(2)) * t79, t181 * t192, t178 * (-t172 * t69 + t177 * t63) - t173 * (t172 * t63 + t177 * t69), 0, 0, t178 * (-t172 * t67 + t177 * t62) - t173 * (t172 * t62 + t177 * t67), 0, t178 * (-t172 * t90 + t177 * t82) - t173 * (t172 * t82 + t177 * t90), 0, t178 * (-t172 * t32 + t177 * t18) - t173 * (-pkin(3) * t46 + t172 * t18 + t177 * t32) + pkin(2) * t46, 0, t181 * (-t171 * t36 + t231), t178 * (-t172 * t40 + t177 * t22) - t173 * (t172 * t22 + t177 * t40), 0, 0, t178 * (-t172 * t39 + t177 * t21) - t173 * (t172 * t21 + t177 * t39), 0, t178 * (-t172 * t55 + t177 * t38) - t173 * (t172 * t38 + t177 * t55), 0, 0, t178 * (-t172 * t4 + t177 * t2) - t173 * (-pkin(3) * t16 + t172 * t2 + t177 * t4) + pkin(2) * t16, t178 * (t1 * t177 - t172 * t3) - t173 * (-pkin(3) * t5 + t1 * t172 + t177 * t3) + pkin(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, 0, t162, -t201, 0, qJDD(2), -t160, 0, 0, 0, -t221, t133, t116, t221, t113, t165, pkin(3) * (t172 * (-t164 - t149) + t237) - t127, pkin(3) * (-t216 + t177 * (-t150 - t164)) - t128, pkin(3) * (t113 * t172 - t177 * t116), -pkin(3) * t187, t94, 0, 0, t92, 0, t103, 0, 0, pkin(3) * (t118 * t177 + t172 * t80) + t185, pkin(3) * (t100 * t177 + t172 * t45) + t198, t61, 0, 0, t60, 0, t81, 0, pkin(3) * (t172 * t47 + t177 * t76) + t183, 0, pkin(3) * (t172 * t29 + t177 * t35) + t199, t20, 0, 0, t19, 0, t37, 0, 0, pkin(3) * (t17 * t172 + t177 * t24) + t184, pkin(3) * (t11 * t177 + t172 * t6) + t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, t133, t116, t221, t113, t165, -t127, -t128, 0, 0, t94, 0, 0, t92, 0, t103, 0, 0, t185, t198, t61, 0, 0, t60, 0, t81, 0, t183, 0, t199, t20, 0, 0, t19, 0, t37, 0, 0, t184, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, 0, 0, t224, 0, t124, 0, 0, 0, 0, t69, 0, 0, t67, 0, t90, 0, t66, 0, 0, t40, 0, 0, t39, 0, t55, 0, 0, t7, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, 0, 0, t102, 0, t106, 0, -t59, 0, 0, t51, 0, 0, -t49, 0, t71, 0, 0, t8, -pkin(6) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t229, 0, 0, t229, 0, t84, 0, 0, 0, 0;];
tauJ_reg = t10;
