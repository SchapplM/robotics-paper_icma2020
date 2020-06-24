% Calculate time derivative of joint inertia matrix for
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
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350OL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:56:05
% EndTime: 2020-06-23 21:56:21
% DurationCPUTime: 9.68s
% Computational Cost: add. (5046->398), mult. (12153->627), div. (0->0), fcn. (11412->10), ass. (0->214)
t141 = sin(qJ(5));
t142 = sin(qJ(4));
t146 = cos(qJ(5));
t239 = qJD(5) * t146;
t147 = cos(qJ(4));
t242 = qJD(4) * t147;
t315 = t141 * t242 + t142 * t239;
t279 = Ifges(5,2) + Ifges(6,3);
t278 = Ifges(7,3) + Ifges(6,2);
t238 = qJD(5) * t147;
t243 = qJD(4) * t146;
t312 = t141 * t238 + t142 * t243;
t140 = sin(qJ(6));
t145 = cos(qJ(6));
t235 = qJD(6) * t145;
t193 = qJD(6) + t243;
t194 = -qJD(6) * t146 - qJD(4);
t241 = qJD(5) * t141;
t252 = t145 * t147;
t53 = t193 * t252 + (t194 * t140 - t145 * t241) * t142;
t253 = t142 * t146;
t98 = t140 * t147 + t145 * t253;
t310 = t140 * t53 + t98 * t235;
t240 = qJD(5) * t142;
t208 = t141 * t240;
t254 = t142 * t145;
t54 = t194 * t254 + (-t193 * t147 + t208) * t140;
t97 = -t140 * t253 + t252;
t309 = t141 * t97 * t235 + (t141 * t54 + t239 * t97) * t140;
t256 = t141 * t142;
t149 = cos(qJ(2));
t246 = qJD(2) * t149;
t143 = sin(qJ(3));
t144 = sin(qJ(2));
t148 = cos(qJ(3));
t101 = -t143 * t149 - t144 * t148;
t303 = qJD(2) + qJD(3);
t74 = t303 * t101;
t166 = t143 * t144 - t148 * t149;
t75 = t303 * t166;
t30 = pkin(3) * t246 - pkin(4) * t75 + pkin(5) * t74;
t129 = t144 * pkin(3) + pkin(2);
t71 = -pkin(4) * t101 - pkin(5) * t166 + t129;
t308 = t30 * t256 + t315 * t71;
t211 = t142 * t242;
t199 = 0.2e1 * t211;
t280 = Ifges(7,1) - Ifges(7,2);
t136 = t142 ^ 2;
t139 = t147 ^ 2;
t247 = t136 + t139;
t203 = t247 * mrSges(5,3);
t306 = -mrSges(4,2) - t203;
t305 = t140 * t280;
t304 = -t146 * t30 + t71 * t241;
t302 = 2 * m(5);
t301 = 2 * m(6);
t300 = 2 * m(7);
t299 = -2 * mrSges(4,1);
t298 = 2 * Ifges(4,4);
t244 = qJD(4) * t142;
t258 = t147 * t74;
t157 = -qJD(5) * t101 - t166 * t244 - t258;
t180 = t166 * t238 + t75;
t17 = t180 * t141 - t157 * t146;
t297 = 0.2e1 * t17;
t296 = 0.2e1 * t129;
t295 = -0.2e1 * t147;
t294 = pkin(6) * mrSges(7,3);
t293 = m(6) * t71;
t292 = pkin(6) * t146;
t251 = t146 * t147;
t217 = t166 * t251;
t66 = t101 * t141 - t217;
t291 = Ifges(6,1) * t66;
t257 = t166 * t142;
t39 = -t140 * t257 - t145 * t66;
t290 = Ifges(7,1) * t39;
t289 = Ifges(7,1) * t98;
t38 = t140 * t66 - t166 * t254;
t288 = Ifges(7,2) * t38;
t287 = Ifges(7,2) * t97;
t162 = t142 * t74 - t166 * t242;
t158 = qJD(6) * t66 + t162;
t179 = -qJD(6) * t257 - t17;
t6 = t158 * t140 + t179 * t145;
t286 = t140 * t6;
t7 = -t179 * t140 + t158 * t145;
t285 = t140 * t7;
t284 = t145 * t6;
t283 = t145 * t7;
t282 = t30 * t71;
t127 = pkin(3) * t143 - pkin(5);
t281 = pkin(5) - t127;
t135 = t141 ^ 2;
t185 = t135 * t211;
t232 = t136 * t282;
t70 = t71 ^ 2;
t277 = t135 * t232 + t70 * t185;
t255 = t141 * t147;
t103 = pkin(4) * t146 + pkin(5) * t255;
t128 = pkin(3) * t148 + pkin(4);
t275 = pkin(3) * qJD(3);
t226 = t148 * t275;
t196 = t147 * t226;
t227 = t143 * t275;
t42 = (-t127 * t238 - t227) * t146 + (-qJD(5) * t128 + t127 * t244 - t196) * t141;
t85 = -pkin(4) * t241 + (-t141 * t244 + t146 * t238) * pkin(5);
t88 = -t127 * t255 + t128 * t146;
t276 = t103 * t42 + t85 * t88;
t274 = mrSges(6,2) * t142;
t273 = mrSges(6,2) * t147;
t272 = Ifges(6,1) * t146;
t138 = t146 ^ 2;
t271 = t138 * Ifges(6,1);
t270 = t139 * Ifges(5,1);
t269 = t140 * t38;
t267 = t145 * t39;
t266 = t145 * t53;
t265 = t145 * t54;
t264 = t145 * t98;
t262 = t146 * t71;
t261 = t147 * t30;
t260 = t147 * t66;
t259 = t147 * t71;
t112 = t127 * t251;
t89 = t141 * t128 + t112;
t248 = t136 - t139;
t245 = qJD(4) * t166;
t237 = qJD(6) * t140;
t236 = qJD(6) * t141;
t234 = 2 * mrSges(7,3);
t233 = 0.2e1 * t136;
t231 = Ifges(7,1) * t266;
t230 = Ifges(7,1) * t264;
t229 = mrSges(7,3) * t266;
t228 = mrSges(7,3) * t264;
t222 = t71 * t256;
t219 = t97 * t237;
t84 = pkin(4) * t239 + t312 * pkin(5);
t216 = -pkin(5) - t292;
t18 = t157 * t141 + t180 * t146;
t215 = t278 * t18;
t210 = t146 * t242;
t207 = t141 * t239;
t204 = t127 - t292;
t202 = t42 * t222 + t308 * t88;
t201 = t308 * t103 + t85 * t222;
t200 = -0.2e1 * t211;
t198 = t142 * t226;
t197 = t146 * t226;
t195 = mrSges(7,3) * t219;
t187 = t140 * t98 * t236;
t183 = t136 * t207;
t182 = t145 * t305;
t29 = pkin(6) * t66 + t259;
t176 = pkin(6) * t166 - t262;
t40 = t176 * t142;
t13 = t140 * t40 + t145 * t29;
t14 = t140 * t29 - t145 * t40;
t165 = -t71 * t244 + t261;
t8 = pkin(6) * t17 + t165;
t9 = t176 * t242 + (-pkin(6) * t74 + t304) * t142;
t3 = -t14 * qJD(6) + t140 * t9 + t145 * t8;
t181 = -t13 * t53 - t3 * t98;
t131 = pkin(6) * t244;
t41 = -t312 * t127 + t128 * t239 - t141 * t227 + t146 * t196;
t37 = t131 + t41;
t77 = -pkin(6) * t147 + t89;
t95 = t204 * t142;
t49 = t140 * t77 + t145 * t95;
t125 = pkin(6) * t208;
t64 = t204 * t242 + t125 + t198;
t11 = qJD(6) * t49 + t140 * t64 - t145 * t37;
t50 = t140 * t95 - t145 * t77;
t12 = -qJD(6) * t50 + t140 * t37 + t145 * t64;
t175 = t11 * t145 - t12 * t140;
t174 = t11 * t140 + t12 * t145;
t173 = t13 * t145 + t14 * t140;
t105 = t216 * t142;
t132 = t141 * pkin(4);
t96 = t132 + (-pkin(5) * t146 - pkin(6)) * t147;
t68 = t105 * t145 + t140 * t96;
t76 = t131 + t84;
t86 = t216 * t242 + t125;
t25 = qJD(6) * t68 + t140 * t86 - t145 * t76;
t69 = t105 * t140 - t145 * t96;
t26 = -qJD(6) * t69 + t140 * t76 + t145 * t86;
t172 = -t140 * t26 + t145 * t25;
t171 = t140 * t25 + t145 * t26;
t170 = -t267 + t269;
t169 = t140 * t50 + t145 * t49;
t168 = t140 * t69 + t145 * t68;
t167 = 0.2e1 * t203;
t163 = t147 * (-pkin(5) * mrSges(6,2) - t272);
t161 = -t279 * t136 - t270;
t160 = -Ifges(5,1) * t258 - t17 * t272;
t159 = (-t278 * t135 - t271) * t142;
t156 = -0.2e1 * Ifges(6,1) * t183 + 0.2e1 * t54 * t287 + 0.2e1 * t53 * t289 + t279 * t200 + (Ifges(5,1) + t271) * t199 + t278 * (0.2e1 * t183 + 0.2e1 * t185);
t2 = t13 * qJD(6) + t140 * t8 - t145 * t9;
t155 = t140 * t2 + t145 * t3 + (-t13 * t140 + t14 * t145) * qJD(6);
t154 = -t173 * qJD(6) - t140 * t3 + t145 * t2;
t153 = t135 * Ifges(6,1) * t240 + Ifges(7,1) * t187 + t309 * Ifges(7,2) + (t187 + t309) * t294 + t278 * (t138 * t240 + t141 * t210);
t152 = -Ifges(6,3) * t244 + (-t219 + t265) * Ifges(7,2) + t310 * Ifges(7,1) + (t265 + t310) * t294;
t65 = t101 * t146 + t166 * t255;
t151 = -t136 * Ifges(5,1) * t245 + Ifges(4,5) * t74 + Ifges(4,6) * t75 + t208 * t291 + t7 * t287 + t54 * t288 + t6 * t289 + t53 * t290 + (t14 * t54 + t2 * t97) * mrSges(7,3) + t279 * (-t139 * t245 + t142 * t258) + t278 * (t18 * t256 + t315 * t65);
t137 = t145 ^ 2;
t134 = t140 ^ 2;
t104 = -pkin(5) * t251 + t132;
t99 = t136 * t127 * t226;
t94 = t315 * Ifges(7,3);
t67 = t103 * t85;
t28 = t88 * t42;
t20 = t139 * t282;
t1 = [t291 * t297 + 0.2e1 * t6 * t290 + 0.2e1 * t7 * t288 + (t13 * t3 + t14 * t2 + t70 * t183 + t277) * t300 + (t138 * t232 + t20 + (t138 - 0.1e1) * t70 * t211 + t277) * t301 + (t20 + t232) * t302 + (-t13 * t6 + t14 * t7 + t2 * t38 - t3 * t39) * t234 + (mrSges(4,2) * t296 + t101 * t298 + t71 * t167) * t74 + 0.2e1 * t65 * t215 + (t129 * t299 + 0.2e1 * (Ifges(4,2) + Ifges(5,3)) * t101) * t75 + (-t30 * t167 - t298 * t75 + (-t279 * t233 - (2 * Ifges(4,1)) - 0.2e1 * t270) * t74 + (-Ifges(5,1) + t279) * t199 * t166) * t166 + (0.2e1 * pkin(2) * mrSges(3,1) + 0.2e1 * pkin(3) * (-mrSges(4,1) * t101 - mrSges(4,2) * t166) + m(4) * pkin(3) * t296 + 0.2e1 * (Ifges(3,2) - Ifges(3,1)) * t144) * t246 + ((t147 * t297 + (t146 * t74 + t166 * t241) * t233 + (-0.2e1 * t66 - 0.4e1 * t217) * t244) * t71 + 0.2e1 * (-t136 * t146 * t166 + t260) * t30) * mrSges(6,2); t151 + (m(6) * (t127 * t261 + t71 * t196 - t41 * t262 + t304 * t89) + (t127 * t17 + t166 * t41 + t66 * t226 - t89 * t74) * mrSges(6,2) + t160) * t142 + ((mrSges(6,2) * t127 - t272) * t260 + (-t248 * t127 - t89 * t251) * t293 - (-t89 * t273 + t161) * t166) * qJD(4) + m(6) * t202 + (t143 * t75 - t148 * t74 + (t101 * t148 - t143 * t166) * qJD(3)) * mrSges(4,3) * pkin(3) + m(7) * (t11 * t14 + t12 * t13 + t2 * t50 + t3 * t49 + t202) + (t11 * t38 - t12 * t39 - t49 * t6 + t50 * t7 + t181) * mrSges(7,3) - Ifges(3,5) * qJD(2) * t144; (t99 + (t127 * t139 * t148 - t128 * t143) * t275) * t302 + t156 + (t41 * t295 + (t127 * t241 - t197) * t233 + (0.2e1 * t89 - 0.4e1 * t112) * t244) * mrSges(6,2) + (t127 ^ 2 * t211 + t41 * t89 + t28 + t99) * t301 + (t11 * t50 + t12 * t49 + t28) * t300 + (t143 * t299 + 0.2e1 * t306 * t148) * t275 + (t11 * t97 - t12 * t98 - t49 * t53 + t50 * t54) * t234; t151 + (m(6) * (-pkin(5) * t261 + t304 * t104 - t84 * t262) + (-pkin(5) * t17 - t104 * t74 + t166 * t84) * mrSges(6,2) + t160) * t142 + (t66 * t163 + (t248 * pkin(5) - t104 * t251) * t293 - (-t104 * t273 + t161) * t166) * qJD(4) + m(6) * t201 + m(7) * (t13 * t26 + t14 * t25 + t2 * t69 + t3 * t68 + t201) + (t25 * t38 - t26 * t39 - t6 * t68 + t69 * t7 + t181) * mrSges(7,3); m(6) * (pkin(5) * t127 * t200 + t104 * t41 + t84 * t89 + t276) + t156 + ((-t41 - t84) * t147 + (-t281 * t241 - t197) * t136 + (0.2e1 * t251 * t281 + t104 + t89) * t244) * mrSges(6,2) + m(7) * (t11 * t69 + t12 * t68 + t25 * t50 + t26 * t49 + t276) + (-mrSges(4,1) * t143 + (-m(6) * t136 * pkin(5) + t306) * t148 + (-t247 * pkin(5) * t148 - pkin(4) * t143) * m(5)) * t275 + ((-t12 - t26) * t98 + (t11 + t25) * t97 + (t50 + t69) * t54 + (-t49 - t68) * t53) * mrSges(7,3); t156 + (t25 * t97 - t26 * t98 - t53 * t68 + t54 * t69) * t234 + (0.2e1 * t104 * t244 + t84 * t295 + (-0.2e1 * t136 * t241 + 0.4e1 * t142 * t210) * pkin(5)) * mrSges(6,2) + (t25 * t69 + t26 * t68 + t67) * t300 + (pkin(5) ^ 2 * t211 + t104 * t84 + t67) * t301; Ifges(5,3) * t75 + (t215 + (mrSges(6,2) * t259 - Ifges(7,1) * t267 + Ifges(7,2) * t269 + t291 + t173 * mrSges(7,3) + (m(7) * t173 + t170 * mrSges(7,3)) * pkin(6)) * qJD(5)) * t146 + (Ifges(6,1) * t17 - t278 * t65 * qJD(5) + (t38 * t235 + t285) * Ifges(7,2) + (t39 * t237 - t284) * Ifges(7,1) + t165 * mrSges(6,2) + t155 * mrSges(7,3) + (m(7) * t155 + (t285 - t284 + (t140 * t39 + t145 * t38) * qJD(6)) * mrSges(7,3)) * pkin(6)) * t141; t153 + (t159 + (t127 * t274 - t230 + t169 * mrSges(7,3) + (m(7) * t169 - t228) * pkin(6)) * t146) * qJD(5) + (-Ifges(6,1) * t210 - t231 + (t127 * t242 + t198) * mrSges(6,2) + ((-t140 * t49 + t145 * t50) * qJD(6) + t174) * mrSges(7,3) + (-t229 + m(7) * (t50 * t235 - t49 * t237 + t174)) * pkin(6)) * t141; t153 + (t159 + (-pkin(5) * t274 - t230 + t168 * mrSges(7,3) + (m(7) * t168 - t228) * pkin(6)) * t146) * qJD(5) + (-t231 + qJD(4) * t163 + ((-t140 * t68 + t145 * t69) * qJD(6) + t171) * mrSges(7,3) + (-t229 + m(7) * (t69 * t235 - t68 * t237 + t171)) * pkin(6)) * t141; -0.2e1 * t135 * qJD(6) * t182 + (0.2e1 * t137 * Ifges(7,1) + 0.2e1 * t134 * Ifges(7,2) + 0.2e1 * Ifges(6,1) - (2 * Ifges(6,2)) - 0.2e1 * Ifges(7,3) + (pkin(6) ^ 2 * t300 + 0.4e1 * t294) * (t134 + t137)) * t207; t162 * Ifges(6,3) + (-t38 * t237 + t283) * Ifges(7,2) + (t39 * t235 + t286) * Ifges(7,1) + (-t142 * t304 + t71 * t210) * mrSges(6,2) + t154 * mrSges(7,3) + (m(7) * t154 + (-t170 * qJD(6) + t283 + t286) * mrSges(7,3)) * pkin(6); -t41 * mrSges(6,2) + (-t169 * qJD(6) + t175) * mrSges(7,3) + (-t195 + m(7) * (-t49 * t235 - t50 * t237 + t175)) * pkin(6) + t152; -t84 * mrSges(6,2) + (-t168 * qJD(6) + t172) * mrSges(7,3) + (-t195 + m(7) * (-t68 * t235 - t69 * t237 + t172)) * pkin(6) + t152; -t182 * t239 + t280 * t236 * (t134 - t137); 0.2e1 * t235 * t305; Ifges(7,3) * t18; t94; t94; -Ifges(7,3) * t241; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
