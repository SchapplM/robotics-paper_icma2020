% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = CloosQRC350OL_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:56:05
% EndTime: 2020-06-23 21:56:32
% DurationCPUTime: 19.26s
% Computational Cost: add. (9608->485), mult. (22902->718), div. (0->0), fcn. (16865->10), ass. (0->240)
t190 = sin(qJ(3));
t194 = cos(qJ(3));
t195 = cos(qJ(2));
t330 = sin(qJ(2));
t163 = -t190 * t195 - t194 * t330;
t159 = t163 * qJD(1);
t248 = t190 * t330;
t278 = qJD(1) * t195;
t160 = -qJD(1) * t248 + t194 * t278;
t136 = pkin(4) * t160 + t159 * pkin(5);
t189 = sin(qJ(4));
t193 = cos(qJ(4));
t318 = pkin(3) * qJD(2);
t255 = t194 * t318;
t113 = -t136 * t189 + t193 * t255;
t188 = sin(qJ(5));
t192 = cos(qJ(5));
t256 = t190 * t318;
t270 = qJD(5) * t192;
t269 = qJD(5) * t193;
t273 = qJD(4) * t192;
t343 = t188 * t269 + t189 * t273;
t296 = -pkin(4) * t270 - pkin(5) * t343 + t192 * t113 - t188 * t256;
t265 = qJD(2) + qJD(3);
t213 = t189 * t160 + t193 * t265;
t350 = qJD(5) + t213;
t206 = t189 * t350;
t279 = t192 * t193;
t120 = t159 * t279 - t160 * t188;
t271 = qJD(5) * t188;
t245 = t189 * t271;
t357 = (-t120 + t245) * pkin(6);
t274 = qJD(4) * t189;
t184 = pkin(6) * t274;
t288 = t159 * t189;
t356 = pkin(6) * t288 + t184 - t296;
t112 = t136 * t193 + t189 * t255;
t328 = pkin(6) * t192;
t249 = -pkin(5) - t328;
t272 = qJD(4) * t193;
t355 = t249 * t272 - t112 + t357;
t140 = t265 * t163;
t128 = t140 * qJD(1);
t164 = t194 * t195 - t248;
t141 = t265 * t164;
t129 = t141 * qJD(1);
t211 = -pkin(5) * t265 + t256;
t277 = qJD(2) * t195;
t262 = pkin(3) * t277;
t339 = pkin(4) * t129 + pkin(5) * t128 + qJD(1) * t262 + qJD(4) * t211;
t182 = t330 * pkin(3) + pkin(2);
t172 = t182 * qJD(1);
t204 = -t159 * pkin(4) + t160 * pkin(5) + t172;
t252 = qJD(3) * t318;
t344 = -qJD(4) * t204 + t194 * t252;
t39 = t189 * t344 + t193 * t339;
t89 = t189 * t211 + t193 * t204;
t216 = t189 * t39 + t89 * t272;
t90 = -t189 * t204 + t193 * t211;
t354 = t90 * t274 - t216;
t353 = t188 * t272 + t189 * t270;
t180 = pkin(3) * t190 - pkin(5);
t181 = pkin(3) * t194 + pkin(4);
t275 = qJD(3) * t194;
t260 = pkin(3) * t275;
t236 = t193 * t260;
t276 = qJD(3) * t190;
t261 = pkin(3) * t276;
t102 = -t180 * t343 + t181 * t270 - t188 * t261 + t192 * t236;
t122 = pkin(3) * t278 + t136;
t352 = -(-pkin(6) * t159 - t192 * t122) * t189 + t184 + t102;
t235 = t189 * t260;
t241 = t180 - t328;
t351 = -t122 * t193 + t241 * t272 + t235 + t357;
t238 = t189 * t265;
t144 = t160 * t193 - t238;
t266 = -t159 - qJD(4);
t106 = -t144 * t188 - t192 * t266;
t105 = qJD(6) + t106;
t313 = Ifges(7,3) * t105;
t226 = Ifges(6,2) * t106 + t313;
t107 = t144 * t192 - t188 * t266;
t79 = -qJD(4) * t213 + t193 * t128;
t33 = -qJD(5) * t107 - t129 * t192 - t188 * t79;
t349 = (Ifges(6,2) + Ifges(7,3)) * t33;
t348 = -Ifges(6,3) * t144 - t226 * t188;
t187 = sin(qJ(6));
t191 = cos(qJ(6));
t185 = t188 * pkin(4);
t158 = t185 + (-pkin(5) * t192 - pkin(6)) * t193;
t166 = t249 * t189;
t218 = t158 * t191 - t166 * t187;
t347 = qJD(6) * t218 + t187 * t356 + t355 * t191;
t132 = t158 * t187 + t166 * t191;
t346 = -qJD(6) * t132 - t355 * t187 + t191 * t356;
t232 = t190 * t252;
t38 = -t339 * t189 + t193 * t344;
t168 = pkin(4) * t265 + t255;
t73 = t188 * t168 + t192 * t90;
t18 = -qJD(5) * t73 - t188 * t38 - t192 * t232;
t285 = t188 * t193;
t72 = t168 * t192 - t188 * t90;
t345 = t18 * (pkin(4) * t192 + pkin(5) * t285) + (-pkin(4) * t271 + (-t188 * t274 + t192 * t269) * pkin(5) + t113 * t188 + t192 * t256) * t72;
t151 = t180 * t279 + t188 * t181;
t145 = -pkin(6) * t193 + t151;
t156 = t241 * t189;
t219 = t145 * t191 - t156 * t187;
t342 = qJD(6) * t219 + t187 * t352 + t351 * t191;
t109 = t145 * t187 + t156 * t191;
t341 = -qJD(6) * t109 - t351 * t187 + t191 * t352;
t17 = t168 * t270 + t192 * t38 + (-qJD(5) * t90 - t232) * t188;
t199 = t107 * pkin(6) + t89;
t80 = -qJD(4) * t238 + t128 * t189 + t160 * t272;
t340 = -pkin(6) * t80 - qJD(6) * t199 + t17;
t289 = t144 * t193;
t292 = t213 * t189;
t220 = -t289 - t292;
t299 = t193 * t80;
t304 = t189 * t79;
t338 = qJD(4) * t220 + t299 - t304;
t337 = m(5) / 0.2e1;
t336 = -m(6) / 0.2e1;
t335 = m(6) / 0.2e1;
t334 = pkin(6) * m(7);
t333 = pkin(2) * mrSges(3,1);
t331 = t160 / 0.2e1;
t329 = pkin(6) * t106;
t68 = t191 * t107 - t187 * t350;
t327 = Ifges(7,1) * t68;
t67 = t187 * t107 + t191 * t350;
t326 = Ifges(7,2) * t67;
t325 = t73 * mrSges(6,2);
t295 = qJD(6) * t67;
t32 = qJD(5) * t106 - t129 * t188 + t192 * t79;
t8 = t187 * t80 - t191 * t32 + t295;
t324 = t8 * t191;
t323 = t89 * mrSges(6,2);
t294 = qJD(6) * t68;
t9 = t187 * t32 + t191 * t80 + t294;
t322 = t9 * t187;
t321 = t90 * mrSges(5,3);
t320 = t18 * (-t180 * t285 + t181 * t192) + t72 * ((-t180 * t269 - t261) * t192 + (-qJD(5) * t181 + t180 * t274 - t236) * t188);
t317 = Ifges(4,4) * t160;
t315 = Ifges(5,3) * t266;
t312 = t107 * Ifges(6,1);
t311 = t112 * t89;
t310 = (t163 * t192 - t164 * t285) * t33;
t308 = t159 * mrSges(4,3);
t307 = t160 * mrSges(4,3);
t306 = t187 * t67;
t303 = t189 * t89;
t302 = t191 * t68;
t301 = t192 * t89;
t300 = t193 * t38;
t298 = t193 * t89;
t297 = t89 * t120;
t293 = t140 * t189;
t291 = t213 * t192;
t287 = t187 * t189;
t286 = t188 * t189;
t283 = t189 * t191;
t282 = t189 * t192;
t281 = t191 * t192;
t280 = t191 * t193;
t268 = qJD(6) * t187;
t267 = qJD(6) * t191;
t264 = t189 * t325;
t263 = -t330 / 0.2e1;
t259 = Ifges(5,2) * t292;
t257 = t72 * t286;
t254 = Ifges(3,2) * t330;
t253 = t330 * Ifges(3,1);
t242 = t192 * t272;
t240 = m(4) * t172 - mrSges(4,1) * t159 + mrSges(4,2) * t160;
t137 = -pkin(4) * t163 + pkin(5) * t164 + t182;
t78 = pkin(4) * t141 + pkin(5) * t140 + t262;
t239 = t78 * t257 + (t18 * t286 + t353 * t72) * t137;
t234 = -qJD(6) * t192 - qJD(4);
t233 = qJD(6) + t273;
t230 = qJD(2) * t263;
t209 = -qJD(5) * t163 - t140 * t193 + t164 * t274;
t228 = -t164 * t269 - t141;
t51 = t188 * t228 - t192 * t209;
t229 = qJD(6) * t164 * t189 - t51;
t227 = -pkin(6) * t164 - t137 * t192;
t57 = -pkin(6) * t350 + t73;
t22 = t187 * t57 + t191 * t199;
t23 = -t187 * t199 + t191 * t57;
t225 = t187 * t23 - t191 * t22;
t223 = -t189 * t90 + t298;
t101 = t227 * t189;
t131 = t163 * t188 + t164 * t279;
t75 = pkin(6) * t131 + t137 * t193;
t222 = t101 * t191 - t187 * t75;
t44 = t101 * t187 + t191 * t75;
t221 = t189 * t144 - t193 * t213;
t215 = t193 * t39 - t274 * t89;
t214 = -Ifges(7,1) * t302 - Ifges(7,2) * t306;
t116 = t233 * t280 + (t187 * t234 - t191 * t271) * t189;
t162 = t187 * t193 + t189 * t281;
t11 = t32 * pkin(6) + t39;
t3 = t191 * t11 + t187 * t340 + t267 * t57;
t84 = t120 * t187 + t159 * t283;
t212 = -t116 * t22 - t162 * t3 + t23 * t84;
t208 = qJD(6) * t131 + t164 * t272 + t293;
t205 = t192 * t350;
t2 = t187 * t11 - t191 * t340 + t268 * t57;
t203 = t187 * t2 + t191 * t3 + (-t187 * t22 - t191 * t23) * qJD(6);
t202 = -t300 + t354;
t200 = t223 * t337 + (-t282 * t73 + t298) * t335;
t117 = t234 * t283 + (-t193 * t233 + t245) * t187;
t123 = Ifges(4,2) * t159 + Ifges(4,6) * t265 + t317;
t153 = Ifges(4,4) * t159;
t124 = Ifges(4,1) * t160 + Ifges(4,5) * t265 + t153;
t161 = -t187 * t282 + t280;
t85 = -t120 * t191 + t159 * t287;
t197 = (t312 + t323) * t245 + Ifges(5,2) * t299 + t255 * t308 - t160 * t315 + (t159 * t285 + t160 * t192 + t353) * t226 + (-t159 * t298 + t354) * mrSges(5,3) + t159 * t264 - t159 * t259 + (-t116 + t85) * t327 - (-Ifges(4,2) * t160 + t124 + t153) * t159 / 0.2e1 + (-t23 * t117 + t2 * t161 + t22 * t85) * mrSges(7,3) + (-t84 + t117) * t326 + t286 * t349 - t172 * (mrSges(4,1) * t160 + mrSges(4,2) * t159) + (t206 * t266 + t299) * Ifges(6,3) + (t264 - t259) * qJD(4) - t160 * (Ifges(4,1) * t159 - t317) / 0.2e1 + (t266 * t289 - t304) * Ifges(5,1) + (-t32 * t282 + (-t120 - t242) * t107) * Ifges(6,1) - t265 * (Ifges(4,5) * t159 - Ifges(4,6) * t160) / 0.2e1 + t288 * t321 + t9 * Ifges(7,2) * t161 + t8 * Ifges(7,1) * t162 + t123 * t331 + Ifges(4,5) * t128 - Ifges(4,6) * t129;
t171 = Ifges(3,1) * t278 + Ifges(3,5) * qJD(2);
t165 = -pkin(5) * t279 + t185;
t147 = mrSges(4,1) * t265 - t307;
t146 = -mrSges(4,2) * t265 + t308;
t98 = -t131 * t191 + t164 * t287;
t97 = t131 * t187 + t164 * t283;
t87 = t144 * t187 + t213 * t281;
t86 = t144 * t191 - t187 * t291;
t63 = -pkin(6) * t291 + t90;
t62 = -pkin(6) * t144 - t301;
t53 = t122 * t257;
t52 = t188 * t209 + t192 * t228;
t47 = t187 * t329 - t191 * t72;
t46 = t187 * t72 + t191 * t329;
t26 = t187 * t63 - t191 * t62;
t25 = t187 * t62 + t191 * t63;
t24 = t227 * t272 + (-pkin(6) * t140 + t137 * t271 - t192 * t78) * t189;
t21 = pkin(6) * t51 - t137 * t274 + t193 * t78;
t20 = -t187 * t229 + t191 * t208;
t19 = t187 * t208 + t191 * t229;
t5 = qJD(6) * t222 + t187 * t24 + t191 * t21;
t4 = qJD(6) * t44 + t187 * t21 - t191 * t24;
t1 = [(-t140 * t194 - t141 * t190) * mrSges(4,3) * t318 + (-t19 * t68 + t8 * t98) * Ifges(7,1) + t141 * t315 + (-t19 * t22 + t2 * t97 - t20 * t23 - t222 * t9 - t3 * t98 + t4 * t67 - t44 * t8 + t5 * t68) * mrSges(7,3) + m(7) * (-t2 * t222 + t22 * t5 - t23 * t4 + t3 * t44 + t239) + t265 * (Ifges(4,5) * t140 - Ifges(4,6) * t141) / 0.2e1 + (Ifges(4,1) * t140 - Ifges(4,4) * t141) * t331 + t159 * (Ifges(4,4) * t140 - Ifges(4,2) * t141) / 0.2e1 + t172 * (mrSges(4,1) * t141 + mrSges(4,2) * t140) + (Ifges(4,1) * t128 - Ifges(4,4) * t129 + (mrSges(5,3) * t39 + t79 * Ifges(5,1)) * t193 + (mrSges(4,2) * t278 + mrSges(4,3) * t276) * t318 + (-t17 * mrSges(6,2) - t38 * mrSges(5,3) + (Ifges(5,2) + Ifges(6,3)) * t80) * t189 + ((-mrSges(5,3) * t89 - Ifges(5,1) * t144) * t189 + (Ifges(5,2) * t213 + Ifges(6,3) * t350 - t321 - t325) * t193) * qJD(4)) * t164 + (m(5) * (-t189 * t38 - t272 * t90 + t215) + m(6) * (-t17 * t282 + t215 + (-t242 + t245) * t73) + (-qJD(4) * t221 + t189 * t80 + t193 * t79) * mrSges(5,3) + (t80 * t282 + t193 * t32 - t206 * t271 + (-t189 * t107 + t193 * t205) * qJD(4)) * mrSges(6,2)) * t137 + (t107 * t51 + t131 * t32) * Ifges(6,1) + (t223 * t140 - t220 * t78) * mrSges(5,3) + (Ifges(4,4) * t128 - (Ifges(4,2) + Ifges(5,3)) * t129 + (-mrSges(4,1) * t278 + mrSges(4,3) * t275) * t318) * t163 + t182 * (mrSges(4,1) * t129 + mrSges(4,2) * t128) + (t20 * t67 + t9 * t97) * Ifges(7,2) + 0.2e1 * t200 * t78 + (t240 * pkin(3) + (0.2e1 * t254 - 0.3e1 / 0.2e1 * t253 + 0.2e1 * t333 + m(4) * pkin(3) * t182) * qJD(1)) * t277 + (t105 * t52 + t310) * Ifges(7,3) + (t106 * t52 + t310) * Ifges(6,2) + (-t73 * t293 + t39 * t131 + t89 * t51 + (t193 * t107 + t189 * t205) * t78) * mrSges(6,2) + t140 * t259 + Ifges(6,3) * t140 * t206 + t171 * t230 + qJD(2) ^ 2 * Ifges(3,5) * t263 + Ifges(5,1) * t140 * t289 + t140 * t124 / 0.2e1 - t141 * t123 / 0.2e1 + m(6) * t239; m(6) * (t102 * t73 + t151 * t17 + t320) + 0.2e1 * (t216 * t335 - m(5) * t202 / 0.2e1) * t180 + t197 + (t220 * t122 - t180 * t338 - t300) * mrSges(5,3) - m(6) * t53 - 0.2e1 * t200 * t122 + (-t240 * t278 + (-t194 * t128 + (qJD(2) * t160 - t129) * t190) * mrSges(4,3) + ((-qJD(2) * mrSges(4,1) - t147 + m(5) * (-qJD(2) * t181 - t168)) * t190 + (-qJD(2) * mrSges(4,2) + t146 + t221 * mrSges(5,3) + m(6) * t303 + m(5) * (t193 * t90 + t303)) * t194) * qJD(3)) * pkin(3) + (-t297 - t151 * t80 + t189 * t180 * t32 - t17 * t193 - t102 * t350 + (-t122 * t206 - t216) * t192 + (t235 + (qJD(4) * t180 - t122) * t193) * t107) * mrSges(6,2) + (-t109 * t8 - t219 * t9 - t341 * t67 + t342 * t68 + t212) * mrSges(7,3) + (Ifges(3,5) * t230 + t330 * t171 / 0.2e1 + (-t254 + t253 / 0.2e1 - t333) * t278) * qJD(1) + (t109 * t3 - t2 * t219 + t22 * t342 + t341 * t23 + t320 - t53) * m(7); (-t132 * t8 - t218 * t9 - t346 * t67 + t347 * t68 + t212) * mrSges(7,3) + t197 + ((-mrSges(4,2) * qJD(3) - t146) * t194 + (-mrSges(4,1) * qJD(3) + t147 + t307) * t190) * t318 + (-t112 * t107 - t297 - t165 * t80 + t296 * qJD(5) + (-pkin(5) * t32 + t160 * t296 - t39 * t192) * t189 + (-t17 + (-pkin(5) * t107 - t301) * qJD(4) + t296 * t265) * t193) * mrSges(6,2) + 0.2e1 * (t202 * t337 + t216 * t336) * pkin(5) + (pkin(5) * t338 - t112 * t144 + t113 * t213 - t300) * mrSges(5,3) + (t132 * t3 - t218 * t2 + t22 * t347 + t346 * t23 + t345) * m(7) + (t165 * t17 - t296 * t73 - t311 + t345) * m(6) + (-pkin(4) * t232 - t113 * t90 + t168 * t256 - t311) * m(5); -t86 * t326 + t87 * t327 - t129 * Ifges(5,3) + t192 * t349 - m(7) * (t22 * t25 - t23 * t26) - m(6) * (-t192 * t73 + t90) * t89 + (-t90 * t107 + t73 * t144) * mrSges(6,2) + (t39 * mrSges(6,2) + t32 * Ifges(6,1) + (t267 * t67 + t322) * Ifges(7,2) + (-t268 * t68 - t324) * Ifges(7,1) + 0.2e1 * (t336 - m(7) / 0.2e1) * t89 * t72 + t203 * t334) * t188 + ((-t225 * t334 - t214 + t312) * t192 + t348) * qJD(5) + (t22 * t87 + t23 * t86 - t25 * t68 - t26 * t67 + ((t302 + t306) * pkin(6) - t225) * t270 + ((t322 - t324 + (-t187 * t68 + t191 * t67) * qJD(6)) * pkin(6) + t203) * t188) * mrSges(7,3) + (t192 * t312 - (-Ifges(5,1) + Ifges(5,2)) * t144 + t348) * t213; t107 * t313 + t80 * Ifges(6,3) + (t191 * t9 - t268 * t67) * Ifges(7,2) + (t187 * t8 - t267 * t68) * Ifges(7,1) + (t350 * t72 - t17) * mrSges(6,2) - m(7) * (t22 * t46 - t23 * t47 - t73 * t72) + (qJD(6) * t225 - t187 * t3 + t191 * t2) * t334 + (-t323 + (-Ifges(6,1) + Ifges(6,2)) * t107 + t214) * t106 + (-t46 * t68 - t47 * t67 + (t2 - t105 * t22 + (t9 - t294) * pkin(6)) * t191 + (-t3 + t105 * t23 + (t8 - t295) * pkin(6)) * t187) * mrSges(7,3); Ifges(7,3) * t33 + (Ifges(7,1) - Ifges(7,2)) * t68 * t67;];
tauc = t1(:);
