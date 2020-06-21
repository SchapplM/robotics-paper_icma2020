% Calculate vector of cutting torques with Newton-Euler for
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = CloosQRC350OL_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:00:53
% EndTime: 2020-06-20 08:01:49
% DurationCPUTime: 33.87s
% Computational Cost: add. (256548->336), mult. (521533->431), div. (0->0), fcn. (394377->12), ass. (0->135)
t281 = qJDD(1) * pkin(2);
t292 = cos(qJ(2));
t294 = qJD(1) ^ 2;
t309 = t292 * t294;
t285 = sin(qJ(3));
t286 = sin(qJ(2));
t291 = cos(qJ(3));
t258 = (-t285 * t286 + t291 * t292) * qJD(1);
t305 = qJD(1) * qJD(2);
t304 = t292 * t305;
t267 = -t286 * qJDD(1) - t304;
t268 = t292 * qJDD(1) - t286 * t305;
t226 = -t258 * qJD(3) + t291 * t267 - t285 * t268;
t257 = (-t285 * t292 - t286 * t291) * qJD(1);
t227 = t257 * qJD(3) + t285 * t267 + t291 * t268;
t244 = t281 + (-t267 + t304) * pkin(3);
t279 = qJD(2) + qJD(3);
t187 = (t257 * t279 + t227) * pkin(5) + (t258 * t279 - t226) * pkin(4) + t244;
t270 = -pkin(2) * t309 + t286 * g(3);
t248 = (-t286 * t309 + qJDD(2)) * pkin(3) + t270;
t269 = -t286 * t294 * pkin(2) - t292 * g(3);
t252 = (-t286 ^ 2 * t294 - qJD(2) ^ 2) * pkin(3) + t269;
t229 = t285 * t248 + t291 * t252;
t237 = -t257 * pkin(4) + t258 * pkin(5);
t276 = t279 ^ 2;
t278 = qJDD(2) + qJDD(3);
t197 = -t276 * pkin(4) - t278 * pkin(5) + t257 * t237 + t229;
t284 = sin(qJ(4));
t290 = cos(qJ(4));
t178 = -t284 * t187 + t290 * t197;
t228 = t291 * t248 - t285 * t252;
t196 = t278 * pkin(4) - t276 * pkin(5) - t258 * t237 + t228;
t283 = sin(qJ(5));
t289 = cos(qJ(5));
t169 = t289 * t178 + t283 * t196;
t240 = t290 * t258 - t284 * t279;
t204 = -t240 * qJD(4) - t284 * t227 - t290 * t278;
t203 = qJDD(5) - t204;
t253 = qJD(4) + t257;
t219 = -t283 * t240 + t289 * t253;
t220 = t289 * t240 + t283 * t253;
t165 = (t219 * t220 - t203) * pkin(6) + t169;
t177 = t290 * t187 + t284 * t197;
t239 = -t284 * t258 - t290 * t279;
t205 = t239 * qJD(4) + t290 * t227 - t284 * t278;
t225 = qJDD(4) + t226;
t185 = t219 * qJD(5) + t289 * t205 + t283 * t225;
t238 = qJD(5) - t239;
t166 = (t219 * t238 + t185) * pkin(6) + t177;
t282 = sin(qJ(6));
t288 = cos(qJ(6));
t163 = t282 * t165 + t288 * t166;
t206 = t282 * t220 + t288 * t238;
t173 = t206 * qJD(6) - t288 * t185 + t282 * t203;
t184 = -t220 * qJD(5) - t283 * t205 + t289 * t225;
t183 = qJDD(6) + t184;
t207 = -t288 * t220 + t282 * t238;
t186 = -t206 * mrSges(7,1) + t207 * mrSges(7,2);
t217 = qJD(6) + t219;
t188 = -t217 * mrSges(7,2) + t206 * mrSges(7,3);
t160 = m(7) * t163 + t183 * mrSges(7,1) - t173 * mrSges(7,3) - t207 * t186 + t217 * t188;
t164 = -t288 * t165 + t282 * t166;
t172 = -t207 * qJD(6) + t282 * t185 + t288 * t203;
t189 = t217 * mrSges(7,1) - t207 * mrSges(7,3);
t161 = m(7) * t164 - t183 * mrSges(7,2) + t172 * mrSges(7,3) + t206 * t186 - t217 * t189;
t154 = t282 * t160 - t288 * t161;
t198 = -t219 * mrSges(6,1) + t220 * mrSges(6,2);
t209 = t238 * mrSges(6,1) - t220 * mrSges(6,3);
t153 = m(6) * t169 - t203 * mrSges(6,2) + t184 * mrSges(6,3) + t219 * t198 - t238 * t209 + t154;
t168 = -t283 * t178 + t289 * t196;
t167 = (-t220 ^ 2 - t238 ^ 2) * pkin(6) + t168;
t208 = -t238 * mrSges(6,2) + t219 * mrSges(6,3);
t158 = m(6) * t168 + m(7) * t167 + t203 * mrSges(6,1) - t172 * mrSges(7,1) + t173 * mrSges(7,2) - t185 * mrSges(6,3) - t206 * t188 + t207 * t189 - t220 * t198 + t238 * t208;
t214 = -t239 * mrSges(5,1) + t240 * mrSges(5,2);
t231 = t253 * mrSges(5,1) - t240 * mrSges(5,3);
t147 = m(5) * t178 - t225 * mrSges(5,2) + t204 * mrSges(5,3) + t289 * t153 - t283 * t158 + t239 * t214 - t253 * t231;
t230 = -t253 * mrSges(5,2) + t239 * mrSges(5,3);
t303 = t288 * t160 + t282 * t161;
t149 = t225 * mrSges(5,1) + t184 * mrSges(6,1) - t185 * mrSges(6,2) - t205 * mrSges(5,3) + t219 * t208 - t220 * t209 - t240 * t214 + t253 * t230 + (-m(5) - m(6)) * t177 - t303;
t140 = t290 * t147 - t284 * t149;
t236 = -t257 * mrSges(4,1) + t258 * mrSges(4,2);
t246 = t279 * mrSges(4,1) - t258 * mrSges(4,3);
t138 = m(4) * t229 - t278 * mrSges(4,2) + t226 * mrSges(4,3) + t257 * t236 - t279 * t246 + t140;
t245 = -t279 * mrSges(4,2) + t257 * mrSges(4,3);
t302 = m(5) * t196 - t204 * mrSges(5,1) + t205 * mrSges(5,2) + t283 * t153 + t289 * t158 - t239 * t230 + t240 * t231;
t145 = m(4) * t228 + t278 * mrSges(4,1) - t227 * mrSges(4,3) - t258 * t236 + t279 * t245 + t302;
t308 = t285 * t138 + t291 * t145;
t307 = qJD(1) * t286;
t306 = qJD(1) * t292;
t139 = -t284 * t147 - t290 * t149;
t180 = Ifges(7,4) * t207 + Ifges(7,2) * t206 + Ifges(7,6) * t217;
t181 = Ifges(7,1) * t207 + Ifges(7,4) * t206 + Ifges(7,5) * t217;
t301 = mrSges(7,1) * t163 - mrSges(7,2) * t164 + Ifges(7,5) * t173 + Ifges(7,6) * t172 + Ifges(7,3) * t183 + t207 * t180 - t206 * t181;
t179 = Ifges(7,5) * t207 + Ifges(7,6) * t206 + Ifges(7,3) * t217;
t155 = -mrSges(7,1) * t167 + mrSges(7,3) * t164 + Ifges(7,4) * t173 + Ifges(7,2) * t172 + Ifges(7,6) * t183 - t207 * t179 + t217 * t181;
t156 = mrSges(7,2) * t167 - mrSges(7,3) * t163 + Ifges(7,1) * t173 + Ifges(7,4) * t172 + Ifges(7,5) * t183 + t206 * t179 - t217 * t180;
t190 = Ifges(6,5) * t220 + Ifges(6,6) * t219 + Ifges(6,3) * t238;
t191 = Ifges(6,4) * t220 + Ifges(6,2) * t219 + Ifges(6,6) * t238;
t143 = pkin(6) * t303 + mrSges(6,2) * t177 - mrSges(6,3) * t168 + Ifges(6,1) * t185 + Ifges(6,4) * t184 + Ifges(6,5) * t203 + t282 * t155 - t288 * t156 + t219 * t190 - t238 * t191;
t192 = Ifges(6,1) * t220 + Ifges(6,4) * t219 + Ifges(6,5) * t238;
t151 = -mrSges(6,1) * t177 + mrSges(6,3) * t169 + Ifges(6,4) * t185 + Ifges(6,2) * t184 + Ifges(6,6) * t203 - t220 * t190 + t238 * t192 + t301;
t211 = Ifges(5,4) * t240 + Ifges(5,2) * t239 + Ifges(5,6) * t253;
t212 = Ifges(5,1) * t240 + Ifges(5,4) * t239 + Ifges(5,5) * t253;
t300 = -mrSges(5,1) * t177 - mrSges(5,2) * t178 + Ifges(5,5) * t205 + Ifges(5,6) * t204 + Ifges(5,3) * t225 + t283 * t143 + t289 * t151 + t240 * t211 - t239 * t212;
t299 = m(4) * t244 - t226 * mrSges(4,1) + t227 * mrSges(4,2) - t257 * t245 + t258 * t246 + t139;
t210 = Ifges(5,5) * t240 + Ifges(5,6) * t239 + Ifges(5,3) * t253;
t135 = mrSges(5,2) * t196 + mrSges(5,3) * t177 + Ifges(5,1) * t205 + Ifges(5,4) * t204 + Ifges(5,5) * t225 + t289 * t143 - t283 * t151 + t239 * t210 - t253 * t211;
t296 = pkin(6) * t154 - mrSges(6,1) * t168 + mrSges(6,2) * t169 - Ifges(6,5) * t185 - Ifges(6,6) * t184 - Ifges(6,3) * t203 - t288 * t155 - t282 * t156 - t220 * t191 + t219 * t192;
t141 = -mrSges(5,1) * t196 + mrSges(5,3) * t178 + Ifges(5,4) * t205 + Ifges(5,2) * t204 + Ifges(5,6) * t225 - t240 * t210 + t253 * t212 + t296;
t233 = Ifges(4,4) * t258 + Ifges(4,2) * t257 + Ifges(4,6) * t279;
t234 = Ifges(4,1) * t258 + Ifges(4,4) * t257 + Ifges(4,5) * t279;
t298 = pkin(4) * t302 - pkin(5) * t140 + mrSges(4,1) * t228 - mrSges(4,2) * t229 + Ifges(4,5) * t227 + Ifges(4,6) * t226 + Ifges(4,3) * t278 - t284 * t135 - t290 * t141 + t258 * t233 - t257 * t234;
t271 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t307;
t272 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t306;
t297 = m(3) * t281 - t267 * mrSges(3,1) + t268 * mrSges(3,2) + t271 * t307 + t272 * t306 + t299;
t255 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t292 - Ifges(3,2) * t286) * qJD(1);
t256 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t292 - Ifges(3,4) * t286) * qJD(1);
t295 = pkin(3) * t308 + mrSges(3,1) * t270 - mrSges(3,2) * t269 + Ifges(3,5) * t268 + Ifges(3,6) * t267 + Ifges(3,3) * qJDD(2) + t255 * t306 + t256 * t307 + t298;
t293 = cos(qJ(1));
t287 = sin(qJ(1));
t265 = (mrSges(3,1) * t286 + mrSges(3,2) * t292) * qJD(1);
t254 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t292 - Ifges(3,6) * t286) * qJD(1);
t232 = Ifges(4,5) * t258 + Ifges(4,6) * t257 + Ifges(4,3) * t279;
t136 = qJDD(1) * mrSges(2,1) - t294 * mrSges(2,2) + t297;
t133 = m(3) * t270 + qJDD(2) * mrSges(3,1) - t268 * mrSges(3,3) + qJD(2) * t271 - t265 * t306 + t308;
t132 = m(3) * t269 - qJDD(2) * mrSges(3,2) + t267 * mrSges(3,3) - qJD(2) * t272 + t291 * t138 - t285 * t145 - t265 * t307;
t131 = -pkin(4) * t139 - mrSges(4,1) * t244 + mrSges(4,3) * t229 + Ifges(4,4) * t227 + Ifges(4,2) * t226 + Ifges(4,6) * t278 - t258 * t232 + t279 * t234 + t300;
t130 = pkin(5) * t139 + mrSges(4,2) * t244 - mrSges(4,3) * t228 + Ifges(4,1) * t227 + Ifges(4,4) * t226 + Ifges(4,5) * t278 + t290 * t135 - t284 * t141 + t257 * t232 - t279 * t233;
t129 = -t294 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t286 * t132 + t292 * t133;
t128 = mrSges(3,2) * t281 - mrSges(3,3) * t270 + Ifges(3,1) * t268 + Ifges(3,4) * t267 + Ifges(3,5) * qJDD(2) - qJD(2) * t255 + t291 * t130 - t285 * t131 - t254 * t307;
t127 = -pkin(3) * t299 - mrSges(3,1) * t281 + mrSges(3,3) * t269 + Ifges(3,4) * t268 + Ifges(3,2) * t267 + Ifges(3,6) * qJDD(2) + qJD(2) * t256 + t285 * t130 + t291 * t131 - t254 * t306;
t126 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + t295 - pkin(2) * (t292 * t132 - t286 * t133) + t294 * Ifges(2,5);
t125 = -mrSges(2,2) * g(3) + Ifges(2,5) * qJDD(1) - t294 * Ifges(2,6) + t292 * t127 + t286 * t128;
t124 = pkin(2) * t297 + Ifges(2,3) * qJDD(1) - t286 * t127 + t292 * t128;
t1 = [-mrSges(1,2) * g(3) + t293 * t125 - t287 * t126 - pkin(1) * (t287 * t129 + t293 * t136), t125, t128, t130, t135, t143, t156; mrSges(1,1) * g(3) + t287 * t125 + t293 * t126 + pkin(1) * (t293 * t129 - t287 * t136), t126, t127, t131, t141, t151, t155; t124, t124, t295, t298, t300, -t296, t301;];
m_new = t1;
