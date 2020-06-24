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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 21:56:19
% EndTime: 2020-06-23 21:56:44
% DurationCPUTime: 20.48s
% Computational Cost: add. (153667->243), mult. (321836->310), div. (0->0), fcn. (233776->12), ass. (0->126)
t244 = sin(qJ(3));
t245 = sin(qJ(2));
t250 = cos(qJ(3));
t251 = cos(qJ(2));
t218 = (-t244 * t245 + t250 * t251) * qJD(1);
t271 = qJD(1) * qJD(2);
t268 = t251 * t271;
t224 = -t245 * qJDD(1) - t268;
t269 = t245 * t271;
t225 = t251 * qJDD(1) - t269;
t190 = -t218 * qJD(3) + t250 * t224 - t244 * t225;
t217 = (-t244 * t251 - t245 * t250) * qJD(1);
t191 = t217 * qJD(3) + t244 * t224 + t250 * t225;
t240 = qJDD(1) * pkin(2);
t209 = t240 + (-t224 + t268) * pkin(3);
t237 = qJD(2) + qJD(3);
t166 = (t217 * t237 + t191) * pkin(5) + (t218 * t237 - t190) * pkin(4) + t209;
t253 = qJD(1) ^ 2;
t276 = t251 * t253;
t227 = -pkin(2) * t276 + t245 * g(3);
t270 = t245 * t276;
t263 = qJDD(2) - t270;
t213 = t263 * pkin(3) + t227;
t226 = -t245 * t253 * pkin(2) - t251 * g(3);
t239 = t245 ^ 2;
t277 = t239 * t253;
t215 = (-qJD(2) ^ 2 - t277) * pkin(3) + t226;
t193 = t244 * t213 + t250 * t215;
t199 = -t217 * pkin(4) + t218 * pkin(5);
t234 = t237 ^ 2;
t236 = qJDD(2) + qJDD(3);
t171 = -t234 * pkin(4) - t236 * pkin(5) + t217 * t199 + t193;
t243 = sin(qJ(4));
t249 = cos(qJ(4));
t162 = -t243 * t166 + t249 * t171;
t192 = t250 * t213 - t244 * t215;
t170 = t236 * pkin(4) - t234 * pkin(5) - t218 * t199 + t192;
t242 = sin(qJ(5));
t248 = cos(qJ(5));
t157 = t248 * t162 + t242 * t170;
t205 = t249 * t218 - t243 * t237;
t173 = -t205 * qJD(4) - t243 * t191 - t249 * t236;
t172 = qJDD(5) - t173;
t216 = qJD(4) + t217;
t183 = -t242 * t205 + t248 * t216;
t184 = t248 * t205 + t242 * t216;
t282 = t183 * t184;
t265 = -t172 + t282;
t155 = t265 * pkin(6) + t157;
t161 = t249 * t166 + t243 * t171;
t204 = -t243 * t218 - t249 * t237;
t174 = t204 * qJD(4) + t249 * t191 - t243 * t236;
t189 = qJDD(4) + t190;
t164 = t183 * qJD(5) + t248 * t174 + t242 * t189;
t201 = qJD(5) - t204;
t281 = t183 * t201;
t266 = t164 + t281;
t156 = t266 * pkin(6) + t161;
t241 = sin(qJ(6));
t247 = cos(qJ(6));
t149 = t241 * t155 + t247 * t156;
t175 = t241 * t184 + t247 * t201;
t159 = t175 * qJD(6) - t247 * t164 + t241 * t172;
t180 = qJD(6) + t183;
t284 = t175 * t180;
t147 = m(7) * t149 + (-t159 + t284) * mrSges(7,3);
t150 = -t247 * t155 + t241 * t156;
t176 = -t247 * t184 + t241 * t201;
t158 = -t176 * qJD(6) + t241 * t164 + t247 * t172;
t283 = t176 * t180;
t148 = m(7) * t150 + (t158 + t283) * mrSges(7,3);
t143 = t241 * t147 - t247 * t148;
t144 = mrSges(7,3) * t150 + Ifges(7,2) * t158 + (Ifges(7,1) - Ifges(7,3)) * t283;
t145 = -mrSges(7,3) * t149 + Ifges(7,1) * t159 + (-Ifges(7,2) + Ifges(7,3)) * t284;
t287 = -pkin(6) * t143 - mrSges(6,2) * t157 + Ifges(6,3) * t172 + t247 * t144 + t241 * t145 - (Ifges(6,1) - Ifges(6,2)) * t282;
t279 = t204 * t216;
t278 = t205 * t216;
t142 = m(6) * t157 + t265 * mrSges(6,2) + t143;
t267 = -t242 * t162 + t248 * t170;
t274 = -t184 ^ 2 - t201 ^ 2;
t152 = m(6) * t267 + m(7) * (t274 * pkin(6) + t267) + (-t175 ^ 2 - t176 ^ 2) * mrSges(7,3) + t274 * mrSges(6,2);
t136 = m(5) * t162 + t248 * t142 - t242 * t152 + (t173 + t278) * mrSges(5,3);
t262 = t247 * t147 + t241 * t148;
t140 = (-m(5) - m(6)) * t161 + (-t174 + t279) * mrSges(5,3) - t266 * mrSges(6,2) - t262;
t130 = t249 * t136 - t243 * t140;
t198 = -t217 * mrSges(4,1) + t218 * mrSges(4,2);
t211 = t237 * mrSges(4,1) - t218 * mrSges(4,3);
t128 = m(4) * t193 - t236 * mrSges(4,2) + t190 * mrSges(4,3) + t217 * t198 - t237 * t211 + t130;
t210 = -t237 * mrSges(4,2) + t217 * mrSges(4,3);
t260 = t242 * t142 + t248 * t152 + m(5) * t170 + (-t204 ^ 2 - t205 ^ 2) * mrSges(5,3);
t138 = m(4) * t192 + t236 * mrSges(4,1) - t191 * mrSges(4,3) - t218 * t198 + t237 * t210 + t260;
t275 = t244 * t128 + t250 * t138;
t273 = qJD(1) * t245;
t272 = qJD(1) * t251;
t264 = -t242 * t174 + t248 * t189;
t129 = -t243 * t136 - t249 * t140;
t261 = Ifges(7,3) * (-t184 * qJD(5) + qJDD(6) + t264) + (-Ifges(7,1) + Ifges(7,2)) * t175 * t176;
t134 = Ifges(6,1) * t164 + mrSges(6,2) * t161 - t247 * t145 + t241 * t144 + pkin(6) * t262 + (-Ifges(6,2) + Ifges(6,3)) * t281;
t154 = Ifges(6,2) * t264 + (-Ifges(6,2) * qJD(5) + (Ifges(6,1) - Ifges(6,3)) * t201) * t184 + t261;
t259 = Ifges(5,3) * t189 + t242 * t134 + t248 * t154 + (-Ifges(5,1) + Ifges(5,2)) * t204 * t205;
t257 = m(4) * t209 - t190 * mrSges(4,1) + t191 * mrSges(4,2) - t217 * t210 + t218 * t211 + t129;
t228 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t272;
t256 = m(3) * t240 - t224 * mrSges(3,1) + t228 * t272 + t257;
t131 = mrSges(5,3) * t161 + Ifges(5,1) * t174 + t248 * t134 - t242 * t154 + (-Ifges(5,2) + Ifges(5,3)) * t279;
t132 = mrSges(5,3) * t162 + Ifges(5,2) * t173 + (Ifges(5,1) - Ifges(5,3)) * t278 - t287;
t195 = Ifges(4,4) * t218 + Ifges(4,2) * t217 + Ifges(4,6) * t237;
t196 = Ifges(4,1) * t218 + Ifges(4,4) * t217 + Ifges(4,5) * t237;
t255 = pkin(4) * t260 - pkin(5) * t130 + mrSges(4,1) * t192 - mrSges(4,2) * t193 + Ifges(4,5) * t191 + Ifges(4,6) * t190 + Ifges(4,3) * t236 - t243 * t131 - t249 * t132 + t218 * t195 - t217 * t196;
t230 = Ifges(3,1) * t272 + Ifges(3,5) * qJD(2);
t254 = pkin(3) * t275 + mrSges(3,1) * t227 + Ifges(3,5) * t225 + Ifges(3,3) * qJDD(2) + t230 * t273 + t255;
t252 = cos(qJ(1));
t246 = sin(qJ(1));
t229 = Ifges(3,5) * t272 + Ifges(3,3) * qJD(2);
t194 = Ifges(4,5) * t218 + Ifges(4,6) * t217 + Ifges(4,3) * t237;
t126 = qJDD(1) * mrSges(2,1) + (-mrSges(3,3) * t239 - mrSges(2,2)) * t253 + t256;
t124 = m(3) * t226 - mrSges(3,1) * t277 + t224 * mrSges(3,3) - qJD(2) * t228 + t250 * t128 - t244 * t138;
t123 = m(3) * t227 + (-t225 - t269) * mrSges(3,3) + t263 * mrSges(3,1) + t275;
t122 = -pkin(4) * t129 - mrSges(4,1) * t209 + mrSges(4,3) * t193 + Ifges(4,4) * t191 + Ifges(4,2) * t190 + Ifges(4,6) * t236 - t218 * t194 + t237 * t196 + t259;
t121 = pkin(5) * t129 + mrSges(4,2) * t209 - mrSges(4,3) * t192 + Ifges(4,1) * t191 + Ifges(4,4) * t190 + Ifges(4,5) * t236 + t249 * t131 - t243 * t132 + t217 * t194 - t237 * t195;
t120 = -t253 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t251 * t123 + t245 * t124;
t119 = -mrSges(3,3) * t227 + Ifges(3,1) * t225 + Ifges(3,5) * qJDD(2) + t250 * t121 - t244 * t122 + (Ifges(3,2) * qJD(2) - t229) * t273;
t118 = -pkin(3) * t257 - mrSges(3,1) * t240 + mrSges(3,3) * t226 + Ifges(3,2) * t224 + qJD(2) * t230 + t244 * t121 + t250 * t122 - t229 * t272;
t117 = Ifges(2,6) * qJDD(1) + (-t251 * Ifges(3,2) * t245 + Ifges(2,5)) * t253 + mrSges(2,1) * g(3) + t254 - pkin(2) * (-t245 * t123 + t251 * t124);
t116 = -mrSges(2,2) * g(3) + Ifges(2,5) * qJDD(1) - t253 * Ifges(2,6) + t251 * t118 + t245 * t119;
t115 = Ifges(2,3) * qJDD(1) + t251 * t119 - t245 * t118 + pkin(2) * (-mrSges(3,3) * t277 + t256);
t1 = [t252 * t116 - t246 * t117 - pkin(1) * (t246 * t120 + t252 * t126), t116, t119, t121, t131, t134, t145; t246 * t116 + t252 * t117 + pkin(1) * (t252 * t120 - t246 * t126), t117, t118, t122, t132, t154, t144; t115, t115, -Ifges(3,2) * t270 + t254, t255, t259, t287, t261;];
m_new = t1;
