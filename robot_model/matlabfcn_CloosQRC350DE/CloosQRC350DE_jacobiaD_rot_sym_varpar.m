% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in CloosQRC350DE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = CloosQRC350DE_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
JaD_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:15:01
	% DurationCPUTime: 1.65s
	% Computational Cost: add. (3360->97), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->96)
	t167 = sin(qJ(1));
	t163 = t167 ^ 2;
	t165 = qJ(2) + qJ(3);
	t160 = sin(t165);
	t156 = 0.1e1 / t160 ^ 2;
	t161 = cos(t165);
	t159 = t161 ^ 2;
	t214 = t156 * t159;
	t152 = t163 * t214 + 0.1e1;
	t149 = 0.1e1 / t152;
	t155 = 0.1e1 / t160;
	t169 = cos(qJ(1));
	t201 = qJD(1) * t169;
	t191 = t161 * t201;
	t162 = qJD(2) + qJD(3);
	t209 = t162 * t167;
	t194 = t156 * t209;
	t123 = (-(-t160 * t209 + t191) * t155 + t159 * t194) * t149;
	t228 = t123 - t209;
	t206 = t167 * t161;
	t148 = atan2(t206, -t160);
	t147 = cos(t148);
	t146 = sin(t148);
	t195 = t146 * t206;
	t133 = -t147 * t160 + t195;
	t130 = 0.1e1 / t133;
	t168 = cos(qJ(4));
	t203 = t168 * t169;
	t193 = t160 * t203;
	t166 = sin(qJ(4));
	t205 = t167 * t166;
	t145 = t193 - t205;
	t139 = 0.1e1 / t145;
	t131 = 0.1e1 / t133 ^ 2;
	t140 = 0.1e1 / t145 ^ 2;
	t227 = t149 - 0.1e1;
	t217 = t147 * t161;
	t118 = (t123 * t167 - t162) * t217 + (t228 * t160 + t191) * t146;
	t226 = t118 * t130 * t131;
	t184 = qJD(1) * t160 + qJD(4);
	t208 = t162 * t169;
	t192 = t161 * t208;
	t127 = -qJD(4) * t193 - t166 * t192 - t168 * t201 + t184 * t205;
	t204 = t167 * t168;
	t207 = t166 * t169;
	t144 = t160 * t207 + t204;
	t138 = t144 ^ 2;
	t137 = t138 * t140 + 0.1e1;
	t220 = t140 * t144;
	t185 = qJD(4) * t160 + qJD(1);
	t128 = -t185 * t207 + (-t184 * t167 + t192) * t168;
	t224 = t128 * t139 * t140;
	t225 = (-t127 * t220 - t138 * t224) / t137 ^ 2;
	t223 = t131 * t161;
	t222 = t131 * t169;
	t221 = t139 * t166;
	t219 = t144 * t168;
	t218 = t146 * t167;
	t216 = t155 * t159;
	t215 = t155 * t161;
	t164 = t169 ^ 2;
	t213 = t159 * t164;
	t212 = t159 * t167;
	t211 = t160 * t162;
	t210 = t161 * t162;
	t202 = qJD(1) * t167;
	t126 = t131 * t213 + 0.1e1;
	t182 = t201 * t212;
	t200 = 0.2e1 * (-t213 * t226 + (-t160 * t164 * t210 - t182) * t131) / t126 ^ 2;
	t199 = 0.2e1 * t226;
	t198 = 0.2e1 * t225;
	t158 = t161 * t159;
	t178 = t162 * (-t155 * t156 * t158 - t215);
	t197 = 0.2e1 * (t156 * t182 + t163 * t178) / t152 ^ 2;
	t196 = t161 * t222;
	t190 = 0.1e1 + t214;
	t189 = t161 * t200;
	t188 = 0.2e1 * t144 * t224;
	t187 = t161 * t197;
	t186 = t167 * t197;
	t183 = t147 * t149 * t216;
	t181 = t190 * t169;
	t180 = t184 * t169;
	t179 = t140 * t219 - t221;
	t177 = t179 * t169;
	t143 = -t160 * t204 - t207;
	t142 = -t160 * t205 + t203;
	t135 = 0.1e1 / t137;
	t134 = t190 * t167 * t149;
	t124 = 0.1e1 / t126;
	t122 = (-t227 * t161 * t146 - t167 * t183) * t169;
	t120 = -t160 * t218 - t217 + (t146 * t160 + t147 * t206) * t134;
	t119 = -t190 * t186 + (qJD(1) * t181 + 0.2e1 * t167 * t178) * t149;
	t116 = t161 * t177 * t198 + (t177 * t211 + (t179 * t202 + ((qJD(4) * t139 + t188) * t168 + (t127 * t168 + (qJD(4) * t144 - t128) * t166) * t140) * t169) * t161) * t135;
	t115 = (t120 * t223 + t130 * t160) * t169 * t200 + ((t130 * t202 + (t120 * t162 + t118) * t222) * t160 + (-t130 * t208 - (t119 * t147 * t167 + t228 * t146 + (-t123 * t218 + t146 * t162 + t147 * t201) * t134) * t196 + (t131 * t202 + t169 * t199) * t120 - ((t119 - t201) * t146 + ((-t134 * t167 + 0.1e1) * t162 + (t134 - t167) * t123) * t147) * t160 * t222) * t161) * t124;
	t1 = [t155 * t169 * t187 + (t162 * t181 + t202 * t215) * t149, t119, t119, 0, 0, 0; (t130 * t189 + (t130 * t211 + (qJD(1) * t122 + t118) * t223) * t124) * t167 + (t131 * t189 * t122 + (-((t187 - t211 + (t123 * t155 * t212 + t211) * t149) * t146 + (t186 * t216 + t123 * t161 + (t158 * t194 + (-t123 + 0.2e1 * t209) * t161) * t149) * t147) * t196 + (t131 * t211 + t161 * t199) * t122 + (-t130 + ((-t163 + t164) * t183 - t227 * t195) * t131) * t161 * qJD(1)) * t124) * t169, t115, t115, 0, 0, 0; (-t139 * t142 + t143 * t220) * t198 + (t143 * t188 - t185 * t139 * t204 + (-t162 * t206 - t180) * t221 + (t143 * t127 - t142 * t128 + t180 * t219 - (t185 * t166 - t168 * t210) * t144 * t167) * t140) * t135, t116, t116, -0.2e1 * t225 + 0.2e1 * (-t127 * t135 * t140 + (-t135 * t224 - t140 * t225) * t144) * t144, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:00
	% EndTime: 2020-06-23 21:15:03
	% DurationCPUTime: 2.99s
	% Computational Cost: add. (6566->147), mult. (9976->316), div. (1522->14), fcn. (12462->11), ass. (0->138)
	t232 = qJ(2) + qJ(3);
	t227 = cos(t232);
	t238 = cos(qJ(1));
	t301 = t227 * t238;
	t228 = qJD(2) + qJD(3);
	t237 = cos(qJ(4));
	t324 = qJD(5) * t237 + t228;
	t226 = sin(t232);
	t289 = t238 * t237;
	t266 = t226 * t289;
	t234 = sin(qJ(4));
	t235 = sin(qJ(1));
	t294 = t235 * t234;
	t215 = t266 - t294;
	t233 = sin(qJ(5));
	t236 = cos(qJ(5));
	t201 = t215 * t236 + t233 * t301;
	t196 = 0.1e1 / t201 ^ 2;
	t200 = t215 * t233 - t236 * t301;
	t292 = t236 * t200;
	t195 = 0.1e1 / t201;
	t295 = t233 * t195;
	t256 = -t196 * t292 + t295;
	t268 = t226 * t294;
	t212 = t268 - t289;
	t304 = t227 * t234;
	t206 = atan2(t212, t304);
	t203 = cos(t206);
	t202 = sin(t206);
	t313 = t202 * t212;
	t186 = t203 * t304 + t313;
	t183 = 0.1e1 / t186;
	t223 = 0.1e1 / t227;
	t229 = 0.1e1 / t234;
	t184 = 0.1e1 / t186 ^ 2;
	t224 = 0.1e1 / t227 ^ 2;
	t230 = 0.1e1 / t234 ^ 2;
	t290 = t238 * t234;
	t293 = t235 * t237;
	t214 = t226 * t290 + t293;
	t211 = t214 ^ 2;
	t181 = t211 * t184 + 0.1e1;
	t261 = qJD(1) * t226 + qJD(4);
	t297 = t228 * t238;
	t269 = t227 * t297;
	t286 = qJD(1) * t238;
	t190 = -qJD(4) * t266 - t234 * t269 - t237 * t286 + t261 * t294;
	t316 = t190 * t184;
	t210 = t212 ^ 2;
	t306 = t224 * t230;
	t207 = t210 * t306 + 0.1e1;
	t204 = 0.1e1 / t207;
	t284 = qJD(4) * t237;
	t305 = t226 * t228;
	t252 = t227 * t284 - t234 * t305;
	t272 = t212 * t306;
	t299 = t228 * t235;
	t251 = t227 * t299 + t261 * t238;
	t262 = qJD(4) * t226 + qJD(1);
	t192 = t251 * t234 + t262 * t293;
	t307 = t223 * t229;
	t275 = t192 * t307;
	t174 = (-t252 * t272 + t275) * t204;
	t250 = t174 * t212 + t252;
	t169 = (-t174 * t304 + t192) * t202 + t250 * t203;
	t185 = t183 * t184;
	t322 = t169 * t185;
	t323 = (-t211 * t322 - t214 * t316) / t181 ^ 2;
	t287 = qJD(1) * t235;
	t248 = qJD(5) * t215 + t226 * t297 + t227 * t287;
	t191 = -t262 * t290 + (-t261 * t235 + t269) * t237;
	t283 = qJD(5) * t227;
	t259 = t238 * t283 + t191;
	t176 = t259 * t233 + t248 * t236;
	t194 = t200 ^ 2;
	t189 = t194 * t196 + 0.1e1;
	t315 = t196 * t200;
	t177 = -t248 * t233 + t259 * t236;
	t197 = t195 * t196;
	t318 = t177 * t197;
	t321 = (t176 * t315 - t194 * t318) / t189 ^ 2;
	t225 = t223 * t224;
	t231 = t229 * t230;
	t320 = (t192 * t272 + (-t224 * t231 * t284 + t225 * t230 * t305) * t210) / t207 ^ 2;
	t319 = t176 * t196;
	t317 = t184 * t214;
	t302 = t227 * t237;
	t254 = -t226 * t233 + t236 * t302;
	t209 = t254 * t238;
	t314 = t200 * t209;
	t312 = t202 * t214;
	t311 = t202 * t227;
	t310 = t203 * t212;
	t309 = t203 * t214;
	t308 = t203 * t226;
	t303 = t227 * t235;
	t300 = t228 * t229;
	t298 = t228 * t237;
	t296 = t230 * t237;
	t291 = t238 * t183;
	t271 = t224 * t226 * t229;
	t257 = t212 * t271 + t235;
	t182 = t257 * t204;
	t288 = t182 - t235;
	t285 = qJD(4) * t234;
	t281 = 0.2e1 * t323;
	t280 = 0.2e1 * t321;
	t279 = 0.2e1 * t185 * t214;
	t278 = t183 * t323;
	t277 = t223 * t320;
	t276 = t184 * t312;
	t274 = t200 * t318;
	t273 = t212 * t307;
	t270 = t226 * t300;
	t265 = t230 * t284;
	t264 = t184 * t281;
	t260 = -0.2e1 * t229 * t277;
	t193 = -qJD(4) * t268 - t234 * t287 + t251 * t237;
	t258 = -t235 * t283 - t193;
	t213 = t226 * t293 + t290;
	t255 = t212 * t296 - t213 * t229;
	t253 = t226 * t236 + t233 * t302;
	t249 = -qJD(5) * t213 - t226 * t299 + t227 * t286;
	t247 = t190 * t307 - (-t223 * t265 + t224 * t270) * t214;
	t208 = t253 * t238;
	t199 = -t213 * t236 - t233 * t303;
	t198 = -t213 * t233 + t236 * t303;
	t187 = 0.1e1 / t189;
	t179 = 0.1e1 / t181;
	t178 = t255 * t223 * t204;
	t173 = (t202 + (t203 * t273 - t202) * t204) * t214;
	t172 = t182 * t310 + (-t288 * t311 - t308) * t234;
	t170 = t203 * t302 + t202 * t213 - (-t202 * t304 + t310) * t178;
	t168 = -0.2e1 * t257 * t320 + (t192 * t271 + t286 + (t223 * t300 + (-t224 * t265 + 0.2e1 * t225 * t270) * t226) * t212) * t204;
	t166 = 0.2e1 * t255 * t277 + (-t255 * t224 * t305 + (-t192 * t296 + t193 * t229 + (-t213 * t296 + (0.2e1 * t231 * t237 ^ 2 + t229) * t212) * qJD(4)) * t223) * t204;
	t165 = (-t195 * t208 + t196 * t314) * t280 + (-t209 * t319 + (-t208 * t196 + 0.2e1 * t197 * t314) * t177 + (-t253 * t195 + t254 * t315) * t287 + (((-t233 * t285 + t324 * t236) * t195 - (-t324 * t233 - t236 * t285) * t315) * t227 + t256 * t226 * (-qJD(5) - t298)) * t238) * t187;
	t164 = t172 * t214 * t264 + (-(t168 * t310 + (-t174 * t313 + t192 * t203) * t182) * t317 + (t169 * t279 + t316) * t172 + (t227 * t291 - (-t182 * t311 + t202 * t303 - t308) * t317) * t284) * t179 + (-0.2e1 * t278 * t301 + ((-t228 * t291 - (t288 * t228 + t174) * t276) * t226 + (-t183 * t287 + (-t238 * t169 - (-t168 + t286) * t312 - (-t288 * t174 - t228) * t309) * t184) * t227) * t179) * t234;
	t1 = [-t247 * t204 + t214 * t260, t168, t168, t166, 0, 0; 0.2e1 * t212 * t278 + (-t192 * t183 + (t169 * t212 + t173 * t190) * t184) * t179 + (t173 * t264 + (0.2e1 * t173 * t322 + (-t190 * t204 + t190 - (-t174 * t204 * t273 + 0.2e1 * t320) * t214) * t184 * t202 + (-(t212 * t260 + t174) * t317 + (-(-t174 + t275) * t214 + t247 * t212) * t184 * t204) * t203) * t179) * t214, t164, t164, (t170 * t317 - t183 * t215) * t281 + (t170 * t316 + t191 * t183 + (t170 * t279 - t215 * t184) * t169 - (-t227 * t285 - t226 * t298 + t166 * t212 - t178 * t192 + (t178 * t304 + t213) * t174) * t184 * t309 - (t193 + (-t166 * t234 - t174 * t237) * t227 + t250 * t178) * t276) * t179, 0, 0; (-t195 * t198 + t199 * t315) * t280 + (0.2e1 * t199 * t274 + t258 * t295 + t249 * t195 * t236 + (t249 * t200 * t233 - t199 * t176 - t198 * t177 - t258 * t292) * t196) * t187, t165, t165, t256 * t214 * t280 + (t256 * t190 + ((-qJD(5) * t195 - 0.2e1 * t274) * t236 + (t176 * t236 + (-qJD(5) * t200 + t177) * t233) * t196) * t214) * t187, -0.2e1 * t321 + 0.2e1 * (t187 * t319 + (-t187 * t318 - t196 * t321) * t200) * t200, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:01
	% EndTime: 2020-06-23 21:15:07
	% DurationCPUTime: 5.90s
	% Computational Cost: add. (17929->212), mult. (25614->390), div. (1243->12), fcn. (30756->13), ass. (0->164)
	t373 = sin(qJ(5));
	t376 = cos(qJ(5));
	t372 = qJ(2) + qJ(3);
	t369 = sin(t372);
	t370 = cos(t372);
	t375 = sin(qJ(1));
	t377 = cos(qJ(4));
	t441 = t375 * t377;
	t374 = sin(qJ(4));
	t378 = cos(qJ(1));
	t443 = t374 * t378;
	t395 = t369 * t441 + t443;
	t438 = qJD(1) * t378;
	t371 = qJD(2) + qJD(3);
	t446 = t371 * t375;
	t472 = -qJD(5) * t395 - t369 * t446 + t370 * t438;
	t409 = qJD(1) * t369 + qJD(4);
	t410 = qJD(4) * t369 + qJD(1);
	t445 = t371 * t377;
	t427 = t370 * t445;
	t435 = qJD(5) * t370;
	t440 = t378 * t377;
	t473 = -t409 * t440 + (t410 * t374 - t427 - t435) * t375;
	t304 = t473 * t373 + t472 * t376;
	t448 = t370 * t376;
	t337 = -t373 * t395 + t375 * t448;
	t335 = t337 ^ 2;
	t444 = t373 * t377;
	t450 = t369 * t376;
	t396 = t370 * t444 + t450;
	t347 = 0.1e1 / t396 ^ 2;
	t329 = t335 * t347 + 0.1e1;
	t456 = t337 * t347;
	t346 = 0.1e1 / t396;
	t412 = qJD(5) + t445;
	t403 = t412 * t373;
	t411 = qJD(5) * t377 + t371;
	t437 = qJD(4) * t374;
	t388 = -t373 * t437 + t411 * t376;
	t474 = t388 * t370;
	t318 = t369 * t403 - t474;
	t476 = t318 * t347;
	t460 = t346 * t476;
	t431 = 0.2e1 * (t304 * t456 + t335 * t460) / t329 ^ 2;
	t323 = 0.1e1 / t329;
	t401 = -t304 * t346 - t318 * t456;
	t284 = t401 * t323;
	t330 = atan2(t337, -t396);
	t321 = sin(t330);
	t322 = cos(t330);
	t406 = t321 * t396 + t322 * t337;
	t278 = t406 * t284 + t304 * t321 + t318 * t322;
	t297 = t321 * t337 - t322 * t396;
	t294 = 0.1e1 / t297;
	t295 = 0.1e1 / t297 ^ 2;
	t478 = t278 * t294 * t295;
	t477 = (t411 * t373 + t376 * t437) * t370;
	t442 = t375 * t374;
	t355 = t369 * t440 - t442;
	t447 = t370 * t378;
	t398 = -t355 * t373 + t376 * t447;
	t420 = 0.2e1 * t398 * t478;
	t352 = t369 * t442 - t440;
	t449 = t370 * t374;
	t393 = t346 * t352 + t449 * t456;
	t475 = t373 * t393;
	t457 = t337 * t346;
	t283 = (t321 + (-t322 * t457 - t321) * t323) * t398;
	t340 = t355 * t376 + t373 * t447;
	t354 = t369 * t443 + t441;
	t368 = pkin(7) * qJ(5) - qJ(6);
	t363 = sin(t368);
	t364 = cos(t368);
	t405 = t340 * t364 + t354 * t363;
	t311 = 0.1e1 / t405;
	t312 = 0.1e1 / t405 ^ 2;
	t316 = t340 * t363 - t354 * t364;
	t310 = t316 ^ 2;
	t300 = t310 * t312 + 0.1e1;
	t436 = qJD(4) * t377;
	t422 = t378 * t436;
	t425 = t371 * t443;
	t325 = -t369 * t422 - t370 * t425 - t377 * t438 + t409 * t442;
	t367 = pkin(7) * qJD(5) - qJD(6);
	t414 = t340 * t367 + t325;
	t439 = qJD(1) * t375;
	t451 = t369 * t371;
	t387 = -qJD(5) * t355 - t370 * t439 - t378 * t451;
	t326 = t395 * qJD(1) + t354 * qJD(4) - t378 * t427;
	t407 = t378 * t435 - t326;
	t303 = t387 * t373 + t407 * t376;
	t416 = t354 * t367 + t303;
	t286 = t416 * t363 + t414 * t364;
	t463 = t312 * t316;
	t429 = t286 * t463;
	t287 = t414 * t363 - t416 * t364;
	t468 = t287 * t311 * t312;
	t470 = (t310 * t468 + t429) / t300 ^ 2;
	t467 = t295 * t398;
	t385 = -t407 * t373 + t387 * t376;
	t466 = t385 * t295;
	t465 = t311 * t363;
	t464 = t311 * t364;
	t462 = t316 * t363;
	t461 = t316 * t364;
	t459 = t321 * t398;
	t458 = t322 * t398;
	t455 = t354 * t373;
	t454 = t354 * t376;
	t453 = t363 * t367;
	t452 = t364 * t367;
	t434 = qJD(5) * t376;
	t336 = t398 ^ 2;
	t293 = t295 * t336 + 0.1e1;
	t433 = 0.2e1 * (-t336 * t478 + t398 * t466) / t293 ^ 2;
	t432 = 0.2e1 * t470;
	t426 = t370 * t443;
	t419 = t346 * t431;
	t418 = -0.2e1 * t316 * t468;
	t417 = -0.2e1 * t337 * t460;
	t305 = -t472 * t373 + t473 * t376;
	t415 = t352 * t367 - t305;
	t327 = t410 * t441 + (t370 * t446 + t409 * t378) * t374;
	t338 = -t375 * t370 * t373 - t376 * t395;
	t413 = t338 * t367 + t327;
	t408 = -t367 * t454 + t326;
	t404 = t369 * t412;
	t397 = -t369 * t373 + t377 * t448;
	t402 = t367 * t426 - t397 * t439 + (-t376 * t404 - t477) * t378;
	t400 = t338 * t346 - t397 * t456;
	t342 = t396 * t375;
	t349 = t369 * t444 - t448;
	t399 = -t342 * t346 + t349 * t456;
	t392 = -t370 * t436 + t374 * t451;
	t390 = qJD(1) * t396;
	t389 = qJD(5) * t455 + t325 * t376 + t355 * t367;
	t344 = t397 * t378;
	t386 = t369 * t425 + t344 * t367 + (t374 * t439 - t422) * t370;
	t343 = t396 * t378;
	t334 = -t344 * t364 - t363 * t426;
	t333 = t344 * t363 - t364 * t426;
	t332 = -t355 * t363 + t364 * t454;
	t331 = -t355 * t364 - t363 * t454;
	t320 = t388 * t369 + t370 * t403;
	t319 = t412 * t450 + t477;
	t315 = -t338 * t364 + t352 * t363;
	t314 = t338 * t363 + t352 * t364;
	t309 = t378 * t390 + (-t373 * t404 + t474) * t375;
	t307 = t316 * pkin(7) - t364 * t398;
	t306 = t405 * pkin(7) + t363 * t398;
	t298 = 0.1e1 / t300;
	t291 = 0.1e1 / t293;
	t290 = t323 * t475;
	t289 = t399 * t323;
	t288 = t400 * t323;
	t282 = (t321 * t352 + t322 * t449) * t373 - t406 * t290;
	t281 = -t406 * t289 - t321 * t342 + t322 * t349;
	t279 = -t406 * t288 + t321 * t338 - t322 * t397;
	t276 = t399 * t431 + (t349 * t417 + t309 * t346 + (-t304 * t349 + t318 * t342 - t320 * t337) * t347) * t323;
	t275 = t400 * t431 + (-t397 * t417 - t305 * t346 + (t304 * t397 - t318 * t338 - t319 * t337) * t347) * t323;
	t274 = t431 * t475 + (-t393 * t434 + (t417 * t449 - t327 * t346 + (-t304 * t449 - t318 * t352 + t392 * t337) * t347) * t373) * t323;
	t273 = (t311 * t333 + t334 * t463) * t432 + (t334 * t418 - t402 * t465 - t386 * t464 + (-t334 * t286 - t333 * t287 - t386 * t462 + t402 * t461) * t312) * t298;
	t272 = (t281 * t467 + t294 * t343) * t433 + (t281 * t420 + (t343 * t278 - t281 * t385 - (t276 * t337 - t289 * t304 + t320 + (-t289 * t396 - t342) * t284) * t458 - (t276 * t396 + t289 * t318 - t309 + (t289 * t337 - t349) * t284) * t459) * t295 + (t318 * t378 + t375 * t390) * t294) * t291;
	t1 = [t398 * t419 + (-t346 * t385 - t398 * t476) * t323, t276, t276, t274, t275, 0; (t283 * t467 + t294 * t337) * t433 + (t283 * t420 - t304 * t294 + (t337 * t278 - 0.2e1 * t283 * t385) * t295 - ((t284 * t323 * t457 + t431) * t321 + (t337 * t419 + t284 + (-t284 + t401) * t323) * t322) * t398 * t467) * t291, t272, t272, (t282 * t467 - t294 * t455) * t433 + ((-t325 * t373 + t354 * t434) * t294 + (-t466 + t420) * t282 + (-t455 * t278 - (t352 * t434 + t274 * t396 + t290 * t318 + t327 * t373 + (t290 * t337 - t373 * t449) * t284) * t459 - (t434 * t449 + t274 * t337 - (t284 * t396 + t304) * t290 + (t284 * t352 - t392) * t373) * t458) * t295) * t291, (t279 * t467 + t294 * t340) * t433 + (t279 * t420 - t303 * t294 + (t340 * t278 - t279 * t385 - (t275 * t337 - t288 * t304 + t319 + (-t288 * t396 + t338) * t284) * t458 - (t275 * t396 + t288 * t318 + t305 + (t288 * t337 + t397) * t284) * t459) * t295) * t291, 0; (t311 * t314 + t315 * t463) * t432 + (-(-t415 * t363 + t413 * t364) * t311 + t315 * t418 + (-t314 * t287 - (t413 * t363 + t415 * t364) * t316 - t315 * t286) * t312) * t298, t273, t273, (t311 * t331 + t332 * t463) * t432 + (t332 * t418 - t408 * t464 - t389 * t465 + (-t332 * t286 - t331 * t287 + t389 * t461 - t408 * t462) * t312) * t298, (t306 * t311 + t307 * t463) * t432 + (-(t363 * t385 + t398 * t452) * t311 + t307 * t418 + (-t306 * t287 - (-t364 * t385 + t398 * t453) * t316 - t307 * t286) * t312 + (-(t303 * t364 - t325 * t363 - t340 * t453 + t354 * t452) * t311 - t429) * pkin(7)) * t298, -0.2e1 * t470 + 0.2e1 * (t286 * t298 * t312 + (t298 * t468 - t312 * t470) * t316) * t316;];
	JaD_rot = t1;
end