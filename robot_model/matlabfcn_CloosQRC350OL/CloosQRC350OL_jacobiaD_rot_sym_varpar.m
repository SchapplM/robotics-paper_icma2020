% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% CloosQRC350OL
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
%   Wie in CloosQRC350OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = CloosQRC350OL_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350OL_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
JaD_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:17
	% DurationCPUTime: 1.69s
	% Computational Cost: add. (3645->97), mult. (3810->204), div. (753->12), fcn. (4455->9), ass. (0->96)
	t162 = qJD(2) + qJD(3);
	t167 = sin(qJ(1));
	t209 = t162 * t167;
	t230 = 0.2e1 * t209;
	t163 = t167 ^ 2;
	t165 = qJ(2) + qJ(3);
	t160 = sin(t165);
	t156 = 0.1e1 / t160 ^ 2;
	t161 = cos(t165);
	t159 = t161 ^ 2;
	t214 = t156 * t159;
	t151 = t163 * t214 + 0.1e1;
	t148 = 0.1e1 / t151;
	t155 = 0.1e1 / t160;
	t169 = cos(qJ(1));
	t201 = qJD(1) * t169;
	t192 = t161 * t201;
	t195 = t156 * t209;
	t122 = (-(t160 * t209 - t192) * t155 - t159 * t195) * t148;
	t229 = t122 + t209;
	t168 = cos(qJ(4));
	t203 = t168 * t169;
	t166 = sin(qJ(4));
	t205 = t167 * t166;
	t144 = t160 * t203 + t205;
	t206 = t167 * t161;
	t147 = atan2(-t206, -t160);
	t146 = cos(t147);
	t145 = sin(t147);
	t196 = t145 * t206;
	t132 = -t146 * t160 - t196;
	t129 = 0.1e1 / t132;
	t138 = 0.1e1 / t144;
	t130 = 0.1e1 / t132 ^ 2;
	t139 = 0.1e1 / t144 ^ 2;
	t228 = t148 - 0.1e1;
	t217 = t146 * t161;
	t117 = (-t122 * t167 - t162) * t217 + (t229 * t160 - t192) * t145;
	t227 = t117 * t129 * t130;
	t178 = t160 * t205 + t203;
	t208 = t162 * t169;
	t193 = t161 * t208;
	t126 = t178 * qJD(1) - qJD(4) * t144 - t166 * t193;
	t204 = t167 * t168;
	t207 = t166 * t169;
	t143 = t160 * t207 - t204;
	t137 = t143 ^ 2;
	t136 = t137 * t139 + 0.1e1;
	t220 = t139 * t143;
	t185 = -qJD(1) * t160 + qJD(4);
	t186 = qJD(4) * t160 - qJD(1);
	t127 = -t186 * t207 + (t185 * t167 + t193) * t168;
	t225 = t127 * t138 * t139;
	t226 = (-t126 * t220 - t137 * t225) / t136 ^ 2;
	t158 = t161 * t159;
	t215 = t155 * t161;
	t181 = t155 * t156 * t158 + t215;
	t212 = t159 * t167;
	t183 = t201 * t212;
	t224 = (-t181 * t163 * t162 + t156 * t183) / t151 ^ 2;
	t223 = t130 * t161;
	t222 = t130 * t169;
	t221 = t138 * t166;
	t219 = t143 * t168;
	t218 = t145 * t167;
	t216 = t155 * t159;
	t164 = t169 ^ 2;
	t213 = t159 * t164;
	t211 = t160 * t162;
	t210 = t161 * t162;
	t202 = qJD(1) * t167;
	t125 = t130 * t213 + 0.1e1;
	t200 = 0.2e1 * (-t213 * t227 + (-t160 * t164 * t210 - t183) * t130) / t125 ^ 2;
	t199 = 0.2e1 * t227;
	t198 = 0.2e1 * t226;
	t197 = t161 * t222;
	t191 = 0.1e1 + t214;
	t190 = t161 * t200;
	t189 = -0.2e1 * t161 * t224;
	t188 = 0.2e1 * t167 * t224;
	t187 = 0.2e1 * t143 * t225;
	t184 = t146 * t148 * t216;
	t182 = t191 * t169;
	t180 = t185 * t169;
	t179 = t139 * t219 - t221;
	t177 = t179 * t169;
	t142 = -t160 * t204 + t207;
	t134 = 0.1e1 / t136;
	t133 = t191 * t167 * t148;
	t123 = 0.1e1 / t125;
	t121 = (t228 * t161 * t145 - t167 * t184) * t169;
	t119 = t160 * t218 - t217 - (t145 * t160 - t146 * t206) * t133;
	t118 = t191 * t188 + (-qJD(1) * t182 + t181 * t230) * t148;
	t115 = t161 * t177 * t198 + (t177 * t211 + (t179 * t202 + ((qJD(4) * t138 + t187) * t168 + (t126 * t168 + (qJD(4) * t143 - t127) * t166) * t139) * t169) * t161) * t134;
	t114 = (t119 * t223 + t129 * t160) * t169 * t200 + ((t129 * t202 + (t119 * t162 + t117) * t222) * t160 + (-t129 * t208 - (-t118 * t146 * t167 + t229 * t145 - (t122 * t218 + t145 * t162 - t146 * t201) * t133) * t197 + (t130 * t202 + t169 * t199) * t119 - ((t118 + t201) * t145 + ((-t133 * t167 + 0.1e1) * t162 + (-t133 + t167) * t122) * t146) * t160 * t222) * t161) * t123;
	t1 = [t155 * t169 * t189 + (-t162 * t182 - t202 * t215) * t148, t118, t118, 0, 0, 0; (t129 * t190 + (t129 * t211 + (qJD(1) * t121 + t117) * t223) * t123) * t167 + (t130 * t190 * t121 + (-((t189 + t211 + (t122 * t155 * t212 - t211) * t148) * t145 + (t188 * t216 - t122 * t161 + (t158 * t195 + (t122 + t230) * t161) * t148) * t146) * t197 + (t130 * t211 + t161 * t199) * t121 + (-t129 + ((-t163 + t164) * t184 + t228 * t196) * t130) * t161 * qJD(1)) * t123) * t169, t114, t114, 0, 0, 0; (t138 * t178 + t142 * t220) * t198 + (t142 * t187 - t186 * t138 * t204 + (-t162 * t206 + t180) * t221 + (t142 * t126 + t178 * t127 - t180 * t219 - (t186 * t166 - t168 * t210) * t143 * t167) * t139) * t134, t115, t115, -0.2e1 * t226 + 0.2e1 * (-t126 * t134 * t139 + (-t134 * t225 - t139 * t226) * t143) * t143, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:19
	% DurationCPUTime: 3.08s
	% Computational Cost: add. (6566->151), mult. (9976->322), div. (1522->14), fcn. (12462->11), ass. (0->138)
	t240 = qJ(2) + qJ(3);
	t235 = cos(t240);
	t246 = cos(qJ(1));
	t309 = t235 * t246;
	t236 = qJD(2) + qJD(3);
	t245 = cos(qJ(4));
	t334 = qJD(5) * t245 + t236;
	t234 = sin(t240);
	t299 = t245 * t246;
	t242 = sin(qJ(4));
	t243 = sin(qJ(1));
	t302 = t243 * t242;
	t220 = t234 * t299 + t302;
	t241 = sin(qJ(5));
	t244 = cos(qJ(5));
	t205 = t220 * t244 + t241 * t309;
	t200 = 0.1e1 / t205 ^ 2;
	t204 = t220 * t241 - t244 * t309;
	t300 = t244 * t204;
	t199 = 0.1e1 / t205;
	t303 = t241 * t199;
	t263 = -t200 * t300 + t303;
	t216 = t234 * t302 + t299;
	t291 = qJD(4) * t246;
	t270 = t245 * t291;
	t293 = qJD(4) * t242;
	t271 = t243 * t293;
	t305 = t236 * t246;
	t274 = t235 * t305;
	t194 = t216 * qJD(1) - t234 * t270 - t242 * t274 - t271;
	t297 = t246 * t242;
	t301 = t243 * t245;
	t219 = t234 * t297 - t301;
	t231 = 0.1e1 / t235;
	t232 = 0.1e1 / t235 ^ 2;
	t238 = 0.1e1 / t242 ^ 2;
	t292 = qJD(4) * t245;
	t272 = t238 * t292;
	t237 = 0.1e1 / t242;
	t308 = t236 * t237;
	t276 = t234 * t308;
	t315 = t231 * t237;
	t333 = (t231 * t272 - t232 * t276) * t219 + t194 * t315;
	t312 = t235 * t242;
	t210 = atan2(-t216, t312);
	t207 = cos(t210);
	t206 = sin(t210);
	t321 = t206 * t216;
	t190 = t207 * t312 - t321;
	t187 = 0.1e1 / t190;
	t188 = 0.1e1 / t190 ^ 2;
	t332 = 0.2e1 * t219;
	t215 = t219 ^ 2;
	t185 = t188 * t215 + 0.1e1;
	t324 = t188 * t219;
	t214 = t216 ^ 2;
	t314 = t232 * t238;
	t211 = t214 * t314 + 0.1e1;
	t208 = 0.1e1 / t211;
	t313 = t234 * t236;
	t259 = t235 * t292 - t242 * t313;
	t278 = t216 * t314;
	t307 = t236 * t243;
	t275 = t235 * t307;
	t294 = qJD(1) * t246;
	t295 = qJD(1) * t243;
	t196 = t243 * t292 * t234 - t245 * t295 + (t294 * t234 + t275 - t291) * t242;
	t281 = t196 * t315;
	t178 = (t259 * t278 - t281) * t208;
	t257 = -t178 * t216 + t259;
	t173 = (-t178 * t312 - t196) * t206 + t257 * t207;
	t189 = t187 * t188;
	t330 = t173 * t189;
	t331 = (-t194 * t324 - t215 * t330) / t185 ^ 2;
	t255 = qJD(5) * t220 + t234 * t305 + t235 * t295;
	t195 = (-qJD(4) * t234 + qJD(1)) * t297 + (t274 + (-qJD(1) * t234 + qJD(4)) * t243) * t245;
	t290 = qJD(5) * t235;
	t266 = t246 * t290 + t195;
	t180 = t266 * t241 + t255 * t244;
	t198 = t204 ^ 2;
	t193 = t198 * t200 + 0.1e1;
	t323 = t200 * t204;
	t181 = -t255 * t241 + t266 * t244;
	t201 = t199 * t200;
	t326 = t181 * t201;
	t329 = (t180 * t323 - t198 * t326) / t193 ^ 2;
	t233 = t231 * t232;
	t239 = t237 * t238;
	t328 = (t196 * t278 + (-t232 * t239 * t292 + t233 * t238 * t313) * t214) / t211 ^ 2;
	t327 = t180 * t200;
	t325 = t188 * t194;
	t310 = t235 * t245;
	t261 = -t234 * t241 + t244 * t310;
	t213 = t261 * t246;
	t322 = t204 * t213;
	t320 = t206 * t219;
	t319 = t206 * t235;
	t318 = t207 * t216;
	t317 = t207 * t219;
	t316 = t207 * t234;
	t311 = t235 * t243;
	t306 = t236 * t245;
	t304 = t238 * t245;
	t298 = t246 * t187;
	t277 = t232 * t234 * t237;
	t264 = t216 * t277 + t243;
	t186 = t264 * t208;
	t296 = -t186 + t243;
	t288 = 0.2e1 * t331;
	t287 = 0.2e1 * t329;
	t286 = t189 * t332;
	t285 = t187 * t331;
	t284 = t231 * t328;
	t283 = t188 * t320;
	t280 = t204 * t326;
	t279 = t216 * t315;
	t269 = t188 * t288;
	t267 = t237 * t284;
	t197 = t220 * qJD(1) - t234 * t271 + t245 * t275 - t270;
	t265 = -t243 * t290 - t197;
	t218 = t234 * t301 - t297;
	t262 = t216 * t304 - t218 * t237;
	t260 = t234 * t244 + t241 * t310;
	t256 = -qJD(5) * t218 - t234 * t307 + t235 * t294;
	t212 = t260 * t246;
	t203 = -t218 * t244 - t241 * t311;
	t202 = -t218 * t241 + t244 * t311;
	t191 = 0.1e1 / t193;
	t183 = 0.1e1 / t185;
	t182 = t262 * t231 * t208;
	t177 = (-t206 + (t207 * t279 + t206) * t208) * t219;
	t176 = t186 * t318 + (-t296 * t319 - t316) * t242;
	t174 = t207 * t310 - t206 * t218 + (-t206 * t312 - t318) * t182;
	t172 = 0.2e1 * t264 * t328 + (-t196 * t277 - t294 + (-t231 * t308 + (t232 * t272 - 0.2e1 * t233 * t276) * t234) * t216) * t208;
	t170 = -0.2e1 * t262 * t284 + (t262 * t232 * t313 + (t196 * t304 - t197 * t237 + (t218 * t304 + (-0.2e1 * t239 * t245 ^ 2 - t237) * t216) * qJD(4)) * t231) * t208;
	t169 = (-t199 * t212 + t200 * t322) * t287 + (-t213 * t327 + (-t212 * t200 + 0.2e1 * t201 * t322) * t181 + (-t260 * t199 + t261 * t323) * t295 + (((-t241 * t293 + t334 * t244) * t199 - (-t334 * t241 - t244 * t293) * t323) * t235 + t263 * t234 * (-qJD(5) - t306)) * t246) * t191;
	t168 = t176 * t219 * t269 + (-(-t172 * t318 - (t178 * t321 - t196 * t207) * t186) * t324 + (t173 * t286 + t325) * t176 + (t235 * t298 - (t186 * t319 - t206 * t311 - t316) * t324) * t292) * t183 + (-0.2e1 * t285 * t309 + ((-t236 * t298 - (t296 * t236 + t178) * t283) * t234 + (-t187 * t295 + (-t246 * t173 - (-t172 - t294) * t320 - (-t296 * t178 - t236) * t317) * t188) * t235) * t183) * t242;
	t1 = [t333 * t208 + t267 * t332, t172, t172, t170, 0, 0; 0.2e1 * t216 * t285 + (-t196 * t187 + (t173 * t216 + t177 * t194) * t188) * t183 + (t177 * t269 + (0.2e1 * t177 * t330 + (t194 * t208 - t194 - (-t178 * t208 * t279 - 0.2e1 * t328) * t219) * t188 * t206 + (-(-0.2e1 * t216 * t267 - t178) * t324 + (-(t178 + t281) * t219 + t333 * t216) * t188 * t208) * t207) * t183) * t219, t168, t168, (t174 * t324 - t187 * t220) * t288 + (t174 * t325 + t195 * t187 + (t174 * t286 - t188 * t220) * t173 - (-t235 * t293 - t234 * t306 - t170 * t216 - t182 * t196 + (-t182 * t312 - t218) * t178) * t188 * t317 - (-t197 + (-t170 * t242 - t178 * t245) * t235 - t257 * t182) * t283) * t183, 0, 0; (-t199 * t202 + t203 * t323) * t287 + (0.2e1 * t203 * t280 + t265 * t303 + t256 * t199 * t244 + (t256 * t204 * t241 - t203 * t180 - t202 * t181 - t265 * t300) * t200) * t191, t169, t169, t263 * t219 * t287 + (t263 * t194 + ((-qJD(5) * t199 - 0.2e1 * t280) * t244 + (t180 * t244 + (-qJD(5) * t204 + t181) * t241) * t200) * t219) * t191, -0.2e1 * t329 + 0.2e1 * (t191 * t327 + (-t191 * t326 - t200 * t329) * t204) * t204, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:17
	% EndTime: 2020-06-20 08:27:23
	% DurationCPUTime: 5.76s
	% Computational Cost: add. (17109->206), mult. (24753->381), div. (1249->12), fcn. (30668->13), ass. (0->162)
	t352 = qJ(2) + qJ(3);
	t349 = sin(t352);
	t359 = cos(qJ(4));
	t360 = cos(qJ(1));
	t424 = t360 * t359;
	t355 = sin(qJ(4));
	t356 = sin(qJ(1));
	t427 = t356 * t355;
	t337 = t349 * t424 + t427;
	t421 = qJD(4) * t359;
	t401 = t360 * t421;
	t422 = qJD(4) * t355;
	t402 = t356 * t422;
	t350 = cos(t352);
	t351 = qJD(2) + qJD(3);
	t431 = t351 * t356;
	t408 = t350 * t431;
	t308 = t337 * qJD(1) - t349 * t402 + t359 * t408 - t401;
	t425 = t360 * t355;
	t426 = t356 * t359;
	t335 = t349 * t426 - t425;
	t354 = sin(qJ(5));
	t358 = cos(qJ(5));
	t318 = t356 * t350 * t354 + t335 * t358;
	t432 = t350 * t360;
	t460 = -qJD(1) * t432 + t349 * t431;
	t281 = t318 * qJD(5) + t308 * t354 + t460 * t358;
	t433 = t350 * t358;
	t380 = -t335 * t354 + t356 * t433;
	t315 = t380 ^ 2;
	t428 = t354 * t359;
	t435 = t349 * t358;
	t378 = t350 * t428 + t435;
	t329 = 0.1e1 / t378 ^ 2;
	t309 = t315 * t329 + 0.1e1;
	t439 = t380 * t329;
	t328 = 0.1e1 / t378;
	t395 = t351 * t359 + qJD(5);
	t386 = t395 * t354;
	t394 = qJD(5) * t359 + t351;
	t368 = -t354 * t422 + t394 * t358;
	t457 = t368 * t350;
	t289 = t349 * t386 - t457;
	t459 = t289 * t329;
	t447 = t328 * t459;
	t413 = 0.2e1 * (-t281 * t439 + t315 * t447) / t309 ^ 2;
	t461 = (t394 * t354 + t358 * t422) * t350;
	t377 = t349 * t427 + t424;
	t434 = t350 * t355;
	t376 = -t328 * t377 - t434 * t439;
	t458 = t354 * t376;
	t379 = -t349 * t354 + t359 * t433;
	t326 = t379 * t360;
	t406 = t349 * t425;
	t423 = qJD(1) * t356;
	t456 = (t355 * t423 - t401) * t350 - qJD(6) * t326 + t351 * t406;
	t282 = t380 * qJD(5) + t308 * t358 - t460 * t354;
	t310 = atan2(-t380, -t378);
	t301 = sin(t310);
	t302 = cos(t310);
	t278 = -t301 * t380 - t302 * t378;
	t275 = 0.1e1 / t278;
	t323 = t337 * t358 + t354 * t432;
	t336 = t406 - t426;
	t353 = sin(qJ(6));
	t357 = cos(qJ(6));
	t300 = -t323 * t357 + t336 * t353;
	t293 = 0.1e1 / t300;
	t276 = 0.1e1 / t278 ^ 2;
	t294 = 0.1e1 / t300 ^ 2;
	t321 = t337 * t354 - t358 * t432;
	t316 = t321 ^ 2;
	t274 = t316 * t276 + 0.1e1;
	t430 = t351 * t360;
	t367 = qJD(5) * t337 + t349 * t430 + t350 * t423;
	t392 = -qJD(1) * t349 + qJD(4);
	t393 = -qJD(4) * t349 + qJD(1);
	t306 = t393 * t425 + (t350 * t430 + t392 * t356) * t359;
	t391 = qJD(5) * t432 + t306;
	t279 = t391 * t354 + t367 * t358;
	t448 = t279 * t276;
	t303 = 0.1e1 / t309;
	t384 = -t281 * t328 + t289 * t439;
	t265 = t384 * t303;
	t389 = t301 * t378 - t302 * t380;
	t259 = t389 * t265 + t301 * t281 + t302 * t289;
	t277 = t275 * t276;
	t453 = t259 * t277;
	t454 = (-t316 * t453 + t321 * t448) / t274 ^ 2;
	t280 = -t367 * t354 + t391 * t358;
	t407 = t350 * t425;
	t305 = t377 * qJD(1) - t349 * t401 - t351 * t407 - t402;
	t418 = qJD(6) * t336;
	t267 = (qJD(6) * t323 - t305) * t357 + (t280 - t418) * t353;
	t298 = t323 * t353 + t336 * t357;
	t292 = t298 ^ 2;
	t286 = t292 * t294 + 0.1e1;
	t445 = t294 * t298;
	t419 = qJD(6) * t298;
	t268 = -t280 * t357 - t305 * t353 + t419;
	t450 = t268 * t293 * t294;
	t452 = (t267 * t445 - t292 * t450) / t286 ^ 2;
	t449 = t276 * t321;
	t446 = t293 * t357;
	t444 = t298 * t353;
	t443 = t298 * t357;
	t442 = t301 * t321;
	t441 = t302 * t321;
	t440 = t380 * t328;
	t437 = t336 * t354;
	t436 = t336 * t358;
	t429 = t353 * t293;
	t420 = qJD(5) * t358;
	t417 = -0.2e1 * t454;
	t416 = 0.2e1 * t454;
	t415 = -0.2e1 * t452;
	t414 = 0.2e1 * t452;
	t412 = -0.2e1 * t277 * t321;
	t411 = t276 * t442;
	t410 = t276 * t441;
	t400 = t259 * t412;
	t399 = t328 * t413;
	t398 = 0.2e1 * t380 * t447;
	t397 = -0.2e1 * t298 * t450;
	t390 = t358 * t418 - t306;
	t297 = t318 * t357 - t353 * t377;
	t388 = -t318 * t353 - t357 * t377;
	t387 = t349 * t395;
	t385 = qJD(6) * t407 + t379 * t423 - (-t358 * t387 - t461) * t360;
	t383 = t294 * t443 + t429;
	t382 = t318 * t328 + t379 * t439;
	t324 = t378 * t356;
	t331 = t349 * t428 - t433;
	t381 = t324 * t328 - t331 * t439;
	t375 = t349 * t351 * t355 - t350 * t421;
	t372 = t301 + (t302 * t440 - t301) * t303;
	t371 = qJD(1) * t378;
	t369 = -qJD(5) * t437 + qJD(6) * t337 - t305 * t358;
	t325 = t378 * t360;
	t314 = -t326 * t357 + t353 * t407;
	t313 = -t326 * t353 - t357 * t407;
	t312 = t337 * t353 + t357 * t436;
	t311 = -t337 * t357 + t353 * t436;
	t307 = t393 * t426 + (t392 * t360 - t408) * t355;
	t291 = t368 * t349 + t350 * t386;
	t290 = t395 * t435 + t461;
	t288 = t360 * t371 + (-t354 * t387 + t457) * t356;
	t284 = 0.1e1 / t286;
	t272 = 0.1e1 / t274;
	t271 = t303 * t458;
	t270 = t381 * t303;
	t269 = t382 * t303;
	t264 = t372 * t321;
	t262 = (-t301 * t377 + t302 * t434) * t354 - t389 * t271;
	t261 = -t389 * t270 + t301 * t324 + t302 * t331;
	t260 = -t389 * t269 + t301 * t318 - t302 * t379;
	t257 = t381 * t413 + (t331 * t398 - t288 * t328 + (-t281 * t331 - t289 * t324 + t291 * t380) * t329) * t303;
	t256 = t382 * t413 + (-t379 * t398 - t282 * t328 + (t281 * t379 - t289 * t318 + t290 * t380) * t329) * t303;
	t255 = t413 * t458 + (-t376 * t420 + (t398 * t434 - t307 * t328 + (-t281 * t434 + t289 * t377 - t375 * t380) * t329) * t354) * t303;
	t254 = (-t293 * t313 - t314 * t445) * t414 + (t314 * t397 + t385 * t429 + t456 * t446 + (t314 * t267 - t313 * t268 + t385 * t443 - t456 * t444) * t294) * t284;
	t253 = (-t261 * t449 + t275 * t325) * t416 + (t261 * t400 + (t325 * t259 + t261 * t279 + (-t257 * t380 - t270 * t281 + t291 + (-t270 * t378 + t324) * t265) * t441 + (t257 * t378 + t270 * t289 + t288 + (-t270 * t380 - t331) * t265) * t442) * t276 + (t289 * t360 + t356 * t371) * t275) * t272;
	t1 = [t321 * t399 + (-t279 * t328 - t321 * t459) * t303, t257, t257, t255, t256, 0; -t380 * t275 * t417 + (t281 * t275 + (t259 * t380 + t264 * t279) * t276) * t272 + (t264 * t276 * t417 + (-0.2e1 * t264 * t453 + (-t265 * t303 * t440 + t413) * t411 + (-t380 * t399 + t265 + (-t265 + t384) * t303) * t410 + t372 * t448) * t272) * t321, t253, t253, (-t262 * t449 - t275 * t437) * t416 + (t262 * t448 + (-t305 * t354 + t336 * t420) * t275 + (t262 * t412 - t276 * t437) * t259 + (-t377 * t420 + t255 * t378 + t271 * t289 + t307 * t354 + (-t271 * t380 - t354 * t434) * t265) * t411 + (t420 * t434 - t255 * t380 - (t265 * t378 + t281) * t271 + (-t265 * t377 - t375) * t354) * t410) * t272, (-t260 * t449 + t275 * t323) * t416 + (t260 * t400 - t280 * t275 + (t323 * t259 + t260 * t279 + (-t256 * t380 - t269 * t281 + t290 + (-t269 * t378 + t318) * t265) * t441 + (t256 * t378 + t269 * t289 + t282 + (-t269 * t380 + t379) * t265) * t442) * t276) * t272, 0; (t293 * t388 - t297 * t445) * t414 + ((t297 * qJD(6) + t282 * t353 - t307 * t357) * t293 + t297 * t397 + (t388 * t268 + (t388 * qJD(6) + t282 * t357 + t307 * t353) * t298 + t297 * t267) * t294) * t284, t254, t254, (-t293 * t311 - t312 * t445) * t414 + (t312 * t397 + t390 * t446 + t369 * t429 + (t312 * t267 - t311 * t268 + t369 * t443 - t390 * t444) * t294) * t284, t383 * t321 * t415 + (t383 * t279 + ((qJD(6) * t293 + t397) * t357 + (t267 * t357 + (-t268 - t419) * t353) * t294) * t321) * t284, t415 + 0.2e1 * (t267 * t294 * t284 + (-t284 * t450 - t294 * t452) * t298) * t298;];
	JaD_rot = t1;
end