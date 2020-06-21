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
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
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
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	JaD_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (9->5), mult. (37->12), div. (15->6), fcn. (22->2), ass. (0->29)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = t1 ^ 2;
	t3 = cos(qJ(1));
	t4 = t3 ^ 2;
	t6 = t2 / t4;
	t7 = 0.1e1 + t6;
	t8 = t7 ^ 2;
	t11 = t1 / t3;
	t16 = t2 * t1 / t4 / t3;
	t19 = 0.2e1 / t8 * (t11 * qJD(1) + t16 * qJD(1));
	t21 = 0.1e1 / t7 * qJD(1);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = (-0.2e1 * t11 * t21 - 0.2e1 * t16 * t21 + t6 * t19 + t19);
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	JaD_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	JaD_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	JaD_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:15
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (3360->111), mult. (3810->319), div. (753->18), fcn. (4455->9), ass. (0->144)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = qJD(1) * t1;
	t3 = qJ(2) + qJ(3);
	t4 = cos(t3);
	t5 = sin(t3);
	t6 = 0.1e1 / t5;
	t8 = t1 ^ 2;
	t9 = t4 ^ 2;
	t11 = t5 ^ 2;
	t12 = 0.1e1 / t11;
	t14 = t8 * t9 * t12 + 0.1e1;
	t15 = 0.1e1 / t14;
	t18 = cos(qJ(1));
	t19 = qJD(2) + qJD(3);
	t20 = t18 * t19;
	t22 = t18 * t9;
	t23 = t12 * t15;
	t26 = t18 * t4;
	t27 = t14 ^ 2;
	t28 = 0.1e1 / t27;
	t30 = t1 * t9;
	t37 = t9 * t4;
	t40 = 0.1e1 / t11 / t5;
	t43 = t30 * t12 * qJD(1) * t18 - t8 * t37 * t40 * t19 - t8 * t4 * t6 * t19;
	t47 = qJD(1) * t18;
	t51 = t1 * t4;
	t52 = t6 * t15;
	t68 = t9 * qJD(1) * t18 * t12 * t15 - 0.2e1 * t37 * t1 * t40 * t15 * t19 - 0.2e1 * t30 * t12 * t28 * t43 - 0.2e1 * t1 * t28 * t43 - 0.2e1 * t51 * t52 * t19 + t47 * t15;
	t69 = atan2(t51, -t5);
	t70 = cos(t69);
	t71 = t70 * t5;
	t72 = sin(t69);
	t73 = t72 * t1;
	t74 = t73 * t4;
	t75 = -t71 + t74;
	t76 = 0.1e1 / t75;
	t78 = t18 ^ 2;
	t79 = t78 * t9;
	t80 = t75 ^ 2;
	t81 = 0.1e1 / t80;
	t83 = t79 * t81 + 0.1e1;
	t84 = 0.1e1 / t83;
	t85 = t4 * t76 * t84;
	t87 = t1 * t19;
	t89 = t5 * t76 * t84;
	t92 = t47 * t4;
	t93 = t87 * t5;
	t101 = -(t92 - t93) * t6 * t15 + t19 * t9 * t1 * t12 * t15;
	t102 = t101 * t72;
	t104 = t70 * t19;
	t106 = t101 * t70;
	t108 = t72 * qJD(1);
	t110 = t19 * t5;
	t112 = t102 * t5 - t104 * t4 + t106 * t51 + t108 * t26 - t73 * t110;
	t113 = t81 * t84 * t112;
	t115 = t83 ^ 2;
	t116 = 0.1e1 / t115;
	t126 = 0.1e1 / t80 / t75;
	t129 = -t78 * t4 * t81 * t19 * t5 - t22 * t81 * qJD(1) * t1 - t79 * t126 * t112;
	t130 = 0.2e1 * t76 * t116 * t129;
	t141 = t15 * t101;
	t155 = t15 * t70;
	t158 = t22 * t6;
	t170 = t72 * t18;
	t172 = t2 * t4 * t15 * t72 + t20 * t5 * t15 * t72 + 0.2e1 * t26 * t28 * t72 * t43 - t26 * t141 * t70 + qJD(1) * t8 * t9 * t52 * t70 + 0.2e1 * t26 * t15 * t70 * t1 * t19 + t18 * t37 * t12 * t155 * t87 + 0.2e1 * t158 * t28 * t70 * t1 * t43 + t158 * t141 * t73 - t79 * t6 * t155 * qJD(1) + t106 * t26 - t108 * t51 - t170 * t110;
	t175 = t4 * t81 * t84;
	t182 = -t158 * t155 * t1 - t26 * t15 * t72 + t170 * t4;
	t186 = t182 * t18;
	t189 = t5 * t81 * t84;
	t191 = t186 * t4;
	t193 = t126 * t84 * t112;
	t197 = 0.2e1 * t81 * t116 * t129;
	t202 = t18 * t5;
	t209 = t1 * t15 + t30 * t23;
	t210 = t209 * t101;
	t212 = t209 * t72;
	t213 = t19 * t4;
	t220 = t209 * t70;
	t223 = t1 * t5;
	t227 = t68 * t72 * t5 + t68 * t70 * t51 + t102 * t4 + t104 * t5 - t106 * t223 - t108 * t202 + t210 * t71 - t210 * t74 + t212 * t213 - t73 * t213 + t220 * t92 - t220 * t93;
	t234 = t212 * t5 + t220 * t51 - t70 * t4 - t73 * t5;
	t238 = t234 * t18;
	t241 = t238 * t4;
	t245 = t234 * qJD(1) * t1 * t175 - t227 * t18 * t175 + t238 * t19 * t189 + t202 * t113 + t202 * t130 + 0.2e1 * t241 * t193 + t241 * t197 + t2 * t89 - t20 * t85;
	t246 = sin(qJ(4));
	t247 = t5 * t246;
	t249 = t4 * t246;
	t251 = cos(qJ(4));
	t252 = qJD(4) * t251;
	t255 = t18 * qJD(4);
	t260 = -t1 * t246 + t202 * t251;
	t261 = 0.1e1 / t260;
	t265 = t1 * t251 + t202 * t246;
	t266 = t265 ^ 2;
	t267 = t260 ^ 2;
	t268 = 0.1e1 / t267;
	t270 = t266 * t268 + 0.1e1;
	t271 = 0.1e1 / t270;
	t275 = t18 * t251 - t223 * t246;
	t277 = t5 * t251;
	t279 = t4 * t251;
	t281 = qJD(4) * t246;
	t284 = t1 * qJD(4);
	t286 = -t2 * t277 + t20 * t279 - t202 * t281 - t47 * t246 - t284 * t251;
	t287 = t271 * t286;
	t290 = t270 ^ 2;
	t291 = 0.1e1 / t290;
	t292 = t265 * t268;
	t298 = -t2 * t247 + t20 * t249 + t202 * t252 - t284 * t246 + t47 * t251;
	t301 = 0.1e1 / t267 / t260;
	t304 = -t266 * t301 * t286 + t292 * t298;
	t305 = 0.2e1 * t291 * t304;
	t314 = t268 * t271;
	t318 = -t18 * t246 - t223 * t251;
	t321 = t318 * t265;
	t323 = t301 * t271 * t286;
	t327 = 0.2e1 * t268 * t291 * t304;
	t330 = t2 * t4;
	t332 = t246 * t261 * t271;
	t334 = t20 * t5;
	t336 = t26 * qJD(4);
	t340 = t26 * t246;
	t347 = t251 * t265 * t314;
	t353 = t26 * t251;
	t363 = t336 * t246 * t265 * t314 + t336 * t251 * t261 * t271 - 0.2e1 * t340 * t261 * t291 * t304 + 0.2e1 * t353 * t265 * t301 * t287 - t353 * t298 * t268 * t271 - t340 * t314 * t286 + t353 * t292 * t305 - t330 * t332 + t330 * t347 - t334 * t332 + t334 * t347;
	t366 = -t265 ^ 2;
	unknown(1,1) = t2 * t4 * t6 * t15 + 0.2e1 * t26 * t6 * t28 * t43 + t22 * t23 * t19 + t20 * t15;
	unknown(1,2) = t68;
	unknown(1,3) = t68;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = t182 * qJD(1) * t1 * t175 - t172 * t18 * t175 + t186 * t19 * t189 + t51 * t113 + t51 * t130 + 0.2e1 * t191 * t193 + t191 * t197 - t47 * t85 + t87 * t89;
	unknown(2,2) = t245;
	unknown(2,3) = t245;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = (-t2 * t251 - t223 * t252 - t255 * t246 - t47 * t247 - t87 * t249) * t261 * t271 - t275 * t268 * t287 - t275 * t261 * t305 - (t2 * t246 + t223 * t281 - t255 * t251 - t47 * t277 - t87 * t279) * t265 * t314 - t318 * t298 * t314 + 0.2e1 * t321 * t323 + t321 * t327;
	unknown(3,2) = t363;
	unknown(3,3) = t363;
	unknown(3,4) = 0.2e1 * t298 * t265 * t314 + 0.2e1 * t366 * t323 + t366 * t327 - t305;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	JaD_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:15
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (6566->182), mult. (9976->458), div. (1522->22), fcn. (12462->11), ass. (0->184)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = qJD(1) * t1;
	t3 = qJ(2) + qJ(3);
	t4 = sin(t3);
	t5 = sin(qJ(4));
	t6 = t4 * t5;
	t8 = cos(qJ(1));
	t9 = qJD(2) + qJD(3);
	t10 = t8 * t9;
	t11 = cos(t3);
	t12 = t11 * t5;
	t14 = t8 * t4;
	t15 = cos(qJ(4));
	t16 = qJD(4) * t15;
	t18 = qJD(1) * t8;
	t20 = t1 * qJD(4);
	t22 = t10 * t12 + t14 * t16 + t18 * t15 - t2 * t6 - t20 * t5;
	t23 = 0.1e1 / t11;
	t24 = t22 * t23;
	t25 = 0.1e1 / t5;
	t26 = t1 * t4;
	t29 = -t8 * t15 + t26 * t5;
	t30 = t29 ^ 2;
	t31 = t11 ^ 2;
	t32 = 0.1e1 / t31;
	t33 = t30 * t32;
	t34 = t5 ^ 2;
	t35 = 0.1e1 / t34;
	t37 = t33 * t35 + 0.1e1;
	t38 = 0.1e1 / t37;
	t39 = t25 * t38;
	t43 = t1 * t15 + t14 * t5;
	t44 = t43 * t32;
	t47 = t38 * t9 * t4;
	t49 = t43 * t23;
	t52 = t38 * qJD(4) * t15;
	t54 = t37 ^ 2;
	t55 = 0.1e1 / t54;
	t57 = t29 * t32;
	t59 = t1 * t9;
	t63 = t8 * qJD(4);
	t65 = t59 * t12 + t2 * t15 + t26 * t16 + t18 * t6 + t63 * t5;
	t66 = t35 * t65;
	t69 = 0.1e1 / t31 / t11;
	t75 = 0.1e1 / t34 / t5;
	t79 = t30 * t69 * t35 * t9 * t4 - t33 * t75 * qJD(4) * t15 + t57 * t66;
	t80 = 0.2e1 * t25 * t55 * t79;
	t88 = t25 * t29 * t38;
	t95 = t4 * t25;
	t99 = t4 ^ 2;
	t110 = -t4 * t35 * t29 * t32 * t38 * t16 + 0.2e1 * t99 * t25 * t29 * t69 * t38 * t9 - 0.2e1 * t95 * t29 * t32 * t55 * t79 + t95 * t65 * t32 * t38 - 0.2e1 * t1 * t55 * t79 + t9 * t23 * t88 + t18 * t38;
	t111 = t4 * t15;
	t113 = t11 * t15;
	t115 = qJD(4) * t5;
	t119 = t18 * t111 + t59 * t113 - t26 * t115 + t63 * t15 - t2 * t5;
	t124 = t26 * t15 + t8 * t5;
	t128 = t124 * t23;
	t134 = t35 * t38;
	t135 = t9 * t4;
	t140 = t23 * t15;
	t143 = t15 ^ 2;
	t154 = 0.2e1 * t23 * t143 * t29 * t75 * t38 * qJD(4) - t32 * t15 * t29 * t134 * t135 + 0.2e1 * t140 * t29 * t35 * t55 * t79 + t124 * t32 * t25 * t47 + t23 * qJD(4) * t88 + t119 * t23 * t39 - t128 * t35 * t52 - t140 * t66 * t38 - t128 * t80;
	t155 = atan2(t29, t12);
	t156 = cos(t155);
	t157 = t156 * t11;
	t158 = t157 * t5;
	t159 = sin(t155);
	t160 = -t159 * t29;
	t161 = t158 - t160;
	t162 = 0.1e1 / t161;
	t164 = t43 ^ 2;
	t165 = t161 ^ 2;
	t166 = 0.1e1 / t165;
	t168 = t164 * t166 + 0.1e1;
	t169 = 0.1e1 / t168;
	t174 = t135 * t5;
	t176 = t11 * qJD(4) * t15;
	t182 = t65 * t23 * t39 - (-t174 + t176) * t29 * t32 * t35 * t38;
	t183 = t182 * t159;
	t185 = t156 * t9;
	t188 = t182 * t156;
	t191 = -t183 * t12 + t157 * t16 + t159 * t65 - t185 * t6 + t188 * t29;
	t192 = t169 * t191;
	t195 = t168 ^ 2;
	t196 = 0.1e1 / t195;
	t200 = 0.1e1 / t165 / t161;
	t203 = -t164 * t200 * t191 + t43 * t166 * t22;
	t204 = 0.2e1 * t196 * t203;
	t211 = t43 * t38;
	t214 = t38 * t156;
	t215 = -t214 * t29;
	t218 = -t156 * t29;
	t224 = t49 * t25;
	t236 = -0.2e1 * t224 * t55 * t156 * t29 * t79 + t49 * t134 * t218 * t16 - t44 * t39 * t218 * t135 + 0.2e1 * t43 * t55 * t159 * t79 + t224 * t38 * t182 * t160 - t22 * t38 * t159 + t224 * t214 * t65 - t24 * t25 * t215 + t159 * t22 - t211 * t188 + t188 * t43;
	t238 = t166 * t169;
	t243 = -t211 * t159 + t159 * t43 - t224 * t215;
	t246 = t243 * t43;
	t248 = t200 * t169 * t191;
	t252 = 0.2e1 * t166 * t196 * t203;
	t257 = t5 * t162 * t169;
	t261 = t8 * t11;
	t266 = t261 * t5;
	t277 = t95 * t57 * t38 + t1 * t38;
	t278 = t277 * t182;
	t280 = t277 * t159;
	t285 = t156 * t4;
	t290 = t277 * t156;
	t292 = t1 * t11;
	t297 = t159 * t1;
	t300 = t159 * qJD(1) * t266 - t110 * t159 * t12 + t110 * t156 * t29 + t188 * t292 * t5 - t185 * t12 - t278 * t158 - t285 * t16 + t278 * t160 + t280 * t174 - t297 * t174 - t280 * t176 + t297 * t176 + t183 * t6 + t290 * t65;
	t307 = -t280 * t12 + t297 * t12 - t285 * t5 + t290 * t29;
	t310 = t307 * t43;
	t314 = t261 * qJD(4) * t15 * t162 * t169 - 0.2e1 * t266 * t162 * t196 * t203 - t10 * t4 * t257 - t2 * t11 * t257 - t266 * t238 * t191 - t307 * t22 * t238 - t300 * t43 * t238 + 0.2e1 * t310 * t248 + t310 * t252;
	t320 = t10 * t113 - t2 * t111 - t14 * t115 - t20 * t15 - t18 * t5;
	t325 = -t1 * t5 + t14 * t15;
	t336 = -t140 * t29 * t35 * t38 + t128 * t39;
	t337 = t336 * t182;
	t339 = t336 * t159;
	t348 = t336 * t156;
	t352 = -t154 * t159 * t12 + t154 * t156 * t29 - t185 * t111 - t183 * t113 - t157 * t115 + t159 * t119 + t188 * t124 - t337 * t158 + t337 * t160 + t339 * t174 - t339 * t176 + t348 * t65;
	t359 = -t339 * t12 + t159 * t124 + t157 * t15 + t348 * t29;
	t362 = t359 * t43;
	t367 = sin(qJ(5));
	t369 = -t124 * qJD(5);
	t370 = cos(qJ(5));
	t372 = t11 * t370;
	t374 = t4 * t370;
	t376 = qJD(5) * t367;
	t381 = t261 * t367 + t325 * t370;
	t382 = 0.1e1 / t381;
	t386 = -t261 * t370 + t325 * t367;
	t387 = t386 ^ 2;
	t388 = t381 ^ 2;
	t389 = 0.1e1 / t388;
	t391 = t387 * t389 + 0.1e1;
	t392 = 0.1e1 / t391;
	t396 = -t124 * t367 + t292 * t370;
	t399 = t325 * qJD(5);
	t401 = t11 * t367;
	t403 = t4 * t367;
	t405 = qJD(5) * t370;
	t407 = -t10 * t403 - t2 * t401 + t261 * t405 + t320 * t370 - t399 * t367;
	t408 = t392 * t407;
	t411 = t391 ^ 2;
	t412 = 0.1e1 / t411;
	t413 = t386 * t389;
	t419 = t10 * t374 + t2 * t372 + t261 * t376 + t320 * t367 + t399 * t370;
	t422 = 0.1e1 / t388 / t381;
	t425 = -t387 * t422 * t407 + t413 * t419;
	t426 = 0.2e1 * t412 * t425;
	t435 = t389 * t392;
	t439 = -t124 * t370 - t292 * t367;
	t442 = t439 * t386;
	t444 = t422 * t392 * t407;
	t448 = 0.2e1 * t389 * t412 * t425;
	t457 = t15 * qJD(5);
	t469 = t261 * t15 * t367 + t14 * t370;
	t491 = t261 * t15 * t370 - t14 * t367;
	t494 = t491 * t386;
	t498 = (-t10 * t111 * t367 - t2 * t113 * t367 - t261 * t115 * t367 + t261 * t457 * t370 + t10 * t372 - t14 * t376 - t2 * t374) * t382 * t392 - t469 * t389 * t408 - t469 * t382 * t426 - (-t10 * t111 * t370 - t2 * t113 * t370 - t261 * t115 * t370 - t261 * t457 * t367 - t10 * t401 - t14 * t405 + t2 * t403) * t386 * t435 - t491 * t419 * t435 + 0.2e1 * t494 * t444 + t494 * t448;
	t502 = -t43 * qJD(5);
	t506 = -t43 * t367;
	t513 = t413 * t392;
	t517 = -t43 * t370;
	t521 = t517 * t386;
	t528 = -t386 ^ 2;
	unknown(1,1) = t44 * t25 * t47 - t49 * t35 * t52 + t24 * t39 - t49 * t80;
	unknown(1,2) = t110;
	unknown(1,3) = t110;
	unknown(1,4) = t154;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t65 * t162 * t169 + t29 * t162 * t204 + t29 * t166 * t192 - t243 * t22 * t238 - t236 * t43 * t238 + 0.2e1 * t246 * t248 + t246 * t252;
	unknown(2,2) = t314;
	unknown(2,3) = t314;
	unknown(2,4) = t320 * t162 * t169 - t325 * t162 * t204 - t325 * t166 * t192 - t359 * t22 * t238 - t352 * t43 * t238 + 0.2e1 * t362 * t248 + t362 * t252;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = (-t119 * t367 + t18 * t372 - t292 * t376 + t369 * t370 - t59 * t374) * t382 * t392 - t396 * t389 * t408 - t396 * t382 * t426 - (-t119 * t370 - t18 * t401 - t292 * t405 - t369 * t367 + t59 * t403) * t386 * t435 - t439 * t419 * t435 + 0.2e1 * t442 * t444 + t442 * t448;
	unknown(3,2) = t498;
	unknown(3,3) = t498;
	unknown(3,4) = -t22 * t367 * t382 * t392 + t502 * t370 * t382 * t392 - 0.2e1 * t506 * t382 * t412 * t425 - t517 * t419 * t389 * t392 + t22 * t370 * t513 + t502 * t367 * t513 - t506 * t435 * t407 + 0.2e1 * t521 * t444 + t521 * t448;
	unknown(3,5) = 0.2e1 * t419 * t386 * t435 + 0.2e1 * t528 * t444 + t528 * t448 - t426;
	unknown(3,6) = 0.0e0;
	JaD_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:15
	% EndTime: 2020-06-19 21:40:16
	% DurationCPUTime: 1.39s
	% Computational Cost: add. (17929->289), mult. (25614->606), div. (1243->18), fcn. (30756->13), ass. (0->236)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = qJD(1) * t1;
	t3 = qJ(2) + qJ(3);
	t4 = sin(t3);
	t5 = cos(qJ(4));
	t6 = t4 * t5;
	t8 = cos(qJ(1));
	t9 = qJD(2) + qJD(3);
	t10 = t8 * t9;
	t11 = cos(t3);
	t12 = t11 * t5;
	t14 = t8 * t4;
	t15 = sin(qJ(4));
	t16 = qJD(4) * t15;
	t18 = qJD(1) * t8;
	t20 = t1 * qJD(4);
	t22 = -t10 * t12 + t14 * t16 + t18 * t15 + t2 * t6 + t20 * t5;
	t23 = sin(qJ(5));
	t27 = t1 * t15 - t14 * t5;
	t29 = cos(qJ(5));
	t31 = t11 * t29;
	t32 = t2 * t31;
	t33 = t4 * t29;
	t34 = t10 * t33;
	t35 = t8 * t11;
	t36 = qJD(5) * t23;
	t37 = t35 * t36;
	t38 = t27 * qJD(5) * t29 + t22 * t23 - t32 - t34 - t37;
	t39 = t12 * t23;
	t40 = -t39 - t33;
	t41 = 0.1e1 / t40;
	t42 = t38 * t41;
	t43 = t1 * t4;
	t46 = -t8 * t15 - t43 * t5;
	t48 = t1 * t11;
	t50 = t46 * t23 + t48 * t29;
	t51 = t50 ^ 2;
	t52 = t40 ^ 2;
	t53 = 0.1e1 / t52;
	t55 = t51 * t53 + 0.1e1;
	t56 = 0.1e1 / t55;
	t59 = t35 * t29;
	t60 = t27 * t23 + t59;
	t61 = t60 * t53;
	t62 = t9 * t4;
	t63 = t5 * t23;
	t65 = t11 * qJD(4);
	t66 = t15 * t23;
	t68 = qJD(5) * t29;
	t70 = t9 * t11;
	t72 = t4 * qJD(5);
	t74 = -t12 * t68 + t72 * t23 - t70 * t29 + t62 * t63 + t65 * t66;
	t75 = t56 * t74;
	t77 = t60 * t41;
	t78 = t55 ^ 2;
	t79 = 0.1e1 / t78;
	t80 = t50 * t53;
	t82 = t1 * t9;
	t86 = t8 * qJD(4);
	t88 = -t82 * t12 + t2 * t15 + t43 * t16 - t18 * t6 - t86 * t5;
	t90 = t46 * qJD(5);
	t95 = t18 * t31 + t88 * t23 + t90 * t29 - t82 * t33 - t48 * t36;
	t98 = 0.1e1 / t52 / t40;
	t101 = -t51 * t98 * t74 + t80 * t95;
	t102 = 0.2e1 * t79 * t101;
	t106 = t6 * t23;
	t108 = t16 * t23;
	t110 = t5 * qJD(5);
	t111 = t110 * t29;
	t116 = t82 * t106 + t48 * t108 - t48 * t111 - t18 * t33 - t18 * t39 - t82 * t31 + t43 * t36;
	t121 = -t43 * t29 - t48 * t63;
	t124 = t121 * t41;
	t133 = -t4 * qJD(4) * t66 + t11 * qJD(5) * t23 + t62 * t29 + t6 * t68 + t70 * t63;
	t135 = t53 * t56;
	t137 = t106 - t31;
	t140 = t137 * t50;
	t142 = t98 * t56 * t74;
	t146 = 0.2e1 * t53 * t79 * t101;
	t148 = t116 * t41 * t56 - t121 * t53 * t75 - t133 * t50 * t135 - t137 * t95 * t135 - t124 * t102 + 0.2e1 * t140 * t142 + t140 * t146;
	t149 = t4 * t15;
	t151 = t11 * t15;
	t153 = qJD(4) * t5;
	t157 = t18 * t149 + t86 * t15 + t82 * t151 + t43 * t153 + t2 * t5;
	t159 = t41 * t56;
	t163 = t43 * t15 - t8 * t5;
	t168 = t163 * t23;
	t176 = t23 * t50 * t135;
	t184 = t151 * t23;
	t194 = -t151 * qJD(5) * t29 * t50 * t135 + t163 * qJD(5) * t29 * t41 * t56 - 0.2e1 * t168 * t41 * t79 * t101 + 0.2e1 * t184 * t50 * t98 * t75 - t184 * t95 * t53 * t56 + t184 * t80 * t102 - t168 * t135 * t74 + t62 * t15 * t176 + t157 * t23 * t159 - t65 * t5 * t176;
	t197 = t11 * t23;
	t199 = t4 * t23;
	t202 = -t18 * t197 + t82 * t199 - t90 * t23 + t88 * t29 - t48 * t68;
	t207 = -t48 * t23 + t46 * t29;
	t210 = t207 * t41;
	t212 = t5 * t29;
	t219 = t65 * t15 * t29 + t12 * t36 + t62 * t212 + t70 * t23 + t72 * t29;
	t222 = t12 * t29;
	t223 = -t222 + t199;
	t226 = t223 * t50;
	t230 = -t219 * t50 * t135 - t223 * t95 * t135 + t202 * t41 * t56 - t207 * t53 * t75 - t210 * t102 + 0.2e1 * t226 * t142 + t226 * t146;
	t231 = atan2(t50, t40);
	t232 = cos(t231);
	t233 = t232 * t40;
	t234 = sin(t231);
	t235 = -t234 * t50;
	t236 = t233 - t235;
	t237 = 0.1e1 / t236;
	t240 = t27 * t23 + t59;
	t241 = t240 ^ 2;
	t242 = t236 ^ 2;
	t243 = 0.1e1 / t242;
	t245 = t241 * t243 + 0.1e1;
	t246 = 0.1e1 / t245;
	t253 = -t74 * t50 * t135 + t95 * t41 * t56;
	t254 = t253 * t234;
	t257 = t253 * t232;
	t260 = t232 * t74 + t234 * t95 - t254 * t40 + t257 * t50;
	t261 = t246 * t260;
	t264 = t245 ^ 2;
	t265 = 0.1e1 / t264;
	t268 = -t27 * qJD(5);
	t270 = t22 * t23 - t268 * t29 - t32 - t34 - t37;
	t273 = 0.1e1 / t242 / t236;
	t276 = t240 * t243 * t270 - t241 * t273 * t260;
	t277 = 0.2e1 * t265 * t276;
	t284 = t60 * t56;
	t286 = t56 * t232;
	t287 = -t286 * t50;
	t290 = -t232 * t50;
	t305 = t243 * t246;
	t310 = -t284 * t234 + t234 * t60 - t77 * t287;
	t313 = t310 * t240;
	t315 = t273 * t246 * t260;
	t319 = 0.2e1 * t243 * t265 * t276;
	t334 = -t14 * t29 - t35 * t63;
	t343 = t124 * t56 - t140 * t135;
	t344 = t343 * t253;
	t346 = t343 * t234;
	t353 = t343 * t232;
	t364 = t234 * t121 + t232 * t137 - t346 * t40 + t353 * t50;
	t367 = t364 * t240;
	t371 = (t10 * t106 - t10 * t31 + t35 * t108 - t35 * t111 + t14 * t36 + t2 * t33 + t2 * t39) * t237 * t246 - t334 * t243 * t261 - t334 * t237 * t277 - (t148 * t232 * t50 - t148 * t234 * t40 + t234 * t116 + t257 * t121 + t232 * t133 - t254 * t137 - t344 * t233 + t344 * t235 - t346 * t74 + t353 * t95) * t240 * t305 - t364 * t270 * t305 + 0.2e1 * t367 * t315 + t367 * t319;
	t377 = -t10 * t151 - t14 * t153 + t2 * t149 + t20 * t15 - t18 * t5;
	t383 = -t1 * t5 - t14 * t15;
	t384 = t383 * qJD(5);
	t388 = t383 * t23;
	t399 = -t184 * t80 * t56 + t168 * t159;
	t400 = t399 * t253;
	t402 = t399 * t234;
	t408 = t232 * t11;
	t417 = t399 * t232;
	t422 = t234 * t163;
	t424 = t408 * t15 * qJD(5) * t29 - t232 * t9 * t149 * t23 + t408 * t153 * t23 + t234 * t157 * t23 + t194 * t232 * t50 - t194 * t234 * t40 + t257 * t168 - t254 * t184 - t400 * t233 + t400 * t235 - t402 * t74 + t417 * t95 + t422 * t68;
	t431 = t422 * t23 - t402 * t40 + t408 * t66 + t417 * t50;
	t434 = t431 * t240;
	t444 = t10 * t199 + t2 * t197 + t22 * t29 + t268 * t23 - t35 * t68;
	t449 = -t35 * t23 + t27 * t29;
	t458 = -t226 * t135 + t210 * t56;
	t459 = t458 * t253;
	t461 = t458 * t234;
	t468 = t458 * t232;
	t479 = t234 * t207 + t232 * t223 - t461 * t40 + t468 * t50;
	t482 = t479 * t240;
	t488 = pkin(7) * qJ(5) - qJ(6);
	t489 = sin(t488);
	t492 = pkin(7) * qJD(5) - qJD(6);
	t493 = t207 * t492;
	t494 = cos(t488);
	t497 = -t163 * t492;
	t502 = t383 * t489 + t449 * t494;
	t503 = 0.1e1 / t502;
	t507 = t383 * t494 - t449 * t489;
	t508 = t507 ^ 2;
	t509 = t502 ^ 2;
	t510 = 0.1e1 / t509;
	t512 = t508 * t510 + 0.1e1;
	t513 = 0.1e1 / t512;
	t517 = t163 * t494 + t207 * t489;
	t520 = -t449 * t492;
	t523 = -t383 * t492;
	t525 = t377 * t489 + t444 * t494 + t520 * t489 - t523 * t494;
	t526 = t513 * t525;
	t529 = t512 ^ 2;
	t530 = 0.1e1 / t529;
	t536 = t377 * t494 - t444 * t489 + t523 * t489 + t520 * t494;
	t539 = 0.1e1 / t509 / t502;
	t542 = t507 * t510 * t536 - t508 * t539 * t525;
	t543 = 0.2e1 * t530 * t542;
	t551 = t510 * t513;
	t555 = t163 * t489 - t207 * t494;
	t558 = t555 * t507;
	t560 = t539 * t513 * t525;
	t564 = 0.2e1 * t510 * t530 * t542;
	t577 = -t10 * t6 * t29 - t35 * t110 * t23 - t35 * t16 * t29 - t10 * t197 - t14 * t68 + t2 * t199 - t2 * t222;
	t581 = -t14 * t23 + t35 * t212;
	t582 = t581 * t492;
	t590 = t15 * t492;
	t599 = -t35 * t15 * t494 + t581 * t489;
	t620 = -t35 * t15 * t489 - t581 * t494;
	t623 = t620 * t507;
	t627 = (t10 * t149 * t494 + t2 * t151 * t494 - t35 * t153 * t494 + t35 * t590 * t489 + t577 * t489 + t582 * t494) * t503 * t513 - t599 * t510 * t526 - t599 * t503 * t543 - (t10 * t149 * t489 + t2 * t151 * t489 - t35 * t153 * t489 - t35 * t590 * t494 + t582 * t489 - t577 * t494) * t507 * t551 - t620 * t536 * t551 + 0.2e1 * t623 * t560 + t623 * t564;
	t628 = t377 * t29;
	t632 = t383 * t29;
	t633 = t492 * t494;
	t636 = -t27 * t492;
	t643 = t27 * t494 + t632 * t489;
	t651 = t492 * t489;
	t660 = t27 * t489 - t632 * t494;
	t663 = t660 * t507;
	t669 = t240 * t492;
	t671 = -t444 * pkin(7);
	t673 = -t449 * pkin(7);
	t675 = -t377 * pkin(7);
	t677 = -t383 * pkin(7);
	t685 = t240 * t489 + t677 * t489 + t673 * t494;
	t702 = -t240 * t494 + t673 * t489 - t677 * t494;
	t705 = t702 * t507;
	t712 = -t507 ^ 2;
	unknown(1,1) = -t77 * t102 + t42 * t56 - t61 * t75;
	unknown(1,2) = t148;
	unknown(1,3) = t148;
	unknown(1,4) = t194;
	unknown(1,5) = t230;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t95 * t237 * t246 + t50 * t243 * t261 + t50 * t237 * t277 - (0.2e1 * t60 * t79 * t234 * t101 + 0.2e1 * t77 * t79 * t290 * t101 - t77 * t56 * t254 * t50 + t61 * t56 * t290 * t74 - t38 * t56 * t234 + t77 * t286 * t95 + t234 * t38 - t284 * t257 + t257 * t60 - t42 * t287) * t240 * t305 - t310 * t270 * t305 + 0.2e1 * t313 * t315 + t313 * t319;
	unknown(2,2) = t371;
	unknown(2,3) = t371;
	unknown(2,4) = -t377 * t23 * t237 * t246 - t384 * t29 * t237 * t246 + 0.2e1 * t388 * t237 * t265 * t276 - t424 * t240 * t305 + t388 * t305 * t260 - t431 * t270 * t305 + 0.2e1 * t434 * t315 + t434 * t319;
	unknown(2,5) = t444 * t237 * t246 - t449 * t243 * t261 - t449 * t237 * t277 - (t230 * t232 * t50 - t230 * t234 * t40 + t234 * t202 + t257 * t207 + t232 * t219 - t254 * t223 - t459 * t233 + t459 * t235 - t461 * t74 + t468 * t95) * t240 * t305 - t479 * t270 * t305 + 0.2e1 * t482 * t315 + t482 * t319;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = (t157 * t494 + t202 * t489 + t497 * t489 + t493 * t494) * t503 * t513 - t517 * t510 * t526 - t517 * t503 * t543 - (t157 * t489 - t202 * t494 + t493 * t489 - t497 * t494) * t507 * t551 - t555 * t536 * t551 + 0.2e1 * t558 * t560 + t558 * t564;
	unknown(3,2) = t627;
	unknown(3,3) = t627;
	unknown(3,4) = (-t384 * t23 * t489 + t22 * t494 + t628 * t489 + t636 * t489 + t632 * t633) * t503 * t513 - t643 * t510 * t526 - t643 * t503 * t543 - (t384 * t23 * t494 + t22 * t489 - t628 * t494 - t636 * t494 + t632 * t651) * t507 * t551 - t660 * t536 * t551 + 0.2e1 * t663 * t560 + t663 * t564;
	unknown(3,5) = (t270 * t489 + t675 * t489 + t669 * t494 + t671 * t494 + t677 * t633 - t673 * t651) * t503 * t513 - t685 * t510 * t526 - t685 * t503 * t543 - (-t270 * t494 + t669 * t489 + t671 * t489 - t675 * t494 + t673 * t633 + t677 * t651) * t507 * t551 - t702 * t536 * t551 + 0.2e1 * t705 * t560 + t705 * t564;
	unknown(3,6) = 0.2e1 * t536 * t507 * t551 + 0.2e1 * t712 * t560 + t712 * t564 - t543;
	JaD_rot = unknown;
end