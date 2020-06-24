% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = CloosQRC350DE_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobigD_rot_sym_varpar: pkin has to be [7x1] (double)');
JgD_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, -qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t29 = qJD(1) * cos(qJ(1));
	t1 = [0, t29, t29, 0, 0, 0; 0, -t31, -t31, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (13->7), mult. (13->8), div. (0->0), fcn. (13->4), ass. (0->9)
	t109 = qJD(2) + qJD(3);
	t110 = qJ(2) + qJ(3);
	t114 = sin(t110) * t109;
	t111 = sin(qJ(1));
	t113 = qJD(1) * t111;
	t112 = cos(qJ(1));
	t106 = qJD(1) * t112;
	t108 = cos(t110);
	t1 = [0, t106, t106, -t108 * t113 - t112 * t114, 0, 0; 0, -t113, -t113, -t108 * t106 + t111 * t114, 0, 0; 0, 0, 0, -t109 * t108, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:00
	% EndTime: 2020-06-23 21:15:00
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->14), mult. (43->24), div. (0->0), fcn. (43->6), ass. (0->15)
	t153 = qJ(2) + qJ(3);
	t150 = sin(t153);
	t156 = cos(qJ(4));
	t163 = t156 * (qJD(4) * t150 + qJD(1));
	t152 = qJD(2) + qJD(3);
	t155 = sin(qJ(1));
	t162 = t152 * t155;
	t157 = cos(qJ(1));
	t161 = t152 * t157;
	t160 = qJD(1) * t155;
	t149 = qJD(1) * t157;
	t158 = -qJD(1) * t150 - qJD(4);
	t154 = sin(qJ(4));
	t151 = cos(t153);
	t1 = [0, t149, t149, -t150 * t161 - t151 * t160, t157 * t163 + (t151 * t161 + t158 * t155) * t154, 0; 0, -t160, -t160, -t151 * t149 + t150 * t162, -t155 * t163 + (-t151 * t162 + t158 * t157) * t154, 0; 0, 0, 0, -t152 * t151, t151 * qJD(4) * t156 - t152 * t150 * t154, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:01
	% EndTime: 2020-06-23 21:15:01
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (130->34), mult. (184->60), div. (0->0), fcn. (185->8), ass. (0->27)
	t253 = qJ(2) + qJ(3);
	t250 = sin(t253);
	t252 = qJD(2) + qJD(3);
	t271 = t252 * t250;
	t256 = sin(qJ(1));
	t270 = t252 * t256;
	t258 = cos(qJ(4));
	t269 = t252 * t258;
	t259 = cos(qJ(1));
	t268 = t258 * t259;
	t267 = t259 * t252;
	t266 = qJD(1) * t256;
	t249 = qJD(1) * t259;
	t255 = sin(qJ(4));
	t265 = qJD(4) * t255;
	t264 = qJD(5) + t269;
	t263 = qJD(4) * t250 + qJD(1);
	t262 = -qJD(1) * t250 - qJD(4);
	t254 = sin(qJ(5));
	t261 = t262 * t254 * t258;
	t251 = cos(t253);
	t260 = -t250 * t267 - t251 * t266;
	t257 = cos(qJ(5));
	t248 = t264 * t254 * t250 + (t254 * t265 + (-qJD(5) * t258 - t252) * t257) * t251;
	t247 = (-t261 + (-qJD(1) * t251 + t255 * qJD(5)) * t257) * t259 + (-(qJD(1) * t255 + t250 * t265 - t251 * t269) * t254 + t257 * t271 + (t250 * t258 * t257 + t251 * t254) * qJD(5)) * t256;
	t246 = -t256 * t261 + (-t264 * t251 + t263 * t255) * t254 * t259 + (-(t250 * t268 - t256 * t255) * qJD(5) + t260) * t257;
	t1 = [0, t249, t249, t260, -pkin(7) * t246 + t263 * t268 + (t251 * t267 + t262 * t256) * t255, t246; 0, -t266, -t266, -t251 * t249 + t250 * t270, -pkin(7) * t247 - t263 * t258 * t256 + (-t251 * t270 + t262 * t259) * t255, t247; 0, 0, 0, -t252 * t251, qJD(4) * t251 * t258 - pkin(7) * t248 - t255 * t271, t248;];
	JgD_rot = t1;
end