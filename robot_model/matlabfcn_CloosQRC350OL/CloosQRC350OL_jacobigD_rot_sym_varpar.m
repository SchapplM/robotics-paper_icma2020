% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% CloosQRC350OL
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = CloosQRC350OL_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350OL_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobigD_rot_sym_varpar: pkin has to be [6x1] (double)');
JgD_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, -qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, -qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t32 = qJD(1) * sin(qJ(1));
	t31 = qJD(1) * cos(qJ(1));
	t1 = [0, -t31, -t31, 0, 0, 0; 0, -t32, -t32, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->9), mult. (13->8), div. (0->0), fcn. (13->4), ass. (0->9)
	t107 = qJD(2) + qJD(3);
	t108 = qJ(2) + qJ(3);
	t113 = sin(t108) * t107;
	t109 = sin(qJ(1));
	t112 = qJD(1) * t109;
	t110 = cos(qJ(1));
	t111 = qJD(1) * t110;
	t106 = cos(t108);
	t1 = [0, -t111, -t111, -t106 * t112 - t110 * t113, 0, 0; 0, -t112, -t112, t106 * t111 - t109 * t113, 0, 0; 0, 0, 0, -t107 * t106, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:51
	% EndTime: 2020-06-23 22:04:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->16), mult. (43->24), div. (0->0), fcn. (43->6), ass. (0->15)
	t156 = qJD(2) + qJD(3);
	t159 = sin(qJ(1));
	t167 = t156 * t159;
	t161 = cos(qJ(1));
	t166 = t156 * t161;
	t165 = qJD(1) * t159;
	t164 = qJD(1) * t161;
	t157 = qJ(2) + qJ(3);
	t154 = sin(t157);
	t163 = qJD(1) * t154 - qJD(4);
	t160 = cos(qJ(4));
	t162 = (qJD(4) * t154 - qJD(1)) * t160;
	t158 = sin(qJ(4));
	t155 = cos(t157);
	t1 = [0, -t164, -t164, -t154 * t166 - t155 * t165, t161 * t162 + (t155 * t166 - t163 * t159) * t158, 0; 0, -t165, -t165, -t154 * t167 + t155 * t164, t159 * t162 + (t155 * t167 + t163 * t161) * t158, 0; 0, 0, 0, -t156 * t155, t155 * qJD(4) * t160 - t156 * t154 * t158, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:51
	% EndTime: 2020-06-23 22:04:51
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (82->33), mult. (112->57), div. (0->0), fcn. (114->8), ass. (0->26)
	t231 = qJ(2) + qJ(3);
	t228 = sin(t231);
	t236 = cos(qJ(4));
	t252 = t228 * t236;
	t230 = qJD(2) + qJD(3);
	t251 = t230 * t228;
	t234 = sin(qJ(1));
	t250 = t230 * t234;
	t249 = t230 * t236;
	t232 = sin(qJ(5));
	t248 = t232 * t236;
	t237 = cos(qJ(1));
	t247 = t237 * t230;
	t246 = qJD(1) * t234;
	t245 = qJD(1) * t237;
	t233 = sin(qJ(4));
	t244 = qJD(4) * t233;
	t243 = qJD(5) + t249;
	t242 = qJD(4) * t228 - qJD(1);
	t241 = qJD(1) * t228 - qJD(4);
	t240 = t242 * t236;
	t239 = t241 * t234;
	t229 = cos(t231);
	t238 = -t228 * t247 - t229 * t246;
	t235 = cos(qJ(5));
	t1 = [0, -t245, -t245, t238, t237 * t240 + (t229 * t247 - t239) * t233, t239 * t248 + (-t243 * t229 + t242 * t233) * t232 * t237 + (-(t234 * t233 + t237 * t252) * qJD(5) + t238) * t235; 0, -t246, -t246, -t228 * t250 + t229 * t245, t234 * t240 + (t229 * t250 + t241 * t237) * t233, (-t241 * t248 + (qJD(1) * t229 + t233 * qJD(5)) * t235) * t237 + (-(qJD(1) * t233 - t228 * t244 + t229 * t249) * t232 - t235 * t251 + (-t229 * t232 - t235 * t252) * qJD(5)) * t234; 0, 0, 0, -t230 * t229, t229 * qJD(4) * t236 - t233 * t251, t243 * t232 * t228 + (t232 * t244 + (-qJD(5) * t236 - t230) * t235) * t229;];
	JgD_rot = t1;
end