% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = CloosQRC350DE_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobig_rot_sym_varpar: pkin has to be [7x1] (double)');
Jg_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0; 0, cos(qJ(1)), 0, 0, 0, 0; -1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t19 = cos(qJ(1));
	t18 = sin(qJ(1));
	t1 = [0, t18, t18, 0, 0, 0; 0, t19, t19, 0, 0, 0; -1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->3), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
	t82 = cos(qJ(1));
	t81 = sin(qJ(1));
	t80 = qJ(2) + qJ(3);
	t79 = cos(t80);
	t1 = [0, t81, t81, t82 * t79, 0, 0; 0, t82, t82, -t81 * t79, 0, 0; -1, 0, 0, -sin(t80), 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:15:00
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (10->5), mult. (9->8), div. (0->0), fcn. (21->6), ass. (0->9)
	t105 = qJ(2) + qJ(3);
	t103 = sin(t105);
	t106 = sin(qJ(4));
	t110 = t103 * t106;
	t109 = cos(qJ(1));
	t108 = cos(qJ(4));
	t107 = sin(qJ(1));
	t104 = cos(t105);
	t1 = [0, t107, t107, t109 * t104, t107 * t108 + t109 * t110, 0; 0, t109, t109, -t107 * t104, -t107 * t110 + t108 * t109, 0; -1, 0, 0, -t103, t104 * t106, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:00
	% EndTime: 2020-06-23 21:15:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->13), mult. (42->21), div. (0->0), fcn. (67->8), ass. (0->19)
	t170 = qJ(2) + qJ(3);
	t169 = cos(t170);
	t173 = sin(qJ(1));
	t182 = t173 * t169;
	t172 = sin(qJ(4));
	t181 = t173 * t172;
	t175 = cos(qJ(4));
	t180 = t173 * t175;
	t176 = cos(qJ(1));
	t179 = t176 * t169;
	t178 = t176 * t172;
	t177 = t176 * t175;
	t174 = cos(qJ(5));
	t171 = sin(qJ(5));
	t168 = sin(t170);
	t167 = -t169 * t175 * t171 - t168 * t174;
	t166 = -(t168 * t177 - t181) * t171 + t174 * t179;
	t165 = -(-t168 * t180 - t178) * t171 - t174 * t182;
	t1 = [0, t173, t173, t179, -pkin(7) * t166 + t168 * t178 + t180, t166; 0, t176, t176, -t182, -pkin(7) * t165 - t168 * t181 + t177, t165; -1, 0, 0, -t168, -pkin(7) * t167 + t169 * t172, t167;];
	Jg_rot = t1;
end