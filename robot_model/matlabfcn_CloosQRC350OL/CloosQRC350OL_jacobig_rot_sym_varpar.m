% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = CloosQRC350OL_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350OL_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobig_rot_sym_varpar: pkin has to be [6x1] (double)');
Jg_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, -sin(qJ(1)), 0, 0, 0, 0; 0, cos(qJ(1)), 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t19 = cos(qJ(1));
	t18 = sin(qJ(1));
	t1 = [0, -t18, -t18, 0, 0, 0; 0, t19, t19, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->4), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->5)
	t81 = cos(qJ(1));
	t80 = sin(qJ(1));
	t79 = qJ(2) + qJ(3);
	t78 = cos(t79);
	t1 = [0, -t80, -t80, t81 * t78, 0, 0; 0, t81, t81, t80 * t78, 0, 0; 1, 0, 0, -sin(t79), 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->6), mult. (9->8), div. (0->0), fcn. (21->6), ass. (0->9)
	t109 = qJ(2) + qJ(3);
	t107 = sin(t109);
	t110 = sin(qJ(4));
	t114 = t107 * t110;
	t113 = cos(qJ(1));
	t112 = cos(qJ(4));
	t111 = sin(qJ(1));
	t108 = cos(t109);
	t1 = [0, -t111, -t111, t113 * t108, -t111 * t112 + t113 * t114, 0; 0, t113, t113, t111 * t108, t111 * t114 + t113 * t112, 0; 1, 0, 0, -t107, t108 * t110, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (22->11), mult. (24->18), div. (0->0), fcn. (44->8), ass. (0->16)
	t157 = qJ(2) + qJ(3);
	t156 = cos(t157);
	t160 = sin(qJ(1));
	t169 = t160 * t156;
	t159 = sin(qJ(4));
	t168 = t160 * t159;
	t162 = cos(qJ(4));
	t167 = t160 * t162;
	t163 = cos(qJ(1));
	t166 = t163 * t156;
	t165 = t163 * t159;
	t164 = t163 * t162;
	t161 = cos(qJ(5));
	t158 = sin(qJ(5));
	t155 = sin(t157);
	t1 = [0, -t160, -t160, t166, t155 * t165 - t167, -(t155 * t164 + t168) * t158 + t161 * t166; 0, t163, t163, t169, t155 * t168 + t164, -(t155 * t167 - t165) * t158 + t161 * t169; 1, 0, 0, -t155, t156 * t159, -t156 * t162 * t158 - t155 * t161;];
	Jg_rot = t1;
end