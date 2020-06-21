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
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
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
	Jg_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
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
	unknown(3,1) = -1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->20)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = t2;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = -1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->20)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = t2;
	unknown(2,3) = t2;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = -1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->3), mult. (2->2), div. (0->0), fcn. (9->4), ass. (0->23)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	t3 = qJ(2) + qJ(3);
	t4 = cos(t3);
	t7 = sin(t3);
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = (t2 * t4);
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = t2;
	unknown(2,3) = t2;
	unknown(2,4) = -(t1 * t4);
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = -1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = -t7;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (10->5), mult. (9->9), div. (0->0), fcn. (21->6), ass. (0->25)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	t3 = qJ(2) + qJ(3);
	t4 = cos(t3);
	t6 = sin(t3);
	t8 = sin(qJ(4));
	t10 = cos(qJ(4));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = (t2 * t4);
	unknown(1,5) = (t2 * t6 * t8 + t1 * t10);
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = t2;
	unknown(2,3) = t2;
	unknown(2,4) = -(t1 * t4);
	unknown(2,5) = (-t1 * t6 * t8 + t10 * t2);
	unknown(2,6) = 0;
	unknown(3,1) = -1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = -t6;
	unknown(3,5) = (t4 * t8);
	unknown(3,6) = 0;
	Jg_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (35->13), mult. (42->23), div. (0->0), fcn. (67->8), ass. (0->34)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	t3 = qJ(2) + qJ(3);
	t4 = cos(t3);
	t5 = t2 * t4;
	t6 = sin(t3);
	t7 = t2 * t6;
	t8 = sin(qJ(4));
	t10 = cos(qJ(4));
	t15 = sin(qJ(5));
	t17 = cos(qJ(5));
	t19 = -(-t1 * t8 + t7 * t10) * t15 + t5 * t17;
	t22 = t1 * t4;
	t23 = t1 * t6;
	t31 = -(-t23 * t10 - t2 * t8) * t15 - t17 * t22;
	t38 = -t4 * t10 * t15 - t6 * t17;
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = t5;
	unknown(1,5) = (-pkin(7) * t19 + t1 * t10 + t7 * t8);
	unknown(1,6) = t19;
	unknown(2,1) = 0;
	unknown(2,2) = t2;
	unknown(2,3) = t2;
	unknown(2,4) = -t22;
	unknown(2,5) = (-pkin(7) * t31 + t2 * t10 - t23 * t8);
	unknown(2,6) = t31;
	unknown(3,1) = -1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = -t6;
	unknown(3,5) = (-pkin(7) * t38 + t4 * t8);
	unknown(3,6) = t38;
	Jg_rot = unknown;
end