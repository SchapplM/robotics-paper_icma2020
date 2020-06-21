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
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
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
	JgD_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
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
	JgD_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->20)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t3 = sin(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (qJD(1) * t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -(qJD(1) * t3);
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
	JgD_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->22)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	unknown(1,1) = 0;
	unknown(1,2) = t2;
	unknown(1,3) = t2;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t4;
	unknown(2,3) = -t4;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	JgD_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (13->7), mult. (13->9), div. (0->0), fcn. (13->4), ass. (0->26)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	t5 = qJ(2) + qJ(3);
	t6 = cos(t5);
	t8 = qJD(2) + qJD(3);
	t10 = sin(t5);
	unknown(1,1) = 0;
	unknown(1,2) = t2;
	unknown(1,3) = t2;
	unknown(1,4) = (-t1 * t8 * t10 - t4 * t6);
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t4;
	unknown(2,3) = -t4;
	unknown(2,4) = (t3 * t8 * t10 - t2 * t6);
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = -(t8 * t6);
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	JgD_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->16), mult. (43->30), div. (0->0), fcn. (43->6), ass. (0->33)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	t5 = qJ(2) + qJ(3);
	t6 = cos(t5);
	t8 = qJD(2) + qJD(3);
	t9 = t1 * t8;
	t10 = sin(t5);
	t13 = sin(qJ(4));
	t14 = t10 * t13;
	t16 = t6 * t13;
	t19 = cos(qJ(4));
	t20 = qJD(4) * t19;
	t27 = t3 * t8;
	unknown(1,1) = 0;
	unknown(1,2) = t2;
	unknown(1,3) = t2;
	unknown(1,4) = (-t9 * t10 - t4 * t6);
	unknown(1,5) = (-t3 * qJD(4) * t13 + t1 * t10 * t20 - t4 * t14 + t9 * t16 + t2 * t19);
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t4;
	unknown(2,3) = -t4;
	unknown(2,4) = (t27 * t10 - t2 * t6);
	unknown(2,5) = (-t1 * qJD(4) * t13 - t3 * t10 * t20 - t2 * t14 - t27 * t16 - t4 * t19);
	unknown(2,6) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = -(t8 * t6);
	unknown(3,5) = (t6 * qJD(4) * t19 - t8 * t10 * t13);
	unknown(3,6) = 0;
	JgD_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:15
	% EndTime: 2020-06-19 21:40:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (130->41), mult. (184->76), div. (0->0), fcn. (185->8), ass. (0->51)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	t5 = qJ(2) + qJ(3);
	t6 = cos(t5);
	t8 = qJD(2) + qJD(3);
	t9 = t1 * t8;
	t10 = sin(t5);
	t13 = sin(qJ(4));
	t14 = t10 * t13;
	t16 = t6 * t13;
	t18 = t1 * t10;
	t19 = cos(qJ(4));
	t20 = qJD(4) * t19;
	t23 = t3 * qJD(4);
	t25 = t10 * t19;
	t27 = t6 * t19;
	t29 = qJD(4) * t13;
	t34 = sin(qJ(5));
	t40 = cos(qJ(5));
	t42 = t6 * t40;
	t44 = t10 * t40;
	t47 = qJD(5) * t34;
	t49 = -(-t2 * t13 - t18 * t29 - t23 * t19 - t4 * t25 + t9 * t27) * t34 - (-t3 * t13 + t18 * t19) * qJD(5) * t40 - t4 * t42 - t9 * t44 - t1 * t6 * t47;
	t53 = t3 * t8;
	t58 = t3 * t10;
	t61 = t1 * qJD(4);
	t79 = -(t4 * t13 - t61 * t19 - t2 * t25 - t53 * t27 + t58 * t29) * t34 - (-t1 * t13 - t58 * t19) * qJD(5) * t40 - t2 * t42 + t53 * t44 + t3 * t6 * t47;
	t82 = t8 * t6;
	t83 = t8 * t10;
	t85 = t6 * qJD(4);
	t96 = t10 * qJD(5) * t34 - t27 * qJD(5) * t40 + t85 * t13 * t34 + t83 * t19 * t34 - t82 * t40;
	unknown(1,1) = 0;
	unknown(1,2) = t2;
	unknown(1,3) = t2;
	unknown(1,4) = (-t9 * t10 - t4 * t6);
	unknown(1,5) = (-pkin(7) * t49 - t23 * t13 - t4 * t14 + t9 * t16 + t18 * t20 + t2 * t19);
	unknown(1,6) = t49;
	unknown(2,1) = 0;
	unknown(2,2) = -t4;
	unknown(2,3) = -t4;
	unknown(2,4) = (t53 * t10 - t2 * t6);
	unknown(2,5) = (-pkin(7) * t79 - t61 * t13 - t2 * t14 - t53 * t16 - t4 * t19 - t58 * t20);
	unknown(2,6) = t79;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = -t82;
	unknown(3,5) = (-pkin(7) * t96 - t83 * t13 + t85 * t19);
	unknown(3,6) = t96;
	JgD_rot = unknown;
end