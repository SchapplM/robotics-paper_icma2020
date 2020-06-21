% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
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
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = CloosQRC350DE_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
JRD_rot=NaN(9,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->54)
	unknown=NaN(9,6);
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
	unknown(4,1) = 0;
	unknown(4,2) = 0;
	unknown(4,3) = 0;
	unknown(4,4) = 0;
	unknown(4,5) = 0;
	unknown(4,6) = 0;
	unknown(5,1) = 0;
	unknown(5,2) = 0;
	unknown(5,3) = 0;
	unknown(5,4) = 0;
	unknown(5,5) = 0;
	unknown(5,6) = 0;
	unknown(6,1) = 0;
	unknown(6,2) = 0;
	unknown(6,3) = 0;
	unknown(6,4) = 0;
	unknown(6,5) = 0;
	unknown(6,6) = 0;
	unknown(7,1) = 0;
	unknown(7,2) = 0;
	unknown(7,3) = 0;
	unknown(7,4) = 0;
	unknown(7,5) = 0;
	unknown(7,6) = 0;
	unknown(8,1) = 0;
	unknown(8,2) = 0;
	unknown(8,3) = 0;
	unknown(8,4) = 0;
	unknown(8,5) = 0;
	unknown(8,6) = 0;
	unknown(9,1) = 0;
	unknown(9,2) = 0;
	unknown(9,3) = 0;
	unknown(9,4) = 0;
	unknown(9,5) = 0;
	unknown(9,6) = 0;
	JRD_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->58)
	unknown=NaN(9,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	unknown(1,1) = -t2;
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = t4;
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = -t4;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t2;
	unknown(5,2) = 0.0e0;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = 0.0e0;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = 0.0e0;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = 0.0e0;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	JRD_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->9), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->66)
	unknown=NaN(9,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(2));
	t5 = sin(qJ(1));
	t6 = t5 * qJD(2);
	t7 = cos(qJ(2));
	t9 = -t2 * t3 - t6 * t7;
	t10 = qJD(1) * t5;
	t12 = t1 * qJD(2);
	t14 = -t10 * t7 - t12 * t3;
	t17 = t10 * t3 - t12 * t7;
	t20 = -t2 * t7 + t6 * t3;
	unknown(1,1) = t9;
	unknown(1,2) = t14;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = t17;
	unknown(2,2) = t20;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -qJD(2) * t7;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = t20;
	unknown(4,2) = t17;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t14;
	unknown(5,2) = -t9;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = qJD(2) * t3;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = -t10;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = -t2;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	JRD_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (60->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->70)
	unknown=NaN(9,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = qJ(2) + qJ(3);
	t4 = sin(t3);
	t6 = sin(qJ(1));
	t7 = qJD(2) + qJD(3);
	t8 = t6 * t7;
	t9 = cos(t3);
	t11 = -t2 * t4 - t8 * t9;
	t12 = qJD(1) * t6;
	t14 = t1 * t7;
	t16 = -t12 * t9 - t14 * t4;
	t19 = t12 * t4 - t14 * t9;
	t22 = -t2 * t9 + t8 * t4;
	t23 = t7 * t9;
	t24 = t7 * t4;
	unknown(1,1) = t11;
	unknown(1,2) = t16;
	unknown(1,3) = t16;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = t19;
	unknown(2,2) = t22;
	unknown(2,3) = t22;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t23;
	unknown(3,3) = -t23;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = t22;
	unknown(4,2) = t19;
	unknown(4,3) = t19;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t16;
	unknown(5,2) = -t11;
	unknown(5,3) = -t11;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t24;
	unknown(6,3) = t24;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = -t12;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = -t2;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	JRD_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:15
	% EndTime: 2020-06-19 21:40:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (164->36), mult. (226->68), div. (0->0), fcn. (226->6), ass. (0->94)
	unknown=NaN(9,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = qJ(2) + qJ(3);
	t4 = sin(t3);
	t5 = cos(qJ(4));
	t6 = t4 * t5;
	t8 = sin(qJ(1));
	t9 = qJD(2) + qJD(3);
	t10 = t8 * t9;
	t11 = cos(t3);
	t12 = t11 * t5;
	t14 = t8 * t4;
	t15 = sin(qJ(4));
	t16 = qJD(4) * t15;
	t18 = qJD(1) * t8;
	t20 = t1 * qJD(4);
	t22 = -t10 * t12 + t14 * t16 + t15 * t18 - t2 * t6 - t20 * t5;
	t24 = t1 * t9;
	t26 = t1 * t11;
	t28 = -t18 * t12 - t26 * t16 - t24 * t6;
	t29 = t4 * t15;
	t31 = t11 * t15;
	t33 = t1 * t4;
	t34 = qJD(4) * t5;
	t37 = t8 * qJD(4);
	t39 = t37 * t15 + t18 * t29 - t2 * t5 - t24 * t31 - t33 * t34;
	t45 = -t24 * t12 + t15 * t2 + t33 * t16 + t18 * t6 + t37 * t5;
	t48 = t8 * t11;
	t50 = t10 * t6 - t2 * t12 + t48 * t16;
	t56 = t10 * t31 + t14 * t34 + t15 * t20 + t18 * t5 + t2 * t29;
	t57 = t9 * t11;
	t59 = t4 * qJD(4);
	t61 = t59 * t15 - t57 * t5;
	t62 = t9 * t4;
	t64 = t11 * qJD(4);
	t70 = t18 * t31 + t24 * t29 - t26 * t34;
	t74 = -t10 * t29 + t2 * t31 + t48 * t34;
	t77 = t57 * t15 + t59 * t5;
	t86 = -t24 * t11 + t18 * t4;
	t92 = t10 * t11 + t2 * t4;
	unknown(1,1) = t22;
	unknown(1,2) = t28;
	unknown(1,3) = t28;
	unknown(1,4) = t39;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = t45;
	unknown(2,2) = t50;
	unknown(2,3) = t50;
	unknown(2,4) = t56;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t61;
	unknown(3,3) = t61;
	unknown(3,4) = t62 * t15 - t64 * t5;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = t56;
	unknown(4,2) = t70;
	unknown(4,3) = t70;
	unknown(4,4) = t45;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t39;
	unknown(5,2) = t74;
	unknown(5,3) = t74;
	unknown(5,4) = -t22;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t77;
	unknown(6,3) = t77;
	unknown(6,4) = t64 * t15 + t62 * t5;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = t10 * t4 - t2 * t11;
	unknown(7,2) = t86;
	unknown(7,3) = t86;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = t18 * t11 + t24 * t4;
	unknown(8,2) = t92;
	unknown(8,3) = t92;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t62;
	unknown(9,3) = t62;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	JRD_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:15
	% EndTime: 2020-06-19 21:40:16
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (454->103), mult. (696->170), div. (0->0), fcn. (708->8), ass. (0->132)
	unknown=NaN(9,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = qJ(2) + qJ(3);
	t4 = sin(t3);
	t5 = cos(qJ(4));
	t6 = t4 * t5;
	t8 = sin(qJ(1));
	t9 = qJD(2) + qJD(3);
	t10 = t8 * t9;
	t11 = cos(t3);
	t12 = t11 * t5;
	t14 = t8 * t4;
	t15 = sin(qJ(4));
	t16 = qJD(4) * t15;
	t18 = qJD(1) * t8;
	t20 = t1 * qJD(4);
	t22 = -t10 * t12 + t14 * t16 + t18 * t15 - t2 * t6 - t20 * t5;
	t23 = cos(qJ(5));
	t28 = (-t1 * t15 - t14 * t5) * qJD(5);
	t29 = sin(qJ(5));
	t31 = t11 * t29;
	t33 = t4 * t29;
	t35 = t8 * t11;
	t36 = qJD(5) * t23;
	t38 = t10 * t33 - t2 * t31 + t22 * t23 - t28 * t29 - t35 * t36;
	t39 = t12 * t23;
	t41 = t1 * t9;
	t42 = t6 * t23;
	t44 = t1 * t11;
	t45 = t16 * t23;
	t47 = t5 * qJD(5);
	t48 = t47 * t29;
	t52 = t1 * t4;
	t54 = t18 * t33 - t18 * t39 - t41 * t31 - t52 * t36 - t41 * t42 - t44 * t45 - t44 * t48;
	t55 = t4 * t15;
	t57 = t11 * t15;
	t59 = qJD(4) * t5;
	t62 = t8 * qJD(4);
	t64 = t62 * t15 + t18 * t55 - t2 * t5 - t41 * t57 - t52 * t59;
	t69 = (-t52 * t15 - t8 * t5) * qJD(5);
	t77 = t41 * t12 - t2 * t15 - t52 * t16 - t18 * t6 - t62 * t5;
	t81 = -t8 * t15 + t52 * t5;
	t82 = t81 * qJD(5);
	t84 = t11 * t23;
	t85 = t18 * t84;
	t86 = t4 * t23;
	t87 = t41 * t86;
	t88 = qJD(5) * t29;
	t89 = t44 * t88;
	t92 = -t81 * qJD(5);
	t94 = t18 * t31;
	t95 = t41 * t33;
	t96 = t44 * t36;
	t105 = t10 * t31 + t10 * t42 + t14 * t36 + t2 * t33 - t2 * t39 + t35 * t45 + t35 * t48;
	t111 = t10 * t57 + t14 * t59 + t20 * t15 + t18 * t5 + t2 * t55;
	t116 = (-t1 * t5 + t14 * t15) * qJD(5);
	t124 = t10 * t86 - t2 * t84 - t22 * t29 - t28 * t23 + t35 * t88;
	t125 = t9 * t11;
	t126 = t5 * t23;
	t128 = t4 * qJD(4);
	t129 = t23 * t15;
	t132 = t9 * t4;
	t134 = t11 * qJD(5);
	t136 = -t125 * t126 + t128 * t129 + t132 * t29 - t134 * t23 + t6 * t88;
	t138 = t11 * qJD(4);
	t142 = t5 * t29;
	t144 = t15 * t29;
	t148 = t4 * qJD(5);
	t151 = t12 * t29;
	t153 = t6 * t29;
	t155 = t16 * t29;
	t157 = t47 * t23;
	t162 = t18 * t151 + t41 * t153 + t44 * t155 - t44 * t157 + t18 * t86 - t41 * t84 + t52 * t88;
	t179 = -t10 * t153 + t10 * t84 - t14 * t88 + t2 * t151 - t35 * t155 + t35 * t157 + t2 * t86;
	t188 = t125 * t142 - t128 * t144 + t132 * t23 + t134 * t29 + t6 * t36;
	t202 = -t18 * t57 - t41 * t55 + t44 * t59;
	t206 = t10 * t55 - t2 * t57 - t35 * t59;
	t209 = -t125 * t15 - t128 * t5;
	unknown(1,1) = t38;
	unknown(1,2) = t54;
	unknown(1,3) = t54;
	unknown(1,4) = t64 * t23 - t69 * t29;
	unknown(1,5) = -t82 * t23 - t77 * t29 - t85 - t87 - t89;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t77 * t23 - t92 * t29 + t94 + t95 - t96;
	unknown(2,2) = t105;
	unknown(2,3) = t105;
	unknown(2,4) = t111 * t23 - t116 * t29;
	unknown(2,5) = t124;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t136;
	unknown(3,3) = t136;
	unknown(3,4) = -t138 * t126 + t132 * t129 + t57 * t88;
	unknown(3,5) = -t12 * t36 - t125 * t23 + t132 * t142 + t138 * t144 + t148 * t29;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = t124;
	unknown(4,2) = t162;
	unknown(4,3) = t162;
	unknown(4,4) = -t69 * t23 - t64 * t29;
	unknown(4,5) = -t77 * t23 + t82 * t29 + t94 + t95 - t96;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t92 * t23 + t77 * t29 + t85 + t87 + t89;
	unknown(5,2) = t179;
	unknown(5,3) = t179;
	unknown(5,4) = -t111 * t29 - t116 * t23;
	unknown(5,5) = -t38;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t188;
	unknown(6,3) = t188;
	unknown(6,4) = -t132 * t144 + t138 * t142 + t57 * t36;
	unknown(6,5) = t12 * t88 + t125 * t29 + t132 * t126 + t138 * t129 + t148 * t23;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = -t111;
	unknown(7,2) = t202;
	unknown(7,3) = t202;
	unknown(7,4) = t77;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = t64;
	unknown(8,2) = t206;
	unknown(8,3) = t206;
	unknown(8,4) = t22;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t209;
	unknown(9,3) = t209;
	unknown(9,4) = -t132 * t5 - t138 * t15;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	JRD_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:16
	% EndTime: 2020-06-19 21:40:17
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (1531->227), mult. (2288->372), div. (0->0), fcn. (2102->10), ass. (0->203)
	unknown=NaN(9,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = qJ(2) + qJ(3);
	t4 = sin(t3);
	t5 = cos(qJ(4));
	t6 = t4 * t5;
	t8 = sin(qJ(1));
	t9 = qJD(2) + qJD(3);
	t10 = t8 * t9;
	t11 = cos(t3);
	t12 = t11 * t5;
	t14 = t8 * t4;
	t15 = sin(qJ(4));
	t16 = qJD(4) * t15;
	t18 = qJD(1) * t8;
	t20 = t1 * qJD(4);
	t22 = -t10 * t12 + t14 * t16 + t15 * t18 - t2 * t6 - t20 * t5;
	t23 = cos(qJ(5));
	t27 = -t1 * t15 - t14 * t5;
	t28 = t27 * qJD(5);
	t29 = sin(qJ(5));
	t31 = t11 * t29;
	t33 = t4 * t29;
	t35 = t8 * t11;
	t36 = qJD(5) * t23;
	t38 = t10 * t33 - t2 * t31 + t22 * t23 - t28 * t29 - t35 * t36;
	t40 = pkin(7) * qJ(5) - qJ(6);
	t41 = cos(t40);
	t45 = t23 * t27 - t29 * t35;
	t47 = pkin(7) * qJD(5) - qJD(6);
	t48 = t45 * t47;
	t49 = sin(t40);
	t51 = t4 * t15;
	t53 = t11 * t15;
	t55 = qJD(4) * t5;
	t59 = -t10 * t53 - t14 * t55 - t15 * t20 - t18 * t5 - t2 * t51;
	t63 = t1 * t5 - t14 * t15;
	t64 = t63 * t47;
	t66 = -t38 * t41 - t41 * t64 + t48 * t49 - t49 * t59;
	t67 = t12 * t23;
	t69 = t1 * t9;
	t70 = t6 * t23;
	t72 = t1 * t11;
	t73 = t16 * t23;
	t75 = t5 * qJD(5);
	t76 = t75 * t29;
	t80 = t1 * t4;
	t82 = t18 * t33 - t18 * t67 - t31 * t69 - t36 * t80 - t69 * t70 - t72 * t73 - t72 * t76;
	t84 = t5 * t23;
	t88 = (-t29 * t80 + t72 * t84) * t47;
	t90 = t53 * t49;
	t92 = t51 * t49;
	t94 = t55 * t49;
	t96 = t15 * t47;
	t97 = t96 * t41;
	t99 = t18 * t90 - t41 * t82 + t49 * t88 + t69 * t92 - t72 * t94 - t72 * t97;
	t104 = t8 * qJD(4);
	t106 = t104 * t15 + t18 * t51 - t2 * t5 - t53 * t69 - t55 * t80;
	t107 = t106 * t23;
	t111 = -t15 * t80 - t5 * t8;
	t112 = t111 * qJD(5);
	t113 = t29 * t41;
	t115 = t111 * t23;
	t116 = t47 * t49;
	t123 = -t104 * t5 + t12 * t69 - t15 * t2 - t16 * t80 - t18 * t6;
	t127 = -t15 * t8 + t5 * t80;
	t128 = t127 * t47;
	t132 = t127 * qJD(5);
	t134 = t11 * t23;
	t135 = t18 * t134;
	t136 = t4 * t23;
	t137 = t69 * t136;
	t138 = qJD(5) * t29;
	t139 = t72 * t138;
	t140 = -t123 * t29 - t132 * t23 - t135 - t137 - t139;
	t145 = (-t127 * t29 + t23 * t72) * t47;
	t149 = t18 * t31;
	t150 = t69 * t33;
	t151 = t72 * t36;
	t152 = t123 * t23 - t132 * t29 - t149 - t150 + t151;
	t153 = t152 * pkin(7);
	t156 = t72 * t29;
	t157 = t127 * t23 + t156;
	t158 = t157 * pkin(7);
	t159 = t47 * t41;
	t161 = -t106 * pkin(7);
	t163 = -t111 * pkin(7);
	t167 = t157 * t47;
	t170 = -t111 * t47;
	t174 = -t127 * qJD(5);
	t176 = -t123 * t23 - t174 * t29 + t149 + t150 - t151;
	t180 = (-t127 * t23 - t156) * t47;
	t183 = t111 * t47;
	t193 = t10 * t31 + t10 * t70 + t14 * t36 + t2 * t33 - t2 * t67 + t35 * t73 + t35 * t76;
	t198 = (t14 * t29 - t35 * t84) * t47;
	t204 = -t10 * t92 - t193 * t41 + t198 * t49 + t2 * t90 + t35 * t94 + t35 * t97;
	t205 = -t59 * t23;
	t207 = -t63 * qJD(5);
	t209 = -t63 * t23;
	t212 = t27 * t47;
	t220 = t10 * t136 - t134 * t2 + t138 * t35 - t22 * t29 - t23 * t28;
	t225 = (-t23 * t35 - t27 * t29) * t47;
	t227 = t38 * pkin(7);
	t229 = t45 * pkin(7);
	t231 = t59 * pkin(7);
	t233 = t63 * pkin(7);
	t240 = -t38 * t49 - t41 * t48 + t41 * t59 - t49 * t64;
	t241 = t9 * t11;
	t243 = t4 * qJD(4);
	t244 = t15 * t23;
	t247 = t9 * t4;
	t249 = t11 * qJD(5);
	t251 = t138 * t6 - t23 * t249 - t241 * t84 + t243 * t244 + t247 * t29;
	t254 = (-t70 - t31) * t47;
	t256 = t15 * t49;
	t258 = t5 * t49;
	t261 = t159 * t51 + t241 * t256 + t243 * t258 - t251 * t41 + t254 * t49;
	t264 = t11 * qJD(4);
	t269 = t23 * t47;
	t276 = t5 * t29;
	t278 = t15 * t29;
	t282 = t4 * qJD(5);
	t284 = -t12 * t36 - t23 * t241 + t247 * t276 + t264 * t278 + t282 * t29;
	t286 = t12 * t29;
	t288 = (-t286 - t136) * t47;
	t295 = -t12 * t138 - t23 * t282 - t241 * t29 - t244 * t264 - t247 * t84;
	t296 = t295 * pkin(7);
	t298 = t67 - t33;
	t299 = t298 * pkin(7);
	t301 = t15 * pkin(7);
	t304 = t5 * pkin(7);
	t307 = pkin(7) * t47;
	t312 = t298 * t47;
	t314 = t15 * t41;
	t316 = t5 * t41;
	t322 = t53 * t41;
	t324 = t51 * t41;
	t326 = t55 * t41;
	t328 = t96 * t49;
	t330 = -t18 * t322 - t324 * t69 + t326 * t72 - t328 * t72 - t41 * t88 - t49 * t82;
	t332 = t29 * t49;
	t361 = t10 * t324 - t193 * t49 - t198 * t41 - t2 * t322 - t326 * t35 + t328 * t35;
	t380 = t116 * t51 - t241 * t314 - t243 * t316 - t251 * t49 - t254 * t41;
	t411 = t6 * t29;
	t413 = t16 * t29;
	t415 = t75 * t23;
	t420 = -t134 * t69 + t136 * t18 + t138 * t80 + t18 * t286 + t411 * t69 + t413 * t72 - t415 * t72;
	t434 = t10 * t134 - t10 * t411 + t136 * t2 - t138 * t14 + t2 * t286 - t35 * t413 + t35 * t415;
	t443 = t23 * t247 + t241 * t276 - t243 * t278 + t249 * t29 + t36 * t6;
	unknown(1,1) = t66;
	unknown(1,2) = t99;
	unknown(1,3) = t99;
	unknown(1,4) = -t107 * t41 + t112 * t113 + t115 * t116 - t123 * t49 - t128 * t41;
	unknown(1,5) = t116 * t163 - t140 * t41 + t145 * t49 + t153 * t49 + t158 * t159 - t161 * t41;
	unknown(1,6) = -t106 * t41 - t152 * t49 - t167 * t41 - t170 * t49;
	unknown(2,1) = -t106 * t49 - t176 * t41 + t180 * t49 - t183 * t41;
	unknown(2,2) = t204;
	unknown(2,3) = t204;
	unknown(2,4) = t113 * t207 + t116 * t209 - t205 * t41 - t212 * t41 - t22 * t49;
	unknown(2,5) = t116 * t233 + t159 * t229 - t220 * t41 + t225 * t49 + t227 * t49 - t231 * t41;
	unknown(2,6) = t240;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t261;
	unknown(3,3) = t261;
	unknown(3,4) = -t138 * t41 * t53 - t244 * t247 * t41 + t264 * t41 * t84 - t269 * t49 * t53 - t12 * t159 + t247 * t258 + t256 * t264;
	unknown(3,5) = t247 * t301 * t41 - t264 * t304 * t41 + t307 * t49 * t53 + t159 * t299 - t284 * t41 + t288 * t49 + t296 * t49;
	unknown(3,6) = -t116 * t53 - t247 * t314 + t264 * t316 - t295 * t49 - t312 * t41;
	unknown(4,1) = t240;
	unknown(4,2) = t330;
	unknown(4,3) = t330;
	unknown(4,4) = -t107 * t49 + t112 * t332 - t115 * t159 + t123 * t41 - t128 * t49;
	unknown(4,5) = t116 * t158 - t140 * t49 - t145 * t41 - t153 * t41 - t159 * t163 - t161 * t49;
	unknown(4,6) = -t106 * t49 + t152 * t41 - t167 * t49 + t170 * t41;
	unknown(5,1) = t106 * t41 - t176 * t49 - t180 * t41 - t183 * t49;
	unknown(5,2) = t361;
	unknown(5,3) = t361;
	unknown(5,4) = -t159 * t209 - t205 * t49 + t207 * t332 - t212 * t49 + t22 * t41;
	unknown(5,5) = t116 * t229 - t159 * t233 - t220 * t49 - t225 * t41 - t227 * t41 - t231 * t49;
	unknown(5,6) = -t66;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t380;
	unknown(6,3) = t380;
	unknown(6,4) = -t138 * t49 * t53 - t244 * t247 * t49 + t264 * t49 * t84 + t269 * t41 * t53 - t116 * t12 - t247 * t316 - t264 * t314;
	unknown(6,5) = t247 * t301 * t49 - t264 * t304 * t49 - t307 * t41 * t53 + t116 * t299 - t284 * t49 - t288 * t41 - t296 * t41;
	unknown(6,6) = t159 * t53 - t247 * t256 + t258 * t264 + t295 * t41 - t312 * t49;
	unknown(7,1) = t220;
	unknown(7,2) = t420;
	unknown(7,3) = t420;
	unknown(7,4) = -t106 * t29 - t112 * t23;
	unknown(7,5) = -t152;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = t123 * t29 - t174 * t23 + t135 + t137 + t139;
	unknown(8,2) = t434;
	unknown(8,3) = t434;
	unknown(8,4) = -t207 * t23 + t29 * t59;
	unknown(8,5) = -t38;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t443;
	unknown(9,3) = t443;
	unknown(9,4) = -t247 * t278 + t264 * t276 + t36 * t53;
	unknown(9,5) = -t295;
	unknown(9,6) = 0.0e0;
	JRD_rot = unknown;
end