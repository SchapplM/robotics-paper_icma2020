% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
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
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = CloosQRC350DE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
JR_rot=NaN(9,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
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
	JR_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->56)
	unknown=NaN(9,6);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = -t1;
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t2;
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
	unknown(4,1) = t2;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t1;
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
	JR_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->62)
	unknown=NaN(9,6);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = cos(qJ(1));
	t5 = cos(qJ(2));
	t6 = t4 * t5;
	t7 = t4 * t2;
	t8 = t1 * t5;
	unknown(1,1) = -t3;
	unknown(1,2) = t6;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t7;
	unknown(2,2) = -t8;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t2;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = -t8;
	unknown(4,2) = -t7;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t6;
	unknown(5,2) = t3;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t5;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = t4;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = -t1;
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
	JR_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->63)
	unknown=NaN(9,6);
	t1 = sin(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = sin(t2);
	t4 = t1 * t3;
	t5 = cos(qJ(1));
	t6 = cos(t2);
	t7 = t5 * t6;
	t8 = t5 * t3;
	t9 = t1 * t6;
	unknown(1,1) = -t4;
	unknown(1,2) = t7;
	unknown(1,3) = t7;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t8;
	unknown(2,2) = -t9;
	unknown(2,3) = -t9;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t3;
	unknown(3,3) = -t3;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = -t9;
	unknown(4,2) = -t8;
	unknown(4,3) = -t8;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t7;
	unknown(5,2) = t4;
	unknown(5,3) = t4;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t6;
	unknown(6,3) = -t6;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = t5;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = -t1;
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
	JR_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (52->21), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->75)
	unknown=NaN(9,6);
	t1 = sin(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = sin(t2);
	t4 = t1 * t3;
	t5 = cos(qJ(4));
	t7 = cos(qJ(1));
	t8 = sin(qJ(4));
	t10 = -t4 * t5 - t7 * t8;
	t11 = cos(t2);
	t12 = t7 * t11;
	t13 = t12 * t5;
	t14 = t7 * t3;
	t17 = -t1 * t5 - t14 * t8;
	t20 = t1 * t8 - t14 * t5;
	t21 = t1 * t11;
	t22 = t21 * t5;
	t25 = t4 * t8 - t7 * t5;
	t26 = t3 * t5;
	t28 = t12 * t8;
	t29 = t21 * t8;
	t30 = t3 * t8;
	unknown(1,1) = t10;
	unknown(1,2) = t13;
	unknown(1,3) = t13;
	unknown(1,4) = t17;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = t20;
	unknown(2,2) = -t22;
	unknown(2,3) = -t22;
	unknown(2,4) = t25;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t26;
	unknown(3,3) = -t26;
	unknown(3,4) = -t11 * t8;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = t25;
	unknown(4,2) = -t28;
	unknown(4,3) = -t28;
	unknown(4,4) = t20;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = -t17;
	unknown(5,2) = t29;
	unknown(5,3) = t29;
	unknown(5,4) = -t10;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t30;
	unknown(6,3) = t30;
	unknown(6,4) = -t11 * t5;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = -t21;
	unknown(7,2) = -t14;
	unknown(7,3) = -t14;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = -t12;
	unknown(8,2) = t4;
	unknown(8,3) = t4;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = -t11;
	unknown(9,3) = -t11;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (106->28), mult. (149->52), div. (0->0), fcn. (226->8), ass. (0->89)
	unknown=NaN(9,6);
	t1 = sin(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = sin(t2);
	t4 = t1 * t3;
	t5 = cos(qJ(4));
	t7 = cos(qJ(1));
	t8 = sin(qJ(4));
	t10 = -t4 * t5 - t7 * t8;
	t11 = cos(qJ(5));
	t13 = cos(t2);
	t14 = t1 * t13;
	t15 = sin(qJ(5));
	t17 = t10 * t11 - t14 * t15;
	t18 = t7 * t13;
	t19 = t5 * t11;
	t21 = t7 * t3;
	t23 = -t21 * t15 + t18 * t19;
	t26 = -t1 * t5 - t21 * t8;
	t30 = -t1 * t8 + t21 * t5;
	t32 = t18 * t11;
	t35 = t18 * t15;
	t39 = -t14 * t19 + t4 * t15;
	t42 = t4 * t8 - t7 * t5;
	t46 = -t10 * t15 - t14 * t11;
	t47 = t3 * t5;
	t50 = -t47 * t11 - t13 * t15;
	t51 = t13 * t8;
	t53 = t13 * t5;
	t57 = t5 * t15;
	t60 = -t21 * t11 - t18 * t57;
	t68 = t4 * t11 + t14 * t57;
	t72 = -t13 * t11 + t47 * t15;
	t77 = t18 * t8;
	t78 = t14 * t8;
	t79 = t3 * t8;
	unknown(1,1) = t17;
	unknown(1,2) = t23;
	unknown(1,3) = t23;
	unknown(1,4) = t26 * t11;
	unknown(1,5) = -t30 * t15 + t32;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t30 * t11 - t35;
	unknown(2,2) = t39;
	unknown(2,3) = t39;
	unknown(2,4) = t42 * t11;
	unknown(2,5) = t46;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t50;
	unknown(3,3) = t50;
	unknown(3,4) = -t51 * t11;
	unknown(3,5) = -t3 * t11 - t53 * t15;
	unknown(3,6) = 0.0e0;
	unknown(4,1) = t46;
	unknown(4,2) = t60;
	unknown(4,3) = t60;
	unknown(4,4) = -t26 * t15;
	unknown(4,5) = -t30 * t11 - t35;
	unknown(4,6) = 0.0e0;
	unknown(5,1) = t30 * t15 - t32;
	unknown(5,2) = t68;
	unknown(5,3) = t68;
	unknown(5,4) = -t42 * t15;
	unknown(5,5) = -t17;
	unknown(5,6) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t72;
	unknown(6,3) = t72;
	unknown(6,4) = t51 * t15;
	unknown(6,5) = -t53 * t11 + t3 * t15;
	unknown(6,6) = 0.0e0;
	unknown(7,1) = -t42;
	unknown(7,2) = t77;
	unknown(7,3) = t77;
	unknown(7,4) = t30;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = t26;
	unknown(8,2) = -t78;
	unknown(8,3) = -t78;
	unknown(8,4) = t10;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = -t79;
	unknown(9,3) = -t79;
	unknown(9,4) = t53;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:15
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (336->59), mult. (483->119), div. (0->0), fcn. (572->10), ass. (0->112)
	unknown=NaN(9,6);
	t1 = sin(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = sin(t2);
	t4 = t1 * t3;
	t5 = cos(qJ(4));
	t7 = cos(qJ(1));
	t8 = sin(qJ(4));
	t10 = -t4 * t5 - t7 * t8;
	t11 = cos(qJ(5));
	t13 = cos(t2);
	t14 = t1 * t13;
	t15 = sin(qJ(5));
	t17 = t10 * t11 - t14 * t15;
	t19 = pkin(7) * qJ(5) - qJ(6);
	t20 = cos(t19);
	t24 = -t4 * t8 + t7 * t5;
	t25 = sin(t19);
	t27 = -t17 * t20 - t24 * t25;
	t28 = t7 * t13;
	t29 = t5 * t11;
	t31 = t7 * t3;
	t33 = -t31 * t15 + t28 * t29;
	t35 = t8 * t25;
	t37 = -t33 * t20 - t28 * t35;
	t40 = -t1 * t5 - t31 * t8;
	t41 = t40 * t11;
	t45 = -t1 * t8 + t31 * t5;
	t49 = t28 * t11;
	t50 = -t45 * t15 + t49;
	t53 = t28 * t15;
	t54 = t45 * t11 + t53;
	t55 = t54 * pkin(7);
	t57 = -t40 * pkin(7);
	t64 = -t45 * t11 - t53;
	t70 = -t14 * t29 + t4 * t15;
	t73 = t14 * t35 - t70 * t20;
	t74 = -t24 * t11;
	t80 = -t10 * t15 - t14 * t11;
	t82 = t17 * pkin(7);
	t84 = t24 * pkin(7);
	t89 = -t17 * t25 + t24 * t20;
	t90 = t3 * t5;
	t93 = -t90 * t11 - t13 * t15;
	t95 = t3 * t8;
	t97 = -t93 * t20 + t95 * t25;
	t98 = t13 * t8;
	t101 = t13 * t5;
	t106 = -t101 * t15 - t3 * t11;
	t110 = t101 * t11 - t3 * t15;
	t111 = t110 * pkin(7);
	t120 = t8 * t20;
	t122 = t28 * t120 - t33 * t25;
	t138 = -t14 * t120 - t70 * t25;
	t148 = -t95 * t20 - t93 * t25;
	t161 = t5 * t15;
	t164 = -t31 * t11 - t28 * t161;
	t170 = t4 * t11 + t14 * t161;
	t174 = -t13 * t11 + t90 * t15;
	unknown(1,1) = t27;
	unknown(1,2) = t37;
	unknown(1,3) = t37;
	unknown(1,4) = -t41 * t20 - t45 * t25;
	unknown(1,5) = -t50 * t20 - t57 * t20 + t55 * t25;
	unknown(1,6) = -t40 * t20 - t54 * t25;
	unknown(2,1) = -t64 * t20 - t40 * t25;
	unknown(2,2) = t73;
	unknown(2,3) = t73;
	unknown(2,4) = -t10 * t25 - t74 * t20;
	unknown(2,5) = -t80 * t20 - t84 * t20 + t82 * t25;
	unknown(2,6) = t89;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t97;
	unknown(3,3) = t97;
	unknown(3,4) = t98 * t11 * t20 - t101 * t25;
	unknown(3,5) = -t98 * pkin(7) * t20 - t106 * t20 + t111 * t25;
	unknown(3,6) = -t110 * t25 + t98 * t20;
	unknown(4,1) = t89;
	unknown(4,2) = t122;
	unknown(4,3) = t122;
	unknown(4,4) = t45 * t20 - t41 * t25;
	unknown(4,5) = -t55 * t20 - t50 * t25 - t57 * t25;
	unknown(4,6) = t54 * t20 - t40 * t25;
	unknown(5,1) = t40 * t20 - t64 * t25;
	unknown(5,2) = t138;
	unknown(5,3) = t138;
	unknown(5,4) = t10 * t20 - t74 * t25;
	unknown(5,5) = -t82 * t20 - t80 * t25 - t84 * t25;
	unknown(5,6) = -t27;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t148;
	unknown(6,3) = t148;
	unknown(6,4) = t98 * t11 * t25 + t101 * t20;
	unknown(6,5) = -t98 * pkin(7) * t25 - t106 * t25 - t111 * t20;
	unknown(6,6) = t110 * t20 + t98 * t25;
	unknown(7,1) = t80;
	unknown(7,2) = t164;
	unknown(7,3) = t164;
	unknown(7,4) = -t40 * t15;
	unknown(7,5) = -t54;
	unknown(7,6) = 0.0e0;
	unknown(8,1) = t45 * t15 - t49;
	unknown(8,2) = t170;
	unknown(8,3) = t170;
	unknown(8,4) = t24 * t15;
	unknown(8,5) = -t17;
	unknown(8,6) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t174;
	unknown(9,3) = t174;
	unknown(9,4) = t98 * t15;
	unknown(9,5) = -t110;
	unknown(9,6) = 0.0e0;
	JR_rot = unknown;
end