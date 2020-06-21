% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in CloosQRC350DE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = CloosQRC350DE_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
Ja_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
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
	Ja_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->2), mult. (6->3), div. (5->2), fcn. (6->2), ass. (0->24)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = t1 ^ 2;
	t3 = cos(qJ(1));
	t4 = t3 ^ 2;
	t6 = t2 / t4;
	t8 = 0.1e1 / (0.1e1 + t6);
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
	unknown(3,1) = (-t6 * t8 - t8);
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Ja_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (324->22), mult. (332->70), div. (79->9), fcn. (489->9), ass. (0->57)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = cos(t2);
	t4 = t1 * t3;
	t5 = sin(t2);
	t6 = 0.1e1 / t5;
	t7 = sin(qJ(1));
	t8 = t7 ^ 2;
	t9 = t3 ^ 2;
	t11 = t5 ^ 2;
	t12 = 0.1e1 / t11;
	t15 = 0.1e1 / (t8 * t9 * t12 + 0.1e1);
	t22 = t9 * t7 * t12 * t15 + t7 * t15;
	t23 = t7 * t3;
	t24 = atan2(t23, -t5);
	t25 = cos(t24);
	t27 = sin(t24);
	t28 = t27 * t7;
	t30 = -t25 * t5 + t28 * t3;
	t32 = t1 ^ 2;
	t34 = t30 ^ 2;
	t35 = 0.1e1 / t34;
	t38 = 0.1e1 / (t32 * t9 * t35 + 0.1e1);
	t39 = 0.1e1 / t30 * t38;
	t53 = t3 * t35 * t38;
	t56 = t1 * t5;
	t67 = -t56 * t39 - (t22 * t25 * t23 + t22 * t27 * t5 - t25 * t3 - t28 * t5) * t1 * t53;
	t68 = t7 * t5;
	t69 = sin(qJ(4));
	t71 = cos(qJ(4));
	t76 = t56 * t71 - t7 * t69;
	t77 = 0.1e1 / t76;
	t81 = t56 * t69 + t7 * t71;
	t82 = t81 ^ 2;
	t83 = t76 ^ 2;
	t84 = 0.1e1 / t83;
	t87 = 0.1e1 / (t82 * t84 + 0.1e1);
	t93 = t84 * t87;
	t103 = -t4 * t71 * t81 * t84 * t87 + t4 * t69 * t77 * t87;
	unknown(1,1) = -t4 * t6 * t15;
	unknown(1,2) = t22;
	unknown(1,3) = t22;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t23 * t39 - (-t1 * t9 * t6 * t15 * t25 * t7 + t27 * t1 * t3 - t4 * t15 * t27) * t1 * t53;
	unknown(2,2) = t67;
	unknown(2,3) = t67;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = (t1 * t71 - t68 * t69) * t77 * t87 - (-t1 * t69 - t68 * t71) * t81 * t93;
	unknown(3,2) = t103;
	unknown(3,3) = t103;
	unknown(3,4) = t81 ^ 2 * t93 + t87;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	Ja_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (600->32), mult. (873->108), div. (144->11), fcn. (1297->11), ass. (0->69)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = sin(t2);
	t4 = t1 * t3;
	t5 = sin(qJ(4));
	t7 = sin(qJ(1));
	t8 = cos(qJ(4));
	t10 = t4 * t5 + t7 * t8;
	t11 = cos(t2);
	t12 = 0.1e1 / t11;
	t13 = t10 * t12;
	t14 = 0.1e1 / t5;
	t15 = t7 * t3;
	t18 = -t1 * t8 + t5 * t15;
	t19 = t18 ^ 2;
	t20 = t11 ^ 2;
	t21 = 0.1e1 / t20;
	t23 = t5 ^ 2;
	t24 = 0.1e1 / t23;
	t27 = 0.1e1 / (t19 * t21 * t24 + 0.1e1);
	t28 = t14 * t27;
	t35 = t3 * t14 * t18 * t21 * t27 + t7 * t27;
	t38 = t1 * t5 + t15 * t8;
	t45 = -t12 * t8 * t18 * t24 * t27 + t38 * t12 * t28;
	t46 = t5 * t11;
	t47 = atan2(t18, t46);
	t48 = cos(t47);
	t49 = t48 * t11;
	t51 = sin(t47);
	t53 = t51 * t18 + t49 * t5;
	t54 = 0.1e1 / t53;
	t56 = t10 ^ 2;
	t57 = t53 ^ 2;
	t58 = 0.1e1 / t57;
	t61 = 0.1e1 / (t56 * t58 + 0.1e1);
	t72 = t58 * t61;
	t75 = t1 * t11;
	t90 = t75 * t5 * t54 * t61 - (t35 * t48 * t18 - t48 * t3 * t5 - t35 * t51 * t46 + t51 * t7 * t46) * t10 * t72;
	t93 = t4 * t8 - t7 * t5;
	t106 = sin(qJ(5));
	t108 = t7 * t11;
	t109 = cos(qJ(5));
	t114 = t75 * t106 + t93 * t109;
	t115 = 0.1e1 / t114;
	t119 = t93 * t106 - t75 * t109;
	t120 = t119 ^ 2;
	t121 = t114 ^ 2;
	t122 = 0.1e1 / t121;
	t125 = 0.1e1 / (t120 * t122 + 0.1e1);
	t131 = t122 * t125;
	t146 = (t75 * t8 * t106 + t4 * t109) * t115 * t125 - (t75 * t8 * t109 - t4 * t106) * t119 * t131;
	unknown(1,1) = t13 * t28;
	unknown(1,2) = t35;
	unknown(1,3) = t35;
	unknown(1,4) = t45;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t18 * t54 * t61 - (t13 * t14 * t27 * t48 * t18 - t10 * t27 * t51 + t51 * t10) * t10 * t72;
	unknown(2,2) = t90;
	unknown(2,3) = t90;
	unknown(2,4) = t93 * t54 * t61 - (t45 * t48 * t18 - t45 * t51 * t46 + t51 * t38 + t49 * t8) * t10 * t72;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = (-t38 * t106 + t108 * t109) * t115 * t125 - (-t108 * t106 - t38 * t109) * t119 * t131;
	unknown(3,2) = t146;
	unknown(3,3) = t146;
	unknown(3,4) = t10 * t109 * t119 * t122 * t125 - t10 * t106 * t115 * t125;
	unknown(3,5) = t119 ^ 2 * t131 + t125;
	unknown(3,6) = 0.0e0;
	Ja_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (1654->58), mult. (2306->159), div. (145->9), fcn. (3132->13), ass. (0->84)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = sin(t2);
	t4 = t1 * t3;
	t5 = cos(qJ(4));
	t7 = sin(qJ(1));
	t8 = sin(qJ(4));
	t10 = -t4 * t5 + t7 * t8;
	t11 = sin(qJ(5));
	t13 = cos(t2);
	t14 = t1 * t13;
	t15 = cos(qJ(5));
	t16 = t14 * t15;
	t17 = t10 * t11 + t16;
	t18 = t13 * t5;
	t21 = -t18 * t11 - t3 * t15;
	t22 = 0.1e1 / t21;
	t23 = t17 * t22;
	t24 = t7 * t3;
	t27 = -t1 * t8 - t24 * t5;
	t29 = t7 * t13;
	t31 = t27 * t11 + t29 * t15;
	t32 = t31 ^ 2;
	t33 = t21 ^ 2;
	t34 = 0.1e1 / t33;
	t37 = 0.1e1 / (t32 * t34 + 0.1e1);
	t39 = t5 * t11;
	t42 = -t24 * t15 - t29 * t39;
	t48 = t3 * t5 * t11 - t13 * t15;
	t50 = t34 * t37;
	t52 = t42 * t22 * t37 - t48 * t31 * t50;
	t55 = -t1 * t5 + t24 * t8;
	t64 = -t13 * t8 * t11 * t31 * t34 * t37 + t55 * t11 * t22 * t37;
	t67 = -t29 * t11 + t27 * t15;
	t72 = t3 * t11 - t18 * t15;
	t75 = t67 * t22 * t37 - t72 * t31 * t50;
	t76 = atan2(t31, t21);
	t77 = cos(t76);
	t79 = sin(t76);
	t81 = t77 * t21 + t79 * t31;
	t82 = 0.1e1 / t81;
	t85 = t10 * t11 + t16;
	t86 = t85 ^ 2;
	t87 = t81 ^ 2;
	t88 = 0.1e1 / t87;
	t91 = 0.1e1 / (t86 * t88 + 0.1e1);
	t101 = t88 * t91;
	t118 = (-t14 * t39 - t4 * t15) * t82 * t91 - (-t52 * t79 * t21 + t52 * t77 * t31 + t79 * t42 + t77 * t48) * t85 * t101;
	t121 = -t4 * t8 - t7 * t5;
	t140 = t10 * t15 - t14 * t11;
	t154 = pkin(7) * qJ(5) - qJ(6);
	t155 = sin(t154);
	t157 = cos(t154);
	t162 = t121 * t155 + t140 * t157;
	t163 = 0.1e1 / t162;
	t167 = t121 * t157 - t140 * t155;
	t168 = t167 ^ 2;
	t169 = t162 ^ 2;
	t170 = 0.1e1 / t169;
	t173 = 0.1e1 / (t168 * t170 + 0.1e1);
	t179 = t170 * t173;
	t185 = t14 * t5 * t15 - t4 * t11;
	t198 = (-t14 * t8 * t157 + t185 * t155) * t163 * t173 - (-t14 * t8 * t155 - t185 * t157) * t167 * t179;
	t199 = t121 * t15;
	t212 = -t140 * pkin(7);
	t214 = -t121 * pkin(7);
	unknown(1,1) = t23 * t37;
	unknown(1,2) = t52;
	unknown(1,3) = t52;
	unknown(1,4) = t64;
	unknown(1,5) = t75;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t31 * t82 * t91 - (t23 * t37 * t77 * t31 - t17 * t37 * t79 + t79 * t17) * t85 * t101;
	unknown(2,2) = t118;
	unknown(2,3) = t118;
	unknown(2,4) = -t121 * t11 * t82 * t91 - (t77 * t13 * t8 * t11 + t79 * t55 * t11 - t64 * t79 * t21 + t64 * t77 * t31) * t85 * t101;
	unknown(2,5) = t140 * t82 * t91 - (-t75 * t79 * t21 + t75 * t77 * t31 + t79 * t67 + t77 * t72) * t85 * t101;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = (t67 * t155 + t55 * t157) * t163 * t173 - (t55 * t155 - t67 * t157) * t167 * t179;
	unknown(3,2) = t198;
	unknown(3,3) = t198;
	unknown(3,4) = (t10 * t157 + t199 * t155) * t163 * t173 - (t10 * t155 - t199 * t157) * t167 * t179;
	unknown(3,5) = (t214 * t155 + t85 * t155 + t212 * t157) * t163 * t173 - (t212 * t155 - t214 * t157 - t85 * t157) * t167 * t179;
	unknown(3,6) = t167 ^ 2 * t179 + t173;
	Ja_rot = unknown;
end