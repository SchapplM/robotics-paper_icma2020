% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% CloosQRC350DE
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = CloosQRC350DE_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'CloosQRC350DE_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
JaD_transl=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
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
	JaD_transl = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->22)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t4 = sin(qJ(1));
	t5 = qJD(1) * t4;
	unknown(1,1) = -t2 * r_i_i_C(1) - t5 * r_i_i_C(2);
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = t5 * r_i_i_C(1) - t2 * r_i_i_C(2);
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
	JaD_transl = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (17->17), mult. (60->32), div. (0->0), fcn. (38->4), ass. (0->30)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(2));
	t4 = t3 * r_i_i_C(1);
	t6 = sin(qJ(1));
	t7 = t6 * qJD(2);
	t8 = cos(qJ(2));
	t9 = t8 * r_i_i_C(1);
	t11 = t8 * r_i_i_C(2);
	t13 = t3 * r_i_i_C(2);
	t15 = qJD(1) * t6;
	t20 = t1 * qJD(2);
	unknown(1,1) = -t2 * pkin(2) - t15 * r_i_i_C(3) - t2 * t11 + t7 * t13 - t2 * t4 - t7 * t9;
	unknown(1,2) = -t20 * t11 + t15 * t13 - t15 * t9 - t20 * t4;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = t15 * pkin(2) - t2 * r_i_i_C(3) + t15 * t11 + t20 * t13 + t15 * t4 - t20 * t9;
	unknown(2,2) = t7 * t11 + t2 * t13 - t2 * t9 + t7 * t4;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -qJD(2) * t8 * r_i_i_C(1) + qJD(2) * t3 * r_i_i_C(2);
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	JaD_transl = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (77->34), mult. (110->44), div. (0->0), fcn. (71->6), ass. (0->49)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = qJ(2) + qJ(3);
	t4 = sin(t3);
	t5 = t4 * r_i_i_C(1);
	t7 = sin(qJ(1));
	t8 = qJD(2) + qJD(3);
	t9 = t7 * t8;
	t10 = cos(t3);
	t11 = t10 * r_i_i_C(1);
	t13 = t10 * r_i_i_C(2);
	t15 = t4 * r_i_i_C(2);
	t17 = qJD(1) * t7;
	t19 = sin(qJ(2));
	t20 = t19 * pkin(3);
	t21 = t20 + pkin(2);
	t23 = t7 * qJD(2);
	t24 = cos(qJ(2));
	t25 = t24 * pkin(3);
	t28 = t17 * t11;
	t29 = t1 * t8;
	t30 = t29 * t5;
	t31 = t17 * t15;
	t32 = t29 * t13;
	t34 = t1 * qJD(2);
	t46 = t2 * t11;
	t47 = t9 * t5;
	t48 = t2 * t15;
	t49 = t9 * t13;
	t55 = t8 * t10 * r_i_i_C(1);
	t57 = t8 * t4 * r_i_i_C(2);
	unknown(1,1) = -t17 * r_i_i_C(3) - t9 * t11 - t2 * t13 + t9 * t15 - t2 * t21 - t2 * t5 - t23 * t25;
	unknown(1,2) = -t17 * t25 - t34 * t20 - t28 - t30 + t31 - t32;
	unknown(1,3) = -t28 - t30 + t31 - t32;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t2 * r_i_i_C(3) - t29 * t11 + t17 * t13 + t29 * t15 + t17 * t21 + t17 * t5 - t34 * t25;
	unknown(2,2) = -t2 * t25 + t23 * t20 - t46 + t47 + t48 + t49;
	unknown(2,3) = -t46 + t47 + t48 + t49;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -qJD(2) * t24 * pkin(3) - t55 + t57;
	unknown(3,3) = -t55 + t57;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	JaD_transl = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:15
	% EndTime: 2020-06-19 21:40:15
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (266->103), mult. (390->125), div. (0->0), fcn. (295->8), ass. (0->108)
	unknown=NaN(3,6);
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
	t24 = t4 * t15;
	t26 = t11 * t15;
	t28 = qJD(4) * t5;
	t32 = t10 * t26 + t14 * t28 + t15 * t20 + t18 * t5 + t2 * t24;
	t34 = t11 * r_i_i_C(3);
	t36 = t4 * r_i_i_C(3);
	t38 = t4 * pkin(4);
	t40 = t11 * pkin(4);
	t42 = t11 * pkin(5);
	t44 = t4 * pkin(5);
	t46 = sin(qJ(2));
	t47 = t46 * pkin(3);
	t48 = t47 + pkin(2);
	t50 = t8 * qJD(2);
	t51 = cos(qJ(2));
	t52 = t51 * pkin(3);
	t55 = t12 * r_i_i_C(1);
	t56 = t18 * t55;
	t57 = t1 * t9;
	t58 = t6 * r_i_i_C(1);
	t59 = t57 * t58;
	t60 = t1 * t11;
	t61 = t16 * r_i_i_C(1);
	t62 = t60 * t61;
	t63 = t26 * r_i_i_C(2);
	t64 = t18 * t63;
	t65 = t24 * r_i_i_C(2);
	t66 = t57 * t65;
	t67 = t28 * r_i_i_C(2);
	t68 = t60 * t67;
	t69 = t18 * t36;
	t70 = t57 * t34;
	t71 = t18 * t40;
	t72 = t57 * t38;
	t73 = t18 * t44;
	t74 = t57 * t42;
	t76 = t1 * qJD(2);
	t78 = -t18 * t52 - t47 * t76 - t56 - t59 - t62 + t64 + t66 - t68 + t69 - t70 - t71 - t72 + t73 - t74;
	t79 = -t56 - t59 - t62 + t64 + t66 - t68 + t69 - t70 - t71 - t72 + t73 - t74;
	t82 = t1 * t4;
	t85 = t8 * qJD(4);
	t87 = t15 * t85 + t18 * t24 - t2 * t5 - t26 * t57 - t28 * t82;
	t94 = -t12 * t57 + t15 * t2 + t16 * t82 + t18 * t6 + t5 * t85;
	t108 = t2 * t55;
	t109 = t10 * t58;
	t110 = t8 * t11;
	t111 = t110 * t61;
	t112 = t2 * t63;
	t113 = t10 * t65;
	t114 = t110 * t67;
	t115 = t2 * t36;
	t116 = t10 * t34;
	t117 = t2 * t40;
	t118 = t10 * t38;
	t119 = t2 * t44;
	t120 = t10 * t42;
	t123 = -t2 * t52 + t47 * t50 - t108 + t109 + t111 + t112 - t113 + t114 + t115 + t116 - t117 + t118 + t119 + t120;
	t124 = -t108 + t109 + t111 + t112 - t113 + t114 + t115 + t116 - t117 + t118 + t119 + t120;
	t128 = t9 * t11;
	t129 = t5 * r_i_i_C(1);
	t130 = t128 * t129;
	t131 = t4 * qJD(4);
	t132 = t15 * r_i_i_C(1);
	t133 = t131 * t132;
	t134 = t15 * r_i_i_C(2);
	t135 = t128 * t134;
	t136 = t5 * r_i_i_C(2);
	t137 = t131 * t136;
	t138 = t9 * t4;
	t139 = t138 * r_i_i_C(3);
	t140 = t128 * pkin(4);
	t141 = t138 * pkin(5);
	t147 = t11 * qJD(4);
	unknown(1,1) = r_i_i_C(1) * t22 + r_i_i_C(2) * t32 + t10 * t36 - t10 * t40 + t10 * t44 - t2 * t34 - t2 * t38 - t2 * t42 - t2 * t48 - t50 * t52;
	unknown(1,2) = t78;
	unknown(1,3) = t79;
	unknown(1,4) = r_i_i_C(1) * t87 + r_i_i_C(2) * t94;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = r_i_i_C(1) * t94 - r_i_i_C(2) * t87 + t18 * t34 + t18 * t38 + t18 * t42 + t18 * t48 + t36 * t57 - t40 * t57 + t44 * t57 - t52 * t76;
	unknown(2,2) = t123;
	unknown(2,3) = t124;
	unknown(2,4) = r_i_i_C(1) * t32 - r_i_i_C(2) * t22;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -pkin(3) * qJD(2) * t51 - t130 + t133 + t135 + t137 + t139 - t140 + t141;
	unknown(3,3) = -t130 + t133 + t135 + t137 + t139 - t140 + t141;
	unknown(3,4) = -t129 * t147 + t132 * t138 + t134 * t147 + t136 * t138;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	JaD_transl = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:15
	% EndTime: 2020-06-19 21:40:15
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (562->169), mult. (860->243), div. (0->0), fcn. (777->10), ass. (0->131)
	unknown=NaN(3,6);
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
	t42 = t11 * t23;
	t44 = t4 * t23;
	t46 = qJD(5) * t29;
	t48 = t10 * t44 - t2 * t42 - t22 * t29 - t28 * t23 + t35 * t46;
	t50 = t4 * t15;
	t52 = t11 * t15;
	t54 = qJD(4) * t5;
	t58 = -t10 * t52 - t14 * t54 - t20 * t15 - t18 * t5 - t2 * t50;
	t60 = t4 * pkin(4);
	t62 = t11 * pkin(4);
	t64 = t11 * pkin(5);
	t66 = t4 * pkin(5);
	t68 = sin(qJ(2));
	t69 = t68 * pkin(3);
	t70 = t69 + pkin(2);
	t72 = t8 * qJD(2);
	t73 = cos(qJ(2));
	t74 = t73 * pkin(3);
	t77 = t12 * t23;
	t79 = t1 * t9;
	t80 = t6 * t23;
	t82 = t1 * t11;
	t83 = t16 * t23;
	t85 = t5 * qJD(5);
	t86 = t85 * t29;
	t90 = t1 * t4;
	t93 = (t18 * t33 - t18 * t77 - t79 * t31 - t90 * t36 - t79 * t80 - t82 * t83 - t82 * t86) * r_i_i_C(1);
	t94 = t12 * t29;
	t96 = t6 * t29;
	t98 = t16 * t29;
	t100 = t85 * t23;
	t106 = (-t82 * t100 + t18 * t44 + t18 * t94 - t79 * t42 + t90 * t46 + t79 * t96 + t82 * t98) * r_i_i_C(2);
	t107 = t52 * r_i_i_C(3);
	t108 = t18 * t107;
	t109 = t50 * r_i_i_C(3);
	t110 = t79 * t109;
	t111 = t54 * r_i_i_C(3);
	t112 = t82 * t111;
	t113 = t18 * t62;
	t114 = t79 * t60;
	t115 = t18 * t66;
	t116 = t79 * t64;
	t118 = t1 * qJD(2);
	t120 = -t118 * t69 - t18 * t74 + t106 - t108 - t110 + t112 - t113 - t114 + t115 - t116 + t93;
	t126 = t8 * qJD(4);
	t128 = t126 * t15 + t18 * t50 - t2 * t5 - t79 * t52 - t90 * t54;
	t134 = (-t90 * t15 - t8 * t5) * qJD(5);
	t135 = t29 * r_i_i_C(1);
	t139 = t23 * r_i_i_C(2);
	t146 = t79 * t12 - t126 * t5 - t2 * t15 - t90 * t16 - t18 * t6;
	t152 = -t8 * t15 + t90 * t5;
	t153 = t152 * qJD(5);
	t155 = t18 * t42;
	t156 = t79 * t44;
	t157 = t82 * t46;
	t162 = t18 * t31;
	t163 = t79 * t33;
	t164 = t82 * t36;
	t169 = -t152 * qJD(5);
	t193 = (t10 * t31 + t10 * t80 + t14 * t36 + t2 * t33 - t2 * t77 + t35 * t83 + t35 * t86) * r_i_i_C(1);
	t202 = (t10 * t42 - t10 * t96 + t35 * t100 - t14 * t46 + t2 * t44 + t2 * t94 - t35 * t98) * r_i_i_C(2);
	t203 = t2 * t107;
	t204 = t10 * t109;
	t205 = t35 * t111;
	t206 = t2 * t62;
	t207 = t10 * t60;
	t208 = t2 * t66;
	t209 = t10 * t64;
	t212 = -t2 * t74 + t72 * t69 + t193 + t202 - t203 + t204 - t205 - t206 + t207 + t208 + t209;
	t219 = (-t1 * t5 + t14 * t15) * qJD(5);
	t229 = t9 * t11;
	t230 = t5 * t23;
	t232 = t4 * qJD(4);
	t233 = t23 * t15;
	t236 = t9 * t4;
	t238 = t11 * qJD(5);
	t241 = (-t229 * t230 - t238 * t23 + t232 * t233 + t236 * t29 + t6 * t46) * r_i_i_C(1);
	t242 = t5 * t29;
	t244 = t15 * t29;
	t250 = (t229 * t242 + t236 * t23 - t232 * t244 + t238 * t29 + t6 * t36) * r_i_i_C(2);
	t251 = t15 * r_i_i_C(3);
	t252 = t229 * t251;
	t253 = t5 * r_i_i_C(3);
	t254 = t232 * t253;
	t255 = t229 * pkin(4);
	t256 = t236 * pkin(5);
	t263 = t11 * qJD(4);
	t281 = t4 * qJD(5);
	unknown(1,1) = t38 * r_i_i_C(1) + t48 * r_i_i_C(2) + t58 * r_i_i_C(3) - t10 * t62 + t10 * t66 - t2 * t60 - t2 * t64 - t2 * t70 - t72 * t74;
	unknown(1,2) = t120;
	unknown(1,3) = t93 + t106 - t108 - t110 + t112 - t113 - t114 + t115 - t116;
	unknown(1,4) = t128 * t23 * r_i_i_C(1) - t128 * t29 * r_i_i_C(2) + t146 * r_i_i_C(3) - t134 * t135 - t134 * t139;
	unknown(1,5) = (-t146 * t29 - t153 * t23 - t155 - t156 - t157) * r_i_i_C(1) + (-t146 * t23 + t153 * t29 + t162 + t163 - t164) * r_i_i_C(2);
	unknown(1,6) = 0.0e0;
	unknown(2,1) = (-t146 * t23 - t169 * t29 + t162 + t163 - t164) * r_i_i_C(1) + (t146 * t29 - t169 * t23 + t155 + t156 + t157) * r_i_i_C(2) + t128 * r_i_i_C(3) + t18 * t60 - t79 * t62 + t18 * t64 + t79 * t66 + t18 * t70 - t118 * t74;
	unknown(2,2) = t212;
	unknown(2,3) = t193 + t202 - t203 + t204 - t205 - t206 + t207 + t208 + t209;
	unknown(2,4) = -t58 * t23 * r_i_i_C(1) + t58 * t29 * r_i_i_C(2) + t22 * r_i_i_C(3) - t219 * t135 - t219 * t139;
	unknown(2,5) = t48 * r_i_i_C(1) - t38 * r_i_i_C(2);
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -qJD(2) * t73 * pkin(3) + t241 + t250 - t252 - t254 - t255 + t256;
	unknown(3,3) = t241 + t250 - t252 - t254 - t255 + t256;
	unknown(3,4) = -t263 * t230 * r_i_i_C(1) + t236 * t233 * r_i_i_C(1) + t52 * t46 * r_i_i_C(1) - t236 * t244 * r_i_i_C(2) + t263 * t242 * r_i_i_C(2) + t52 * t36 * r_i_i_C(2) - t236 * t253 - t263 * t251;
	unknown(3,5) = (-t12 * t36 - t229 * t23 + t236 * t242 + t263 * t244 + t281 * t29) * r_i_i_C(1) + (t12 * t46 + t229 * t29 + t281 * t23 + t236 * t230 + t263 * t233) * r_i_i_C(2);
	unknown(3,6) = 0.0e0;
	JaD_transl = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:16
	% EndTime: 2020-06-19 21:40:16
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (1852->311), mult. (2762->469), div. (0->0), fcn. (2475->12), ass. (0->199)
	unknown=NaN(3,6);
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
	t45 = t27 * t23 - t35 * t29;
	t47 = pkin(7) * qJD(5) - qJD(6);
	t48 = t45 * t47;
	t49 = sin(t40);
	t51 = t4 * t15;
	t53 = t11 * t15;
	t55 = qJD(4) * t5;
	t59 = -t10 * t53 - t14 * t55 - t20 * t15 - t18 * t5 - t2 * t51;
	t63 = t1 * t5 - t14 * t15;
	t64 = t63 * t47;
	t66 = -t38 * t41 - t64 * t41 + t48 * t49 - t59 * t49;
	t72 = -t38 * t49 - t48 * t41 + t59 * t41 - t64 * t49;
	t76 = t11 * t23;
	t78 = t4 * t23;
	t80 = qJD(5) * t29;
	t82 = t10 * t78 - t2 * t76 - t22 * t29 - t28 * t23 + t35 * t80;
	t85 = t4 * pkin(4);
	t87 = t11 * pkin(4);
	t89 = t11 * pkin(5);
	t91 = t4 * pkin(5);
	t93 = sin(qJ(2));
	t94 = t93 * pkin(3);
	t95 = t94 + pkin(2);
	t97 = t8 * qJD(2);
	t98 = cos(qJ(2));
	t99 = t98 * pkin(3);
	t102 = t12 * t23;
	t104 = t1 * t9;
	t105 = t6 * t23;
	t107 = t1 * t11;
	t108 = t16 * t23;
	t110 = t5 * qJD(5);
	t111 = t110 * t29;
	t115 = t1 * t4;
	t117 = -t18 * t102 - t104 * t105 - t104 * t31 - t107 * t108 - t107 * t111 - t115 * t36 + t18 * t33;
	t119 = t5 * t23;
	t123 = (t107 * t119 - t115 * t29) * t47;
	t125 = t53 * t49;
	t127 = t51 * t49;
	t129 = t55 * t49;
	t131 = t15 * t47;
	t132 = t131 * t41;
	t135 = (t104 * t127 - t107 * t129 - t107 * t132 - t117 * t41 + t123 * t49 + t18 * t125) * r_i_i_C(1);
	t138 = t53 * t41;
	t140 = t51 * t41;
	t142 = t55 * t41;
	t144 = t131 * t49;
	t147 = (-t104 * t140 + t107 * t142 - t107 * t144 - t117 * t49 - t123 * t41 - t18 * t138) * r_i_i_C(2);
	t148 = t12 * t29;
	t150 = t6 * t29;
	t152 = t16 * t29;
	t154 = t110 * t23;
	t159 = t104 * t150 - t104 * t76 + t107 * t152 - t107 * t154 + t115 * t80 + t18 * t148 + t18 * t78;
	t160 = t159 * r_i_i_C(3);
	t161 = t159 * pkin(6);
	t162 = t18 * t87;
	t163 = t104 * t85;
	t164 = t18 * t91;
	t165 = t104 * t89;
	t167 = t1 * qJD(2);
	t175 = t8 * qJD(4);
	t177 = -t104 * t53 - t115 * t55 + t175 * t15 + t18 * t51 - t2 * t5;
	t178 = t177 * t23;
	t182 = -t115 * t15 - t8 * t5;
	t183 = t182 * qJD(5);
	t184 = t29 * t41;
	t186 = t182 * t23;
	t187 = t47 * t49;
	t194 = t104 * t12 - t115 * t16 - t2 * t15 - t175 * t5 - t18 * t6;
	t198 = t115 * t5 - t8 * t15;
	t199 = t198 * t47;
	t204 = t29 * t49;
	t206 = t47 * t41;
	t212 = t177 * t29;
	t214 = t23 * r_i_i_C(3);
	t217 = t23 * pkin(6);
	t221 = t198 * qJD(5);
	t223 = t18 * t76;
	t224 = t104 * t78;
	t225 = t107 * t80;
	t226 = -t194 * t29 - t221 * t23 - t223 - t224 - t225;
	t231 = (t107 * t23 - t198 * t29) * t47;
	t235 = t18 * t31;
	t236 = t104 * t33;
	t237 = t107 * t36;
	t238 = t194 * t23 - t221 * t29 - t235 - t236 + t237;
	t239 = t238 * pkin(7);
	t242 = t107 * t29;
	t243 = t198 * t23 + t242;
	t244 = t243 * pkin(7);
	t246 = -t177 * pkin(7);
	t248 = -t182 * pkin(7);
	t264 = t243 * t47;
	t267 = -t182 * t47;
	t279 = -t198 * qJD(5);
	t281 = -t194 * t23 - t279 * t29 + t235 + t236 - t237;
	t285 = (-t198 * t23 - t242) * t47;
	t288 = t182 * t47;
	t300 = t194 * t29 - t279 * t23 + t223 + t224 + t225;
	t317 = t10 * t105 + t10 * t31 - t2 * t102 + t35 * t108 + t35 * t111 + t14 * t36 + t2 * t33;
	t322 = (-t35 * t119 + t14 * t29) * t47;
	t329 = (-t10 * t127 + t2 * t125 + t35 * t129 + t35 * t132 - t317 * t41 + t322 * t49) * r_i_i_C(1);
	t337 = (t10 * t140 - t2 * t138 - t35 * t142 + t35 * t144 - t317 * t49 - t322 * t41) * r_i_i_C(2);
	t345 = -t10 * t150 + t10 * t76 - t14 * t80 + t2 * t148 - t35 * t152 + t35 * t154 + t2 * t78;
	t346 = t345 * r_i_i_C(3);
	t347 = t345 * pkin(6);
	t348 = t2 * t87;
	t349 = t10 * t85;
	t350 = t2 * t91;
	t351 = t10 * t89;
	t356 = -t59 * t23;
	t358 = -t63 * qJD(5);
	t360 = -t63 * t23;
	t363 = t27 * t47;
	t374 = -t59 * t29;
	t384 = (-t35 * t23 - t27 * t29) * t47;
	t386 = t38 * pkin(7);
	t388 = t45 * pkin(7);
	t390 = t59 * pkin(7);
	t392 = t63 * pkin(7);
	t410 = t9 * t11;
	t412 = t4 * qJD(4);
	t413 = t15 * t23;
	t416 = t9 * t4;
	t418 = t11 * qJD(5);
	t420 = -t410 * t119 - t418 * t23 + t416 * t29 + t412 * t413 + t6 * t80;
	t423 = (-t105 - t31) * t47;
	t425 = t15 * t49;
	t427 = t5 * t49;
	t431 = (t51 * t206 - t420 * t41 + t410 * t425 + t412 * t427 + t423 * t49) * r_i_i_C(1);
	t434 = t15 * t41;
	t436 = t5 * t41;
	t440 = (t51 * t187 - t423 * t41 - t410 * t434 - t412 * t436 - t420 * t49) * r_i_i_C(2);
	t441 = t5 * t29;
	t443 = t15 * t29;
	t448 = t416 * t23 + t418 * t29 + t6 * t36 + t410 * t441 - t412 * t443;
	t449 = t448 * r_i_i_C(3);
	t450 = t448 * pkin(6);
	t451 = t410 * pkin(4);
	t452 = t416 * pkin(5);
	t459 = t11 * qJD(4);
	t464 = t23 * t47;
	t502 = t4 * qJD(5);
	t504 = -t12 * t36 - t410 * t23 + t502 * t29 + t416 * t441 + t459 * t443;
	t507 = (-t148 - t78) * t47;
	t514 = -t416 * t119 - t12 * t80 - t502 * t23 - t410 * t29 - t459 * t413;
	t515 = t514 * pkin(7);
	t517 = t102 - t33;
	t518 = t517 * pkin(7);
	t520 = t15 * pkin(7);
	t523 = t5 * pkin(7);
	t526 = pkin(7) * t47;
	t547 = t517 * t47;
	unknown(1,1) = t82 * pkin(6) + t66 * r_i_i_C(1) + t72 * r_i_i_C(2) + t82 * r_i_i_C(3) - t10 * t87 + t10 * t91 - t2 * t85 - t2 * t89 - t2 * t95 - t97 * t99;
	unknown(1,2) = -t167 * t94 - t18 * t99 + t135 + t147 + t160 + t161 - t162 - t163 + t164 - t165;
	unknown(1,3) = t135 + t147 + t160 + t161 - t162 - t163 + t164 - t165;
	unknown(1,4) = (-t178 * t41 + t183 * t184 + t186 * t187 - t194 * t49 - t199 * t41) * r_i_i_C(1) + (-t178 * t49 + t183 * t204 - t186 * t206 + t194 * t41 - t199 * t49) * r_i_i_C(2) - t212 * r_i_i_C(3) - t183 * t214 - t212 * pkin(6) - t183 * t217;
	unknown(1,5) = (t248 * t187 + t244 * t206 - t226 * t41 + t231 * t49 + t239 * t49 - t246 * t41) * r_i_i_C(1) + (t244 * t187 - t248 * t206 - t226 * t49 - t231 * t41 - t239 * t41 - t246 * t49) * r_i_i_C(2) - t238 * r_i_i_C(3) - t238 * pkin(6);
	unknown(1,6) = (-t177 * t41 - t238 * t49 - t264 * t41 - t267 * t49) * r_i_i_C(1) + (-t177 * t49 + t238 * t41 - t264 * t49 + t267 * t41) * r_i_i_C(2);
	unknown(2,1) = (-t177 * t49 - t281 * t41 + t285 * t49 - t288 * t41) * r_i_i_C(1) + (t177 * t41 - t281 * t49 - t285 * t41 - t288 * t49) * r_i_i_C(2) + t300 * r_i_i_C(3) + t300 * pkin(6) + t18 * t85 - t104 * t87 + t18 * t89 + t104 * t91 + t18 * t95 - t167 * t99;
	unknown(2,2) = -t2 * t99 + t97 * t94 + t329 + t337 + t346 + t347 - t348 + t349 + t350 + t351;
	unknown(2,3) = t329 + t337 + t346 + t347 - t348 + t349 + t350 + t351;
	unknown(2,4) = (t358 * t184 + t360 * t187 - t22 * t49 - t356 * t41 - t363 * t41) * r_i_i_C(1) + (t358 * t204 - t360 * t206 + t22 * t41 - t356 * t49 - t363 * t49) * r_i_i_C(2) - t374 * r_i_i_C(3) - t358 * t214 - t374 * pkin(6) - t358 * t217;
	unknown(2,5) = (t392 * t187 + t388 * t206 + t384 * t49 + t386 * t49 - t390 * t41 - t82 * t41) * r_i_i_C(1) + (t388 * t187 - t392 * t206 - t384 * t41 - t386 * t41 - t390 * t49 - t82 * t49) * r_i_i_C(2) - t38 * r_i_i_C(3) - t38 * pkin(6);
	unknown(2,6) = t72 * r_i_i_C(1) - t66 * r_i_i_C(2);
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -qJD(2) * t98 * pkin(3) + t431 + t440 + t449 + t450 - t451 + t452;
	unknown(3,3) = t431 + t440 + t449 + t450 - t451 + t452;
	unknown(3,4) = (t459 * t119 * t41 - t416 * t413 * t41 - t53 * t80 * t41 - t53 * t464 * t49 - t12 * t206 + t416 * t427 + t459 * t425) * r_i_i_C(1) + (t459 * t119 * t49 + t53 * t464 * t41 - t416 * t413 * t49 - t53 * t80 * t49 - t12 * t187 - t416 * t436 - t459 * t434) * r_i_i_C(2) - t416 * t443 * r_i_i_C(3) + t459 * t441 * r_i_i_C(3) + t53 * t36 * r_i_i_C(3) - t416 * t443 * pkin(6) + t459 * t441 * pkin(6) + t53 * t36 * pkin(6);
	unknown(3,5) = (t416 * t520 * t41 - t459 * t523 * t41 + t53 * t526 * t49 + t518 * t206 - t504 * t41 + t507 * t49 + t515 * t49) * r_i_i_C(1) + (-t53 * t526 * t41 + t416 * t520 * t49 - t459 * t523 * t49 + t518 * t187 - t507 * t41 - t515 * t41 - t504 * t49) * r_i_i_C(2) - t514 * r_i_i_C(3) - t514 * pkin(6);
	unknown(3,6) = (-t53 * t187 - t547 * t41 - t416 * t434 + t459 * t436 - t514 * t49) * r_i_i_C(1) + (t53 * t206 + t514 * t41 - t416 * t425 + t459 * t427 - t547 * t49) * r_i_i_C(2);
	JaD_transl = unknown;
end