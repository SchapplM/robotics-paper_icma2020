% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% CloosQRC350DE
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = CloosQRC350DE_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'CloosQRC350DE_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350DE_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
Ja_transl=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.10s
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
	Ja_transl = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->20)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t3 = cos(qJ(1));
	unknown(1,1) = -t1 * r_i_i_C(1) + t3 * r_i_i_C(2);
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
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
	Ja_transl = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (9->9), mult. (22->18), div. (0->0), fcn. (22->4), ass. (0->26)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t5 = cos(qJ(2));
	t6 = t1 * t5;
	t8 = cos(qJ(1));
	t12 = t8 * t5;
	t14 = t8 * t2;
	unknown(1,1) = -t1 * pkin(2) - t3 * r_i_i_C(1) - t6 * r_i_i_C(2) + t8 * r_i_i_C(3);
	unknown(1,2) = t12 * r_i_i_C(1) - t14 * r_i_i_C(2);
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t8 * pkin(2) - t14 * r_i_i_C(1) - t12 * r_i_i_C(2) - t1 * r_i_i_C(3);
	unknown(2,2) = -t6 * r_i_i_C(1) + t3 * r_i_i_C(2);
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t2 * r_i_i_C(1) - t5 * r_i_i_C(2);
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (33->17), mult. (39->23), div. (0->0), fcn. (39->6), ass. (0->37)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = sin(t2);
	t4 = t1 * t3;
	t6 = cos(t2);
	t7 = t1 * t6;
	t9 = cos(qJ(1));
	t11 = sin(qJ(2));
	t12 = t11 * pkin(3);
	t13 = t12 + pkin(2);
	t16 = t9 * t6;
	t17 = t16 * r_i_i_C(1);
	t18 = t9 * t3;
	t19 = t18 * r_i_i_C(2);
	t20 = cos(qJ(2));
	t30 = t7 * r_i_i_C(1);
	t31 = t4 * r_i_i_C(2);
	t36 = t3 * r_i_i_C(1);
	t37 = t6 * r_i_i_C(2);
	unknown(1,1) = -t4 * r_i_i_C(1) - t7 * r_i_i_C(2) + t9 * r_i_i_C(3) - t1 * t13;
	unknown(1,2) = t9 * t20 * pkin(3) + t17 - t19;
	unknown(1,3) = t17 - t19;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t18 * r_i_i_C(1) - t16 * r_i_i_C(2) - t1 * r_i_i_C(3) - t9 * t13;
	unknown(2,2) = -t1 * t20 * pkin(3) - t30 + t31;
	unknown(2,3) = -t30 + t31;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t36 - t37 - t12;
	unknown(3,3) = -t36 - t37;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (96->46), mult. (119->56), div. (0->0), fcn. (127->8), ass. (0->54)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = qJ(2) + qJ(3);
	t3 = sin(t2);
	t4 = t1 * t3;
	t5 = cos(qJ(4));
	t7 = cos(qJ(1));
	t8 = sin(qJ(4));
	t10 = -t4 * t5 - t7 * t8;
	t14 = t4 * t8 - t7 * t5;
	t16 = cos(t2);
	t17 = t1 * t16;
	t21 = sin(qJ(2));
	t22 = t21 * pkin(3);
	t23 = t22 + pkin(2);
	t26 = t7 * t16;
	t27 = t5 * r_i_i_C(1);
	t28 = t26 * t27;
	t29 = t8 * r_i_i_C(2);
	t30 = t26 * t29;
	t31 = t7 * t3;
	t32 = t31 * r_i_i_C(3);
	t33 = t26 * pkin(4);
	t34 = t31 * pkin(5);
	t35 = cos(qJ(2));
	t42 = -t1 * t5 - t31 * t8;
	t46 = t1 * t8 - t31 * t5;
	t56 = t17 * t27;
	t57 = t17 * t29;
	t58 = t4 * r_i_i_C(3);
	t59 = t17 * pkin(4);
	t60 = t4 * pkin(5);
	t69 = t3 * t5 * r_i_i_C(1);
	t71 = t3 * t8 * r_i_i_C(2);
	t72 = t16 * r_i_i_C(3);
	t73 = t3 * pkin(4);
	t74 = t16 * pkin(5);
	unknown(1,1) = -t4 * pkin(4) - t17 * pkin(5) + t10 * r_i_i_C(1) + t14 * r_i_i_C(2) - t17 * r_i_i_C(3) - t1 * t23;
	unknown(1,2) = t7 * t35 * pkin(3) + t28 - t30 - t32 + t33 - t34;
	unknown(1,3) = t28 - t30 - t32 + t33 - t34;
	unknown(1,4) = t42 * r_i_i_C(1) + t46 * r_i_i_C(2);
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(2,1) = -t31 * pkin(4) - t26 * pkin(5) + t46 * r_i_i_C(1) - t42 * r_i_i_C(2) - t26 * r_i_i_C(3) - t7 * t23;
	unknown(2,2) = -t1 * t35 * pkin(3) - t56 + t57 + t58 - t59 + t60;
	unknown(2,3) = -t56 + t57 + t58 - t59 + t60;
	unknown(2,4) = t14 * r_i_i_C(1) - t10 * r_i_i_C(2);
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t69 + t71 - t72 - t73 - t74 - t22;
	unknown(3,3) = -t69 + t71 - t72 - t73 - t74;
	unknown(3,4) = -t16 * t8 * r_i_i_C(1) - t16 * t5 * r_i_i_C(2);
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (163->66), mult. (225->98), div. (0->0), fcn. (263->10), ass. (0->64)
	unknown=NaN(3,6);
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
	t21 = -t10 * t15 - t14 * t11;
	t25 = -t4 * t8 + t7 * t5;
	t29 = sin(qJ(2));
	t30 = t29 * pkin(3);
	t31 = t30 + pkin(2);
	t34 = t7 * t13;
	t35 = t5 * t11;
	t37 = t7 * t3;
	t40 = (-t37 * t15 + t34 * t35) * r_i_i_C(1);
	t41 = t5 * t15;
	t45 = (-t37 * t11 - t34 * t41) * r_i_i_C(2);
	t46 = t8 * r_i_i_C(3);
	t47 = t34 * t46;
	t48 = t34 * pkin(4);
	t49 = t37 * pkin(5);
	t50 = cos(qJ(2));
	t57 = -t1 * t5 - t37 * t8;
	t64 = -t1 * t8 + t37 * t5;
	t68 = t34 * t11;
	t72 = t34 * t15;
	t90 = (-t14 * t35 + t4 * t15) * r_i_i_C(1);
	t94 = (t4 * t11 + t14 * t41) * r_i_i_C(2);
	t95 = t14 * t46;
	t96 = t14 * pkin(4);
	t97 = t4 * pkin(5);
	t111 = t3 * t5;
	t115 = (-t111 * t11 - t13 * t15) * r_i_i_C(1);
	t119 = (-t13 * t11 + t111 * t15) * r_i_i_C(2);
	t121 = t3 * t8 * r_i_i_C(3);
	t122 = t3 * pkin(4);
	t123 = t13 * pkin(5);
	t126 = t13 * t8;
	t131 = t13 * t5;
	unknown(1,1) = -t4 * pkin(4) - t14 * pkin(5) + t17 * r_i_i_C(1) + t21 * r_i_i_C(2) + t25 * r_i_i_C(3) - t1 * t31;
	unknown(1,2) = t7 * t50 * pkin(3) + t40 + t45 + t47 + t48 - t49;
	unknown(1,3) = t40 + t45 + t47 + t48 - t49;
	unknown(1,4) = t57 * t11 * r_i_i_C(1) - t57 * t15 * r_i_i_C(2) + t64 * r_i_i_C(3);
	unknown(1,5) = (-t64 * t15 + t68) * r_i_i_C(1) + (-t64 * t11 - t72) * r_i_i_C(2);
	unknown(1,6) = 0.0e0;
	unknown(2,1) = (-t64 * t11 - t72) * r_i_i_C(1) + (t64 * t15 - t68) * r_i_i_C(2) + t57 * r_i_i_C(3) - t37 * pkin(4) - t34 * pkin(5) - t7 * t31;
	unknown(2,2) = -t1 * t50 * pkin(3) + t90 + t94 - t95 - t96 + t97;
	unknown(2,3) = t90 + t94 - t95 - t96 + t97;
	unknown(2,4) = -t25 * t11 * r_i_i_C(1) + t25 * t15 * r_i_i_C(2) + t10 * r_i_i_C(3);
	unknown(2,5) = t21 * r_i_i_C(1) - t17 * r_i_i_C(2);
	unknown(2,6) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t115 + t119 - t121 - t122 - t123 - t30;
	unknown(3,3) = t115 + t119 - t121 - t122 - t123;
	unknown(3,4) = -t126 * t11 * r_i_i_C(1) + t126 * t15 * r_i_i_C(2) + t131 * r_i_i_C(3);
	unknown(3,5) = (-t3 * t11 - t131 * t15) * r_i_i_C(1) + (-t131 * t11 + t3 * t15) * r_i_i_C(2);
	unknown(3,6) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-06-19 21:40:14
	% EndTime: 2020-06-19 21:40:14
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (460->122), mult. (645->187), div. (0->0), fcn. (703->12), ass. (0->95)
	unknown=NaN(3,6);
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
	t31 = -t17 * t25 + t24 * t20;
	t35 = -t10 * t15 - t14 * t11;
	t40 = sin(qJ(2));
	t41 = t40 * pkin(3);
	t42 = t41 + pkin(2);
	t45 = t7 * t13;
	t46 = t5 * t11;
	t48 = t7 * t3;
	t50 = -t48 * t15 + t45 * t46;
	t52 = t8 * t25;
	t55 = (-t50 * t20 - t45 * t52) * r_i_i_C(1);
	t57 = t8 * t20;
	t60 = (-t50 * t25 + t45 * t57) * r_i_i_C(2);
	t61 = t5 * t15;
	t64 = -t48 * t11 - t45 * t61;
	t65 = t64 * r_i_i_C(3);
	t66 = t64 * pkin(6);
	t67 = t45 * pkin(4);
	t68 = t48 * pkin(5);
	t69 = cos(qJ(2));
	t76 = -t1 * t5 - t48 * t8;
	t77 = t76 * t11;
	t81 = -t1 * t8 + t48 * t5;
	t89 = t76 * t15;
	t94 = t45 * t11;
	t95 = -t81 * t15 + t94;
	t98 = t45 * t15;
	t99 = t81 * t11 + t98;
	t100 = t99 * pkin(7);
	t102 = -t76 * pkin(7);
	t124 = -t81 * t11 - t98;
	t134 = t81 * t15 - t94;
	t143 = -t14 * t46 + t4 * t15;
	t147 = (t14 * t52 - t143 * t20) * r_i_i_C(1);
	t151 = (-t14 * t57 - t143 * t25) * r_i_i_C(2);
	t154 = t4 * t11 + t14 * t61;
	t155 = t154 * r_i_i_C(3);
	t156 = t154 * pkin(6);
	t157 = t14 * pkin(4);
	t158 = t4 * pkin(5);
	t163 = -t24 * t11;
	t172 = -t24 * t15;
	t177 = t17 * pkin(7);
	t179 = t24 * pkin(7);
	t194 = t3 * t5;
	t197 = -t194 * t11 - t13 * t15;
	t199 = t3 * t8;
	t202 = (-t197 * t20 + t199 * t25) * r_i_i_C(1);
	t206 = (-t197 * t25 - t199 * t20) * r_i_i_C(2);
	t209 = -t13 * t11 + t194 * t15;
	t210 = t209 * r_i_i_C(3);
	t211 = t209 * pkin(6);
	t212 = t3 * pkin(4);
	t213 = t13 * pkin(5);
	t216 = t13 * t8;
	t219 = t13 * t5;
	t235 = -t3 * t11 - t219 * t15;
	t239 = t219 * t11 - t3 * t15;
	t240 = t239 * pkin(7);
	unknown(1,1) = -t4 * pkin(4) - t14 * pkin(5) + t35 * pkin(6) + t27 * r_i_i_C(1) + t31 * r_i_i_C(2) + t35 * r_i_i_C(3) - t1 * t42;
	unknown(1,2) = t7 * t69 * pkin(3) + t55 + t60 + t65 + t66 + t67 - t68;
	unknown(1,3) = t55 + t60 + t65 + t66 + t67 - t68;
	unknown(1,4) = (-t77 * t20 - t81 * t25) * r_i_i_C(1) + (t81 * t20 - t77 * t25) * r_i_i_C(2) - t89 * r_i_i_C(3) - t89 * pkin(6);
	unknown(1,5) = (t100 * t25 - t102 * t20 - t95 * t20) * r_i_i_C(1) + (-t100 * t20 - t102 * t25 - t95 * t25) * r_i_i_C(2) - t99 * r_i_i_C(3) - t99 * pkin(6);
	unknown(1,6) = (-t76 * t20 - t99 * t25) * r_i_i_C(1) + (t99 * t20 - t76 * t25) * r_i_i_C(2);
	unknown(2,1) = (-t124 * t20 - t76 * t25) * r_i_i_C(1) + (-t124 * t25 + t76 * t20) * r_i_i_C(2) + t134 * r_i_i_C(3) + t134 * pkin(6) - t48 * pkin(4) - t45 * pkin(5) - t7 * t42;
	unknown(2,2) = -t1 * t69 * pkin(3) + t147 + t151 + t155 + t156 - t157 + t158;
	unknown(2,3) = t147 + t151 + t155 + t156 - t157 + t158;
	unknown(2,4) = (-t10 * t25 - t163 * t20) * r_i_i_C(1) + (t10 * t20 - t163 * t25) * r_i_i_C(2) - t172 * r_i_i_C(3) - t172 * pkin(6);
	unknown(2,5) = (t177 * t25 - t179 * t20 - t35 * t20) * r_i_i_C(1) + (-t177 * t20 - t179 * t25 - t35 * t25) * r_i_i_C(2) - t17 * r_i_i_C(3) - t17 * pkin(6);
	unknown(2,6) = t31 * r_i_i_C(1) - t27 * r_i_i_C(2);
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t202 + t206 + t210 + t211 - t212 - t213 - t41;
	unknown(3,3) = t202 + t206 + t210 + t211 - t212 - t213;
	unknown(3,4) = (t216 * t11 * t20 - t219 * t25) * r_i_i_C(1) + (t216 * t11 * t25 + t219 * t20) * r_i_i_C(2) + t216 * t15 * r_i_i_C(3) + t216 * t15 * pkin(6);
	unknown(3,5) = (-t216 * pkin(7) * t20 - t235 * t20 + t240 * t25) * r_i_i_C(1) + (-t216 * pkin(7) * t25 - t240 * t20 - t235 * t25) * r_i_i_C(2) - t239 * r_i_i_C(3) - t239 * pkin(6);
	unknown(3,6) = (t216 * t20 - t239 * t25) * r_i_i_C(1) + (t239 * t20 + t216 * t25) * r_i_i_C(2);
	Ja_transl = unknown;
end