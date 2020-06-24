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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0, 0; -t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (9->6), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = -pkin(2) + t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [t4 * r_i_i_C(3) + t5 * t2, t7 * t4, 0, 0, 0, 0; -t2 * r_i_i_C(3) + t5 * t4, -t7 * t2, 0, 0, 0, 0; 0, t6, 0, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (33->12), mult. (39->16), div. (0->0), fcn. (39->6), ass. (0->14)
	t6 = qJ(2) + qJ(3);
	t4 = sin(t6);
	t5 = cos(t6);
	t12 = -t4 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t17 = t12 - sin(qJ(2)) * pkin(3);
	t16 = cos(qJ(2)) * pkin(3);
	t15 = r_i_i_C(1) * t5;
	t14 = r_i_i_C(2) * t4;
	t11 = -pkin(2) + t17;
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t2 = t10 * t15;
	t1 = t8 * t14;
	t3 = [t10 * r_i_i_C(3) + t11 * t8, t2 + (-t14 + t16) * t10, -t10 * t14 + t2, 0, 0, 0; -t8 * r_i_i_C(3) + t10 * t11, t1 + (-t15 - t16) * t8, -t15 * t8 + t1, 0, 0, 0; 0, t17, t12, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (96->28), mult. (119->39), div. (0->0), fcn. (127->8), ass. (0->28)
	t14 = qJ(2) + qJ(3);
	t12 = sin(t14);
	t13 = cos(t14);
	t15 = sin(qJ(4));
	t35 = r_i_i_C(2) * t15;
	t39 = pkin(5) + r_i_i_C(3);
	t40 = t12 * t39 + t13 * t35;
	t37 = sin(qJ(2)) * pkin(3);
	t36 = cos(qJ(2)) * pkin(3);
	t20 = cos(qJ(1));
	t18 = cos(qJ(4));
	t29 = t18 * t20;
	t34 = (pkin(4) * t20 + r_i_i_C(1) * t29) * t13;
	t32 = t15 * t20;
	t17 = sin(qJ(1));
	t31 = t17 * t15;
	t30 = t17 * t18;
	t28 = t40 * t17;
	t26 = t39 * t13;
	t25 = -r_i_i_C(1) * t18 - pkin(4);
	t24 = t25 * t13;
	t22 = -t12 * pkin(4) - pkin(2) - t26 - t37;
	t21 = -t26 + (t25 + t35) * t12;
	t4 = -t12 * t29 + t31;
	t3 = t12 * t32 + t30;
	t2 = t12 * t30 + t32;
	t1 = t12 * t31 - t29;
	t5 = [-r_i_i_C(1) * t2 + r_i_i_C(2) * t1 + t22 * t17, (-t40 + t36) * t20 + t34, -t20 * t40 + t34, -r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t22 * t20, (t24 - t36) * t17 + t28, t17 * t24 + t28, r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t21 - t37, t21, (-r_i_i_C(1) * t15 - r_i_i_C(2) * t18) * t13, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (163->46), mult. (225->73), div. (0->0), fcn. (263->10), ass. (0->39)
	t20 = qJ(2) + qJ(3);
	t18 = sin(t20);
	t19 = cos(t20);
	t25 = cos(qJ(5));
	t21 = sin(qJ(5));
	t26 = cos(qJ(4));
	t43 = t21 * t26;
	t30 = t18 * t25 + t19 * t43;
	t40 = t25 * t26;
	t31 = -t18 * t21 + t19 * t40;
	t53 = r_i_i_C(1) * t31 - r_i_i_C(2) * t30;
	t50 = cos(qJ(2)) * pkin(3);
	t49 = pkin(5) * t18;
	t48 = t19 * pkin(5);
	t47 = sin(qJ(2)) * pkin(3);
	t46 = t19 * t21;
	t45 = t19 * t25;
	t28 = cos(qJ(1));
	t44 = t19 * t28;
	t22 = sin(qJ(4));
	t24 = sin(qJ(1));
	t42 = t24 * t22;
	t41 = t24 * t26;
	t39 = t28 * t22;
	t38 = t28 * t26;
	t37 = (t49 - t53) * t24;
	t36 = -r_i_i_C(3) * t22 - pkin(4);
	t35 = -t18 * pkin(4) - pkin(2) - t47;
	t34 = t19 * r_i_i_C(3) * t39 + pkin(4) * t44 + t53 * t28;
	t33 = t36 * t19;
	t32 = t25 * r_i_i_C(1) - t21 * r_i_i_C(2);
	t29 = t36 * t18 - t48 + (t18 * t43 - t45) * r_i_i_C(2) + (-t18 * t40 - t46) * r_i_i_C(1);
	t12 = t18 * t38 - t42;
	t11 = -t18 * t39 - t41;
	t10 = -t18 * t41 - t39;
	t9 = t18 * t42 - t38;
	t2 = t10 * t25 - t24 * t46;
	t1 = -t10 * t21 - t24 * t45;
	t3 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t9 * r_i_i_C(3) + (t35 - t48) * t24, (-t49 + t50) * t28 + t34, -t28 * t49 + t34, t12 * r_i_i_C(3) + t32 * t11, (-t12 * t21 + t25 * t44) * r_i_i_C(1) + (-t12 * t25 - t21 * t44) * r_i_i_C(2), 0; t11 * r_i_i_C(3) - t32 * t12 + ((-t21 * r_i_i_C(1) - t25 * r_i_i_C(2) - pkin(5)) * t19 + t35) * t28, (t33 - t50) * t24 + t37, t24 * t33 + t37, t10 * r_i_i_C(3) + t32 * t9, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; 0, t29 - t47, t29, (r_i_i_C(3) * t26 - t32 * t22) * t19, -r_i_i_C(1) * t30 - r_i_i_C(2) * t31, 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (460->71), mult. (645->122), div. (0->0), fcn. (703->12), ass. (0->61)
	t45 = qJ(2) + qJ(3);
	t43 = sin(t45);
	t44 = cos(t45);
	t50 = cos(qJ(5));
	t46 = sin(qJ(5));
	t51 = cos(qJ(4));
	t77 = t46 * t51;
	t55 = t43 * t50 + t44 * t77;
	t88 = pkin(6) + r_i_i_C(3);
	t97 = t88 * t55;
	t42 = pkin(7) * qJ(5) - qJ(6);
	t37 = sin(t42);
	t38 = cos(t42);
	t64 = r_i_i_C(1) * t38 + r_i_i_C(2) * t37;
	t68 = t88 * t46;
	t96 = t64 * t50 + t68;
	t47 = sin(qJ(4));
	t53 = cos(qJ(1));
	t72 = t53 * t47;
	t49 = sin(qJ(1));
	t74 = t49 * t51;
	t30 = -t43 * t74 - t72;
	t80 = t44 * t49;
	t10 = t30 * t50 - t46 * t80;
	t71 = t53 * t51;
	t75 = t49 * t47;
	t29 = t43 * t75 - t71;
	t2 = t10 * t38 - t29 * t37;
	t62 = t10 * t37 + t29 * t38;
	t95 = r_i_i_C(1) * t62 - t2 * r_i_i_C(2);
	t32 = t43 * t71 - t75;
	t78 = t44 * t53;
	t12 = t32 * t50 + t46 * t78;
	t31 = t43 * t72 + t74;
	t94 = t12 * t37 - t31 * t38;
	t93 = t12 * t38 + t31 * t37;
	t73 = t50 * t51;
	t28 = -t43 * t46 + t44 * t73;
	t81 = t44 * t47;
	t91 = (-t28 * t37 + t38 * t81) * r_i_i_C(1) + (t28 * t38 + t37 * t81) * r_i_i_C(2);
	t90 = r_i_i_C(1) * t94 - r_i_i_C(2) * t93;
	t87 = cos(qJ(2)) * pkin(3);
	t86 = pkin(5) * t43;
	t85 = sin(qJ(2)) * pkin(3);
	t82 = t43 * t47;
	t79 = t44 * t50;
	t76 = t47 * t50;
	t70 = t44 * t75;
	t69 = t44 * t72;
	t20 = t28 * t49;
	t67 = (t20 * t37 - t38 * t70) * r_i_i_C(2) + (t20 * t38 + t37 * t70) * r_i_i_C(1) + (t86 + t97) * t49;
	t22 = t28 * t53;
	t66 = pkin(4) * t78 + (-t22 * t37 + t38 * t69) * r_i_i_C(2) + (-t22 * t38 - t37 * t69) * r_i_i_C(1) - t53 * t97;
	t65 = -pkin(4) * t43 - pkin(5) * t44;
	t63 = -r_i_i_C(1) * t37 + r_i_i_C(2) * t38;
	t59 = -pkin(2) - t85 + t65;
	t56 = -t32 * t46 + t50 * t78;
	t26 = -t43 * t73 - t44 * t46;
	t54 = t65 + (-t26 * t37 - t38 * t82) * r_i_i_C(2) + (-t26 * t38 + t37 * t82) * r_i_i_C(1) + t88 * (t43 * t77 - t79);
	t9 = -t30 * t46 - t49 * t79;
	t1 = [-t2 * r_i_i_C(1) - r_i_i_C(2) * t62 + t49 * t59 + t88 * t9, (-t86 + t87) * t53 + t66, -t53 * t86 + t66, t31 * t96 + t63 * t32, pkin(7) * t90 - t88 * t12 - t64 * t56, -t90; r_i_i_C(1) * t93 + r_i_i_C(2) * t94 + t59 * t53 - t88 * t56, (-pkin(4) * t44 - t87) * t49 + t67, -pkin(4) * t80 + t67, -t29 * t96 + t63 * t30, pkin(7) * t95 - t88 * t10 - t64 * t9, -t95; 0, t54 - t85, t54, ((-t37 * t51 + t38 * t76) * r_i_i_C(1) + (t37 * t76 + t38 * t51) * r_i_i_C(2) + t47 * t68) * t44, -pkin(7) * t91 - t88 * t28 + t64 * t55, t91;];
	Ja_transl = t1;
end