% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% CloosQRC350OL
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = CloosQRC350OL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'CloosQRC350OL_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350OL_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
Ja_transl=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (9->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->8)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = pkin(2) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [-t4 * r_i_i_C(3) - t5 * t2, t7 * t4, 0, 0, 0, 0; -t2 * r_i_i_C(3) + t5 * t4, t7 * t2, 0, 0, 0, 0; 0, t6, 0, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (33->12), mult. (39->16), div. (0->0), fcn. (39->6), ass. (0->14)
	t6 = qJ(2) + qJ(3);
	t4 = sin(t6);
	t5 = cos(t6);
	t12 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t17 = -t12 - sin(qJ(2)) * pkin(3);
	t16 = r_i_i_C(1) * t5;
	t15 = r_i_i_C(2) * t4;
	t13 = cos(qJ(2)) * pkin(3) - t15;
	t11 = pkin(2) - t17;
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t2 = t10 * t16;
	t1 = t8 * t16;
	t3 = [-t10 * r_i_i_C(3) - t11 * t8, t13 * t10 + t2, -t10 * t15 + t2, 0, 0, 0; -t8 * r_i_i_C(3) + t11 * t10, t13 * t8 + t1, -t8 * t15 + t1, 0, 0, 0; 0, t17, -t12, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (96->26), mult. (119->37), div. (0->0), fcn. (127->8), ass. (0->28)
	t17 = cos(qJ(4));
	t39 = r_i_i_C(1) * t17 + pkin(4);
	t13 = qJ(2) + qJ(3);
	t11 = sin(t13);
	t12 = cos(t13);
	t36 = pkin(5) + r_i_i_C(3);
	t24 = t36 * t12;
	t33 = sin(qJ(2)) * pkin(3);
	t38 = t11 * pkin(4) + pkin(2) + t24 + t33;
	t37 = t39 * t12;
	t16 = sin(qJ(1));
	t35 = t37 * t16;
	t19 = cos(qJ(1));
	t34 = t37 * t19;
	t14 = sin(qJ(4));
	t30 = r_i_i_C(2) * t14;
	t29 = t14 * t19;
	t28 = t16 * t14;
	t27 = t16 * t17;
	t26 = t17 * t19;
	t22 = -t36 * t11 - t12 * t30;
	t21 = cos(qJ(2)) * pkin(3) + t22;
	t20 = -t24 + (t30 - t39) * t11;
	t4 = t11 * t26 + t28;
	t3 = -t11 * t29 + t27;
	t2 = -t11 * t27 + t29;
	t1 = t11 * t28 + t26;
	t5 = [r_i_i_C(1) * t2 + r_i_i_C(2) * t1 - t38 * t16, t21 * t19 + t34, t22 * t19 + t34, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t38 * t19, t21 * t16 + t35, t22 * t16 + t35, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t20 - t33, t20, (-r_i_i_C(1) * t14 - r_i_i_C(2) * t17) * t12, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (163->45), mult. (225->73), div. (0->0), fcn. (263->10), ass. (0->38)
	t22 = qJ(2) + qJ(3);
	t21 = cos(t22);
	t20 = sin(t22);
	t27 = cos(qJ(5));
	t23 = sin(qJ(5));
	t28 = cos(qJ(4));
	t47 = t23 * t28;
	t34 = -t20 * t27 - t21 * t47;
	t44 = t27 * t28;
	t35 = -t20 * t23 + t21 * t44;
	t24 = sin(qJ(4));
	t53 = r_i_i_C(3) * t24;
	t55 = r_i_i_C(1) * t35 + r_i_i_C(2) * t34 + t21 * t53;
	t54 = pkin(5) * t20;
	t52 = t21 * pkin(5);
	t51 = sin(qJ(2)) * pkin(3);
	t26 = sin(qJ(1));
	t50 = t21 * t26;
	t49 = t21 * t27;
	t30 = cos(qJ(1));
	t48 = t21 * t30;
	t46 = t26 * t24;
	t45 = t26 * t28;
	t43 = t30 * t24;
	t42 = t30 * t28;
	t40 = -t20 * pkin(4) - pkin(2) - t51;
	t39 = pkin(4) * t50 + t55 * t26;
	t38 = pkin(4) * t48 + t55 * t30;
	t37 = cos(qJ(2)) * pkin(3) - t54;
	t36 = t27 * r_i_i_C(1) - t23 * r_i_i_C(2);
	t31 = -t52 + (t20 * t47 - t49) * r_i_i_C(2) + (-t20 * t44 - t21 * t23) * r_i_i_C(1) + (-pkin(4) - t53) * t20;
	t12 = t20 * t42 + t46;
	t11 = t20 * t43 - t45;
	t10 = t20 * t45 - t43;
	t9 = -t20 * t46 - t42;
	t2 = t12 * t27 + t23 * t48;
	t1 = -t12 * t23 + t27 * t48;
	t3 = [t9 * r_i_i_C(3) - t36 * t10 + ((-t23 * r_i_i_C(1) - t27 * r_i_i_C(2) - pkin(5)) * t21 + t40) * t26, t37 * t30 + t38, -t30 * t54 + t38, t12 * r_i_i_C(3) - t36 * t11, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t11 * r_i_i_C(3) + (-t40 + t52) * t30, t37 * t26 + t39, -t26 * t54 + t39, t10 * r_i_i_C(3) + t36 * t9, (-t10 * t23 + t26 * t49) * r_i_i_C(1) + (-t10 * t27 - t23 * t50) * r_i_i_C(2), 0; 0, t31 - t51, t31, (r_i_i_C(3) * t28 - t36 * t24) * t21, t34 * r_i_i_C(1) - t35 * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 22:04:50
	% EndTime: 2020-06-23 22:04:50
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (342->64), mult. (503->117), div. (0->0), fcn. (621->12), ass. (0->58)
	t83 = pkin(6) + r_i_i_C(3);
	t41 = qJ(2) + qJ(3);
	t39 = sin(t41);
	t40 = cos(t41);
	t48 = cos(qJ(5));
	t43 = sin(qJ(5));
	t49 = cos(qJ(4));
	t74 = t43 * t49;
	t54 = -t39 * t48 - t40 * t74;
	t86 = t54 * t83;
	t44 = sin(qJ(4));
	t51 = cos(qJ(1));
	t68 = t51 * t44;
	t46 = sin(qJ(1));
	t70 = t46 * t49;
	t30 = t39 * t70 - t68;
	t78 = t40 * t46;
	t10 = t30 * t48 + t43 * t78;
	t67 = t51 * t49;
	t71 = t46 * t44;
	t29 = t39 * t71 + t67;
	t42 = sin(qJ(6));
	t47 = cos(qJ(6));
	t85 = t10 * t42 + t29 * t47;
	t84 = t10 * t47 - t29 * t42;
	t82 = pkin(5) * t39;
	t81 = sin(qJ(2)) * pkin(3);
	t77 = t40 * t48;
	t76 = t40 * t51;
	t75 = t42 * t44;
	t73 = t44 * t47;
	t72 = t44 * t48;
	t69 = t48 * t49;
	t66 = t40 * t75;
	t65 = t40 * t73;
	t64 = t40 * t68;
	t63 = t83 * t43;
	t28 = -t39 * t43 + t40 * t69;
	t20 = t28 * t46;
	t62 = (t20 * t42 + t46 * t65) * r_i_i_C(2) + pkin(4) * t78 + (-t20 * t47 + t46 * t66) * r_i_i_C(1) + t46 * t86;
	t22 = t28 * t51;
	t61 = pkin(4) * t76 + (t22 * t42 + t47 * t64) * r_i_i_C(2) + (-t22 * t47 + t42 * t64) * r_i_i_C(1) + t51 * t86;
	t60 = cos(qJ(2)) * pkin(3) - t82;
	t59 = t39 * pkin(4) + t40 * pkin(5);
	t58 = t47 * r_i_i_C(1) - t42 * r_i_i_C(2);
	t57 = t42 * r_i_i_C(1) + t47 * r_i_i_C(2);
	t56 = pkin(2) + t81 + t59;
	t55 = -t30 * t43 + t46 * t77;
	t26 = -t39 * t69 - t40 * t43;
	t53 = -t59 + (t26 * t42 - t39 * t73) * r_i_i_C(2) + (-t26 * t47 - t39 * t75) * r_i_i_C(1) + t83 * (t39 * t74 - t77);
	t52 = t58 * t48 + t63;
	t32 = t39 * t67 + t71;
	t31 = t39 * t68 - t70;
	t14 = t32 * t48 + t43 * t76;
	t13 = -t32 * t43 + t48 * t76;
	t2 = t14 * t47 - t31 * t42;
	t1 = t14 * t42 + t31 * t47;
	t3 = [t84 * r_i_i_C(1) - t85 * r_i_i_C(2) - t56 * t46 - t83 * t55, t60 * t51 + t61, -t51 * t82 + t61, t52 * t31 + t57 * t32, -t58 * t13 - t83 * t14, t1 * r_i_i_C(1) + t2 * r_i_i_C(2); -t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t83 * t13 + t56 * t51, t60 * t46 + t62, -t46 * t82 + t62, t52 * t29 + t57 * t30, -t83 * t10 - t58 * t55, t85 * r_i_i_C(1) + t84 * r_i_i_C(2); 0, t53 - t81, t53, ((t42 * t49 + t47 * t72) * r_i_i_C(1) + (-t42 * t72 + t47 * t49) * r_i_i_C(2) + t44 * t63) * t40, -t83 * t28 - t58 * t54, (t28 * t42 + t65) * r_i_i_C(1) + (t28 * t47 - t66) * r_i_i_C(2);];
	Ja_transl = t3;
end