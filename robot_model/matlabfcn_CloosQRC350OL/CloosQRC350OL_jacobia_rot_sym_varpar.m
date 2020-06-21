% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% CloosQRC350OL
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
%   Wie in CloosQRC350OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = CloosQRC350OL_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350OL_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
Ja_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (358->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t61 = qJ(2) + qJ(3);
	t59 = sin(t61);
	t60 = cos(t61);
	t63 = sin(qJ(1));
	t71 = t63 * t60;
	t54 = atan2(-t71, -t59);
	t52 = sin(t54);
	t53 = cos(t54);
	t45 = -t52 * t71 - t53 * t59;
	t44 = 0.1e1 / t45 ^ 2;
	t65 = cos(qJ(1));
	t77 = t44 * t65 ^ 2;
	t64 = cos(qJ(4));
	t67 = t65 * t64;
	t62 = sin(qJ(4));
	t70 = t63 * t62;
	t51 = t59 * t67 + t70;
	t49 = 0.1e1 / t51 ^ 2;
	t68 = t65 * t62;
	t69 = t63 * t64;
	t50 = t59 * t68 - t69;
	t76 = t49 * t50;
	t75 = t52 * t59;
	t58 = t60 ^ 2;
	t73 = 0.1e1 / t59 ^ 2 * t58;
	t55 = 0.1e1 / (t63 ^ 2 * t73 + 0.1e1);
	t74 = t55 * t63;
	t72 = t60 * t65;
	t66 = t50 ^ 2 * t49 + 0.1e1;
	t56 = 0.1e1 / t59;
	t48 = 0.1e1 / t51;
	t47 = 0.1e1 / t66;
	t46 = (-0.1e1 - t73) * t74;
	t43 = 0.1e1 / t45;
	t42 = 0.1e1 / (t58 * t77 + 0.1e1);
	t41 = (t48 * t62 - t64 * t76) * t47 * t72;
	t40 = (-t59 * t43 - (t63 * t75 - t53 * t60 + (-t53 * t71 + t75) * t46) * t60 * t44) * t65 * t42;
	t1 = [t56 * t55 * t72, t46, t46, 0, 0, 0; (-t43 * t71 - (-t53 * t56 * t58 * t74 + (t55 - 0.1e1) * t60 * t52) * t60 * t77) * t42, t40, t40, 0, 0, 0; ((-t59 * t70 - t67) * t48 - (-t59 * t69 + t68) * t76) * t47, t41, t41, t66 * t47, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (600->33), mult. (873->84), div. (144->11), fcn. (1297->11), ass. (0->49)
	t86 = qJ(2) + qJ(3);
	t83 = cos(t86);
	t88 = sin(qJ(4));
	t101 = t83 * t88;
	t82 = sin(t86);
	t91 = cos(qJ(4));
	t92 = cos(qJ(1));
	t94 = t92 * t91;
	t89 = sin(qJ(1));
	t97 = t89 * t88;
	t73 = t82 * t97 + t94;
	t72 = atan2(-t73, t101);
	t69 = sin(t72);
	t70 = cos(t72);
	t63 = t70 * t101 - t69 * t73;
	t62 = 0.1e1 / t63 ^ 2;
	t95 = t92 * t88;
	t96 = t89 * t91;
	t76 = t82 * t95 - t96;
	t106 = t62 * t76;
	t77 = t82 * t94 + t97;
	t87 = sin(qJ(5));
	t90 = cos(qJ(5));
	t98 = t83 * t92;
	t68 = t77 * t90 + t87 * t98;
	t66 = 0.1e1 / t68 ^ 2;
	t67 = t77 * t87 - t90 * t98;
	t105 = t66 * t67;
	t104 = t70 * t73;
	t103 = t76 ^ 2 * t62;
	t80 = 0.1e1 / t83;
	t84 = 0.1e1 / t88;
	t102 = t80 * t84;
	t100 = t83 * t89;
	t99 = t83 * t91;
	t93 = t67 ^ 2 * t66 + 0.1e1;
	t85 = 0.1e1 / t88 ^ 2;
	t81 = 0.1e1 / t83 ^ 2;
	t75 = t82 * t96 - t95;
	t71 = 0.1e1 / (t73 ^ 2 * t81 * t85 + 0.1e1);
	t65 = 0.1e1 / t68;
	t64 = 0.1e1 / t93;
	t61 = 0.1e1 / t63;
	t60 = (-t73 * t81 * t82 * t84 - t89) * t71;
	t59 = 0.1e1 / (0.1e1 + t103);
	t58 = (t73 * t85 * t91 - t75 * t84) * t80 * t71;
	t57 = ((t82 * t90 + t87 * t99) * t65 - (-t82 * t87 + t90 * t99) * t105) * t64 * t92;
	t56 = (t60 * t104 * t106 + (t61 * t98 - (-t70 * t82 + (-t60 * t83 - t100) * t69) * t106) * t88) * t59;
	t1 = [-t76 * t71 * t102, t60, t60, t58, 0, 0; (-t73 * t61 - (-t69 + (t102 * t104 + t69) * t71) * t103) * t59, t56, t56, (t77 * t61 - (t70 * t99 - t69 * t75 + (-t69 * t101 - t104) * t58) * t106) * t59, 0, 0; ((t90 * t100 - t75 * t87) * t65 - (-t87 * t100 - t75 * t90) * t105) * t64, t57, t57, (t90 * t105 - t87 * t65) * t76 * t64, t93 * t64, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (1532->48), mult. (2174->116), div. (145->9), fcn. (3100->13), ass. (0->60)
	t129 = qJ(2) + qJ(3);
	t127 = sin(t129);
	t136 = cos(qJ(4));
	t137 = cos(qJ(1));
	t142 = t137 * t136;
	t132 = sin(qJ(4));
	t133 = sin(qJ(1));
	t145 = t133 * t132;
	t122 = t127 * t142 + t145;
	t131 = sin(qJ(5));
	t135 = cos(qJ(5));
	t128 = cos(t129);
	t147 = t128 * t137;
	t111 = t122 * t135 + t131 * t147;
	t143 = t137 * t132;
	t144 = t133 * t136;
	t121 = t127 * t143 - t144;
	t130 = sin(qJ(6));
	t134 = cos(qJ(6));
	t101 = -t111 * t134 + t121 * t130;
	t98 = 0.1e1 / t101 ^ 2;
	t99 = t111 * t130 + t121 * t134;
	t154 = t98 * t99;
	t109 = t122 * t131 - t135 * t147;
	t120 = t127 * t144 - t143;
	t148 = t128 * t135;
	t106 = t120 * t131 - t133 * t148;
	t146 = t131 * t136;
	t138 = t127 * t135 + t128 * t146;
	t105 = atan2(t106, -t138);
	t102 = sin(t105);
	t103 = cos(t105);
	t95 = t102 * t106 - t103 * t138;
	t94 = 0.1e1 / t95 ^ 2;
	t153 = t109 * t94;
	t152 = t109 ^ 2 * t94;
	t115 = 0.1e1 / t138 ^ 2;
	t151 = t106 * t115;
	t150 = t121 * t135;
	t149 = t128 * t132;
	t141 = t99 ^ 2 * t98 + 0.1e1;
	t140 = t128 * t143;
	t107 = t133 * t128 * t131 + t120 * t135;
	t139 = -t127 * t131 + t136 * t148;
	t119 = -t127 * t145 - t142;
	t116 = t127 * t146 - t148;
	t114 = 0.1e1 / t138;
	t113 = t139 * t137;
	t112 = t138 * t133;
	t104 = 0.1e1 / (t106 ^ 2 * t115 + 0.1e1);
	t97 = 0.1e1 / t101;
	t96 = 0.1e1 / t141;
	t93 = 0.1e1 / t95;
	t92 = 0.1e1 / (0.1e1 + t152);
	t91 = (-t114 * t119 - t149 * t151) * t131 * t104;
	t90 = (-t112 * t114 - t116 * t151) * t104;
	t89 = (-t107 * t114 + t139 * t151) * t104;
	t88 = ((-t113 * t130 - t134 * t140) * t97 + (-t113 * t134 + t130 * t140) * t154) * t96;
	t87 = (((t106 * t90 + t116) * t103 + (t138 * t90 + t112) * t102) * t153 - t138 * t93 * t137) * t92;
	t1 = [-t109 * t114 * t104, t90, t90, t91, t89, 0; (t106 * t93 + (t102 + (-t103 * t106 * t114 - t102) * t104) * t152) * t92, t87, t87, (t121 * t131 * t93 + ((t106 * t91 + t131 * t149) * t103 + (t119 * t131 + t138 * t91) * t102) * t153) * t92, (-t111 * t93 + ((t106 * t89 - t139) * t103 + (t138 * t89 + t107) * t102) * t153) * t92, 0; ((t107 * t130 - t119 * t134) * t97 + (t107 * t134 + t119 * t130) * t154) * t96, t88, t88, ((-t122 * t134 + t130 * t150) * t97 + (t122 * t130 + t134 * t150) * t154) * t96, (t130 * t97 + t134 * t154) * t96 * t109, t141 * t96;];
	Ja_rot = t1;
end