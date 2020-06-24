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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (324->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t62 = qJ(2) + qJ(3);
	t60 = sin(t62);
	t61 = cos(t62);
	t64 = sin(qJ(1));
	t72 = t64 * t61;
	t55 = atan2(t72, -t60);
	t53 = sin(t55);
	t54 = cos(t55);
	t46 = t53 * t72 - t54 * t60;
	t45 = 0.1e1 / t46 ^ 2;
	t66 = cos(qJ(1));
	t78 = t45 * t66 ^ 2;
	t65 = cos(qJ(4));
	t68 = t66 * t65;
	t63 = sin(qJ(4));
	t71 = t64 * t63;
	t52 = t60 * t68 - t71;
	t50 = 0.1e1 / t52 ^ 2;
	t69 = t66 * t63;
	t70 = t64 * t65;
	t51 = t60 * t69 + t70;
	t77 = t50 * t51;
	t76 = t53 * t60;
	t59 = t61 ^ 2;
	t75 = 0.1e1 / t60 ^ 2 * t59;
	t74 = t61 * t66;
	t56 = 0.1e1 / (t64 ^ 2 * t75 + 0.1e1);
	t73 = t64 * t56;
	t67 = t51 ^ 2 * t50 + 0.1e1;
	t57 = 0.1e1 / t60;
	t49 = 0.1e1 / t52;
	t48 = 0.1e1 / t67;
	t47 = (0.1e1 + t75) * t73;
	t44 = 0.1e1 / t46;
	t43 = 0.1e1 / (t59 * t78 + 0.1e1);
	t42 = (t49 * t63 - t65 * t77) * t48 * t74;
	t41 = (-t60 * t44 - (-t64 * t76 - t54 * t61 + (t54 * t72 + t76) * t47) * t61 * t45) * t66 * t43;
	t1 = [-t57 * t56 * t74, t47, t47, 0, 0, 0; (-t44 * t72 - (-t54 * t57 * t59 * t73 + (-t56 + 0.1e1) * t61 * t53) * t61 * t78) * t43, t41, t41, 0, 0, 0; ((-t60 * t71 + t68) * t49 - (-t60 * t70 - t69) * t77) * t48, t42, t42, t67 * t48, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (600->31), mult. (873->84), div. (144->11), fcn. (1297->11), ass. (0->49)
	t82 = qJ(2) + qJ(3);
	t78 = sin(t82);
	t87 = cos(qJ(4));
	t88 = cos(qJ(1));
	t90 = t88 * t87;
	t84 = sin(qJ(4));
	t85 = sin(qJ(1));
	t93 = t85 * t84;
	t71 = t78 * t93 - t90;
	t79 = cos(t82);
	t97 = t79 * t84;
	t70 = atan2(t71, t97);
	t67 = sin(t70);
	t68 = cos(t70);
	t61 = t67 * t71 + t68 * t97;
	t60 = 0.1e1 / t61 ^ 2;
	t91 = t88 * t84;
	t92 = t85 * t87;
	t73 = t78 * t91 + t92;
	t102 = t60 * t73;
	t74 = t78 * t90 - t93;
	t83 = sin(qJ(5));
	t86 = cos(qJ(5));
	t94 = t79 * t88;
	t66 = t74 * t86 + t83 * t94;
	t64 = 0.1e1 / t66 ^ 2;
	t65 = t74 * t83 - t86 * t94;
	t101 = t64 * t65;
	t100 = t68 * t71;
	t99 = t73 ^ 2 * t60;
	t76 = 0.1e1 / t79;
	t80 = 0.1e1 / t84;
	t98 = t76 * t80;
	t96 = t79 * t85;
	t95 = t79 * t87;
	t89 = t65 ^ 2 * t64 + 0.1e1;
	t81 = 0.1e1 / t84 ^ 2;
	t77 = 0.1e1 / t79 ^ 2;
	t72 = t78 * t92 + t91;
	t69 = 0.1e1 / (t71 ^ 2 * t77 * t81 + 0.1e1);
	t63 = 0.1e1 / t66;
	t62 = 0.1e1 / t89;
	t59 = 0.1e1 / t61;
	t58 = (t71 * t77 * t78 * t80 + t85) * t69;
	t57 = 0.1e1 / (0.1e1 + t99);
	t56 = (-t71 * t81 * t87 + t72 * t80) * t76 * t69;
	t55 = ((t78 * t86 + t83 * t95) * t63 - (-t78 * t83 + t86 * t95) * t101) * t62 * t88;
	t54 = (-t58 * t100 * t102 + (t59 * t94 - (-t68 * t78 + (-t58 * t79 + t96) * t67) * t102) * t84) * t57;
	t1 = [t73 * t69 * t98, t58, t58, t56, 0, 0; (-t71 * t59 - (t67 + (t98 * t100 - t67) * t69) * t99) * t57, t54, t54, (t74 * t59 - (t68 * t95 + t67 * t72 + (-t67 * t97 + t100) * t56) * t102) * t57, 0, 0; ((-t72 * t83 + t86 * t96) * t63 - (-t72 * t86 - t83 * t96) * t101) * t62, t55, t55, (t86 * t101 - t83 * t63) * t73 * t62, t89 * t62, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:15:00
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (1654->49), mult. (2306->119), div. (145->9), fcn. (3132->13), ass. (0->62)
	t138 = qJ(2) + qJ(3);
	t136 = sin(t138);
	t140 = sin(qJ(4));
	t144 = cos(qJ(1));
	t153 = t144 * t140;
	t141 = sin(qJ(1));
	t143 = cos(qJ(4));
	t154 = t141 * t143;
	t129 = -t136 * t154 - t153;
	t139 = sin(qJ(5));
	t137 = cos(t138);
	t142 = cos(qJ(5));
	t158 = t137 * t142;
	t116 = t129 * t139 + t141 * t158;
	t156 = t139 * t143;
	t145 = t136 * t142 + t137 * t156;
	t115 = atan2(t116, -t145);
	t112 = sin(t115);
	t113 = cos(t115);
	t106 = t112 * t116 - t113 * t145;
	t105 = 0.1e1 / t106 ^ 2;
	t152 = t144 * t143;
	t155 = t141 * t140;
	t131 = t136 * t152 - t155;
	t157 = t137 * t144;
	t147 = -t131 * t139 + t142 * t157;
	t164 = t105 * t147;
	t119 = t131 * t142 + t139 * t157;
	t130 = t136 * t153 + t154;
	t135 = pkin(7) * qJ(5) - qJ(6);
	t132 = sin(t135);
	t133 = cos(t135);
	t148 = t119 * t133 + t130 * t132;
	t109 = 0.1e1 / t148 ^ 2;
	t110 = t119 * t132 - t130 * t133;
	t163 = t109 * t110;
	t162 = t113 * t116;
	t124 = 0.1e1 / t145 ^ 2;
	t161 = t116 * t124;
	t160 = t130 * t142;
	t159 = t137 * t140;
	t151 = t137 * t153;
	t150 = t110 ^ 2 * t109 + 0.1e1;
	t149 = t112 * t145 + t162;
	t146 = -t136 * t139 + t143 * t158;
	t128 = t136 * t155 - t152;
	t125 = t136 * t156 - t158;
	t123 = 0.1e1 / t145;
	t122 = t146 * t144;
	t121 = t145 * t141;
	t117 = -t141 * t137 * t139 + t129 * t142;
	t114 = 0.1e1 / (t116 ^ 2 * t124 + 0.1e1);
	t108 = 0.1e1 / t148;
	t107 = 0.1e1 / t150;
	t104 = 0.1e1 / t106;
	t103 = 0.1e1 / (t105 * t147 ^ 2 + 0.1e1);
	t102 = (-t123 * t128 - t159 * t161) * t139 * t114;
	t101 = (t121 * t123 - t125 * t161) * t114;
	t100 = (-t117 * t123 + t146 * t161) * t114;
	t99 = (-(t122 * t132 - t133 * t151) * t108 - (-t122 * t133 - t132 * t151) * t163) * t107;
	t98 = (-(t149 * t101 - t112 * t121 + t113 * t125) * t164 - t145 * t104 * t144) * t103;
	t1 = [-t147 * t123 * t114, t101, t101, t102, t100, 0; (-t116 * t104 - (t112 + (-t123 * t162 - t112) * t114) * t147 * t164) * t103, t98, t98, (t130 * t139 * t104 - ((t112 * t128 + t113 * t159) * t139 + t149 * t102) * t164) * t103, (-t119 * t104 - (t100 * t149 + t112 * t117 - t113 * t146) * t164) * t103, 0; (-(t117 * t132 + t128 * t133) * t108 - (-t117 * t133 + t128 * t132) * t163) * t107, t99, t99, (-(-t131 * t133 - t132 * t160) * t108 - (-t131 * t132 + t133 * t160) * t163) * t107, (-(pkin(7) * t148 + t132 * t147) * t108 - (pkin(7) * t110 - t133 * t147) * t163) * t107, t150 * t107;];
	Ja_rot = t1;
end