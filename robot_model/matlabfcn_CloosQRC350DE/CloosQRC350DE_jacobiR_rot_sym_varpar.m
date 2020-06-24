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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t1 = [-t9, 0, 0, 0, 0, 0; -t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(2));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(2));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t15, t12, 0, 0, 0, 0; -t14, -t13, 0, 0, 0, 0; 0, -t8, 0, 0, 0, 0; -t13, -t14, 0, 0, 0, 0; -t12, t15, 0, 0, 0, 0; 0, -t10, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t24 = qJ(2) + qJ(3);
	t23 = cos(t24);
	t25 = sin(qJ(1));
	t28 = t25 * t23;
	t22 = sin(t24);
	t26 = cos(qJ(1));
	t27 = t26 * t22;
	t21 = t26 * t23;
	t20 = t25 * t22;
	t1 = [-t20, t21, t21, 0, 0, 0; -t27, -t28, -t28, 0, 0, 0; 0, -t22, -t22, 0, 0, 0; -t28, -t27, -t27, 0, 0, 0; -t21, t20, t20, 0, 0, 0; 0, -t23, -t23, 0, 0, 0; t26, 0, 0, 0, 0, 0; -t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:14:59
	% EndTime: 2020-06-23 21:14:59
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (52->21), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t93 = qJ(2) + qJ(3);
	t91 = sin(t93);
	t96 = cos(qJ(4));
	t105 = t91 * t96;
	t94 = sin(qJ(4));
	t95 = sin(qJ(1));
	t104 = t95 * t94;
	t103 = t95 * t96;
	t97 = cos(qJ(1));
	t102 = t97 * t91;
	t101 = t97 * t94;
	t100 = t97 * t96;
	t92 = cos(t93);
	t99 = t92 * t103;
	t98 = t92 * t101;
	t90 = t95 * t91;
	t89 = t91 * t94;
	t88 = t92 * t100;
	t87 = t92 * t104;
	t86 = -t91 * t100 + t104;
	t85 = t91 * t101 + t103;
	t84 = t91 * t103 + t101;
	t83 = t91 * t104 - t100;
	t1 = [-t84, t88, t88, -t85, 0, 0; t86, -t99, -t99, t83, 0, 0; 0, -t105, -t105, -t92 * t94, 0, 0; t83, -t98, -t98, t86, 0, 0; t85, t87, t87, t84, 0, 0; 0, t89, t89, -t92 * t96, 0, 0; -t95 * t92, -t102, -t102, 0, 0, 0; -t97 * t92, t90, t90, 0, 0, 0; 0, -t92, -t92, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:00
	% EndTime: 2020-06-23 21:15:00
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (106->29), mult. (149->41), div. (0->0), fcn. (226->8), ass. (0->38)
	t127 = qJ(2) + qJ(3);
	t125 = sin(t127);
	t129 = sin(qJ(4));
	t148 = t125 * t129;
	t126 = cos(t127);
	t128 = sin(qJ(5));
	t147 = t126 * t128;
	t131 = cos(qJ(5));
	t146 = t126 * t131;
	t133 = cos(qJ(1));
	t145 = t126 * t133;
	t132 = cos(qJ(4));
	t144 = t128 * t132;
	t130 = sin(qJ(1));
	t143 = t130 * t129;
	t142 = t130 * t132;
	t141 = t131 * t132;
	t140 = t133 * t129;
	t139 = t133 * t132;
	t138 = t126 * t143;
	t122 = t125 * t139 - t143;
	t137 = -t122 * t128 + t131 * t145;
	t136 = -t122 * t131 - t128 * t145;
	t135 = -t125 * t128 + t126 * t141;
	t134 = t125 * t131 + t126 * t144;
	t123 = t126 * t140;
	t121 = -t125 * t140 - t142;
	t120 = -t125 * t142 - t140;
	t119 = t125 * t143 - t139;
	t118 = -t125 * t141 - t147;
	t117 = t125 * t144 - t146;
	t116 = t135 * t133;
	t115 = t134 * t133;
	t114 = t135 * t130;
	t113 = t134 * t130;
	t112 = t120 * t131 - t130 * t147;
	t111 = -t120 * t128 - t130 * t146;
	t1 = [t112, t116, t116, t121 * t131, t137, 0; t136, -t114, -t114, t119 * t131, t111, 0; 0, t118, t118, -t129 * t146, -t134, 0; t111, -t115, -t115, -t121 * t128, t136, 0; -t137, t113, t113, -t119 * t128, -t112, 0; 0, t117, t117, t129 * t147, -t135, 0; -t119, t123, t123, t122, 0, 0; t121, -t138, -t138, t120, 0, 0; 0, -t148, -t148, t126 * t132, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-23 21:15:00
	% EndTime: 2020-06-23 21:15:01
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (336->49), mult. (483->92), div. (0->0), fcn. (572->10), ass. (0->57)
	t215 = qJ(2) + qJ(3);
	t213 = sin(t215);
	t220 = cos(qJ(4));
	t221 = cos(qJ(1));
	t230 = t221 * t220;
	t217 = sin(qJ(4));
	t218 = sin(qJ(1));
	t234 = t218 * t217;
	t207 = t213 * t230 - t234;
	t216 = sin(qJ(5));
	t219 = cos(qJ(5));
	t214 = cos(t215);
	t237 = t214 * t221;
	t194 = t207 * t219 + t216 * t237;
	t231 = t221 * t217;
	t233 = t218 * t220;
	t206 = t213 * t231 + t233;
	t212 = pkin(7) * qJ(5) - qJ(6);
	t208 = sin(t212);
	t209 = cos(t212);
	t225 = t194 * t208 - t206 * t209;
	t226 = t194 * t209 + t206 * t208;
	t243 = t208 * t219;
	t242 = t209 * t219;
	t241 = t213 * t217;
	t240 = t214 * t216;
	t239 = t214 * t217;
	t238 = t214 * t219;
	t236 = t216 * t220;
	t235 = t217 * t219;
	t232 = t219 * t220;
	t229 = t214 * t234;
	t228 = t214 * t231;
	t205 = -t213 * t233 - t231;
	t192 = t205 * t219 - t218 * t240;
	t204 = t213 * t234 - t230;
	t184 = t192 * t209 - t204 * t208;
	t227 = t192 * t208 + t204 * t209;
	t203 = -t213 * t216 + t214 * t232;
	t224 = -t203 * t208 + t209 * t239;
	t223 = t203 * t209 + t208 * t239;
	t193 = -t207 * t216 + t219 * t237;
	t222 = t213 * t219 + t214 * t236;
	t201 = -t213 * t232 - t240;
	t200 = t213 * t236 - t238;
	t199 = t203 * t221;
	t198 = t222 * t221;
	t197 = t203 * t218;
	t196 = t222 * t218;
	t191 = -t205 * t216 - t218 * t238;
	t190 = -t201 * t209 + t208 * t241;
	t189 = -t201 * t208 - t209 * t241;
	t188 = -t199 * t209 - t208 * t228;
	t187 = -t199 * t208 + t209 * t228;
	t186 = t197 * t209 + t208 * t229;
	t185 = t197 * t208 - t209 * t229;
	t1 = [-t184, t188, t188, t206 * t242 - t207 * t208, t225 * pkin(7) - t193 * t209, -t225; t226, t186, t186, -t204 * t242 - t205 * t208, t227 * pkin(7) - t191 * t209, -t227; 0, t190, t190, (-t208 * t220 + t209 * t235) * t214, -t224 * pkin(7) + t209 * t222, t224; -t227, t187, t187, t206 * t243 + t207 * t209, -t226 * pkin(7) - t193 * t208, t226; t225, t185, t185, -t204 * t243 + t205 * t209, -t184 * pkin(7) - t191 * t208, t184; 0, t189, t189, (t208 * t235 + t209 * t220) * t214, -t223 * pkin(7) + t208 * t222, t223; t191, -t198, -t198, t206 * t216, -t194, 0; -t193, t196, t196, -t204 * t216, -t192, 0; 0, t200, t200, t216 * t239, -t203, 0;];
	JR_rot = t1;
end