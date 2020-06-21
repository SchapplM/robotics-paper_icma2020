% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% CloosQRC350OL
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = CloosQRC350OL_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'CloosQRC350OL_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobiR_rot_sym_varpar: pkin has to be [6x1] (double)');
JR_rot=NaN(9,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t15 = t10 * t9;
	t12 = cos(qJ(1));
	t14 = t12 * t9;
	t11 = cos(qJ(2));
	t13 = t10 * t11;
	t8 = t12 * t11;
	t1 = [-t15, t8, 0, 0, 0, 0; t14, t13, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; -t13, -t14, 0, 0, 0, 0; t8, -t15, 0, 0, 0, 0; 0, -t11, 0, 0, 0, 0; -t12, 0, 0, 0, 0, 0; -t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:15
	% EndTime: 2020-06-20 08:27:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t24 = qJ(2) + qJ(3);
	t22 = sin(t24);
	t25 = sin(qJ(1));
	t28 = t25 * t22;
	t26 = cos(qJ(1));
	t27 = t26 * t22;
	t23 = cos(t24);
	t21 = t26 * t23;
	t20 = t25 * t23;
	t1 = [-t28, t21, t21, 0, 0, 0; t27, t20, t20, 0, 0, 0; 0, -t22, -t22, 0, 0, 0; -t20, -t27, -t27, 0, 0, 0; t21, -t28, -t28, 0, 0, 0; 0, -t23, -t23, 0, 0, 0; -t26, 0, 0, 0, 0, 0; -t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (53->22), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t91 = qJ(2) + qJ(3);
	t89 = sin(t91);
	t94 = cos(qJ(4));
	t104 = t89 * t94;
	t93 = sin(qJ(1));
	t103 = t93 * t89;
	t92 = sin(qJ(4));
	t102 = t93 * t92;
	t101 = t93 * t94;
	t95 = cos(qJ(1));
	t100 = t95 * t89;
	t99 = t95 * t92;
	t98 = t95 * t94;
	t90 = cos(t91);
	t97 = t90 * t102;
	t96 = t90 * t99;
	t88 = t89 * t92;
	t87 = t90 * t98;
	t86 = t90 * t101;
	t85 = t89 * t98 + t102;
	t84 = -t89 * t99 + t101;
	t83 = -t89 * t101 + t99;
	t82 = t89 * t102 + t98;
	t1 = [t83, t87, t87, t84, 0, 0; t85, t86, t86, -t82, 0, 0; 0, -t104, -t104, -t90 * t92, 0, 0; t82, -t96, -t96, -t85, 0, 0; t84, -t97, -t97, t83, 0, 0; 0, t88, t88, -t90 * t94, 0, 0; -t93 * t90, -t100, -t100, 0, 0, 0; t95 * t90, -t103, -t103, 0, 0, 0; 0, -t90, -t90, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:16
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (104->21), mult. (149->41), div. (0->0), fcn. (226->8), ass. (0->38)
	t133 = qJ(2) + qJ(3);
	t131 = sin(t133);
	t135 = sin(qJ(4));
	t153 = t131 * t135;
	t132 = cos(t133);
	t134 = sin(qJ(5));
	t152 = t132 * t134;
	t137 = cos(qJ(5));
	t151 = t132 * t137;
	t139 = cos(qJ(1));
	t150 = t132 * t139;
	t138 = cos(qJ(4));
	t149 = t134 * t138;
	t136 = sin(qJ(1));
	t148 = t136 * t135;
	t147 = t136 * t138;
	t146 = t137 * t138;
	t145 = t139 * t135;
	t144 = t139 * t138;
	t124 = t131 * t147 - t145;
	t143 = -t124 * t134 + t136 * t151;
	t142 = -t124 * t137 - t136 * t152;
	t141 = -t131 * t134 + t132 * t146;
	t140 = -t131 * t137 - t132 * t149;
	t128 = t132 * t145;
	t127 = t132 * t148;
	t126 = t131 * t144 + t148;
	t125 = t131 * t145 - t147;
	t123 = -t131 * t148 - t144;
	t122 = -t131 * t146 - t152;
	t121 = t131 * t149 - t151;
	t120 = t141 * t139;
	t119 = t140 * t139;
	t118 = t141 * t136;
	t117 = t140 * t136;
	t116 = t126 * t137 + t134 * t150;
	t115 = -t126 * t134 + t137 * t150;
	t1 = [t142, t120, t120, -t125 * t137, t115, 0; t116, t118, t118, t123 * t137, t143, 0; 0, t122, t122, -t135 * t151, t140, 0; -t143, t119, t119, t125 * t134, -t116, 0; t115, t117, t117, -t123 * t134, t142, 0; 0, t121, t121, t135 * t152, -t141, 0; t123, t128, t128, t126, 0, 0; t125, t127, t127, t124, 0, 0; 0, -t153, -t153, t132 * t138, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-06-20 08:27:16
	% EndTime: 2020-06-20 08:27:17
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (221->40), mult. (341->83), div. (0->0), fcn. (490->10), ass. (0->54)
	t199 = qJ(2) + qJ(3);
	t197 = sin(t199);
	t202 = sin(qJ(4));
	t207 = cos(qJ(1));
	t212 = t207 * t202;
	t203 = sin(qJ(1));
	t206 = cos(qJ(4));
	t215 = t203 * t206;
	t192 = t197 * t215 - t212;
	t205 = cos(qJ(5));
	t198 = cos(t199);
	t201 = sin(qJ(5));
	t223 = t198 * t201;
	t179 = t192 * t205 + t203 * t223;
	t211 = t207 * t206;
	t216 = t203 * t202;
	t191 = t197 * t216 + t211;
	t200 = sin(qJ(6));
	t204 = cos(qJ(6));
	t227 = t179 * t200 + t191 * t204;
	t226 = t179 * t204 - t191 * t200;
	t222 = t198 * t205;
	t221 = t198 * t207;
	t220 = t200 * t202;
	t219 = t200 * t205;
	t218 = t201 * t206;
	t217 = t202 * t204;
	t214 = t204 * t205;
	t213 = t205 * t206;
	t210 = t198 * t220;
	t209 = t198 * t217;
	t208 = t198 * t212;
	t178 = -t192 * t201 + t203 * t222;
	t190 = -t197 * t201 + t198 * t213;
	t189 = -t197 * t205 - t198 * t218;
	t194 = t197 * t211 + t216;
	t193 = t197 * t212 - t215;
	t188 = -t197 * t213 - t223;
	t187 = t197 * t218 - t222;
	t186 = t190 * t207;
	t185 = t189 * t207;
	t184 = t190 * t203;
	t183 = t189 * t203;
	t182 = t194 * t205 + t201 * t221;
	t181 = -t194 * t201 + t205 * t221;
	t177 = -t188 * t204 - t197 * t220;
	t176 = t188 * t200 - t197 * t217;
	t175 = -t186 * t204 + t200 * t208;
	t174 = t186 * t200 + t204 * t208;
	t173 = -t184 * t204 + t203 * t210;
	t172 = t184 * t200 + t203 * t209;
	t171 = t182 * t204 - t193 * t200;
	t170 = t182 * t200 + t193 * t204;
	t1 = [t226, t175, t175, t193 * t214 + t194 * t200, -t181 * t204, t170; -t171, t173, t173, t191 * t214 + t192 * t200, -t178 * t204, t227; 0, t177, t177, (t200 * t206 + t202 * t214) * t198, -t189 * t204, t190 * t200 + t209; -t227, t174, t174, -t193 * t219 + t194 * t204, t181 * t200, t171; t170, t172, t172, -t191 * t219 + t192 * t204, t178 * t200, t226; 0, t176, t176, (-t202 * t219 + t204 * t206) * t198, t189 * t200, t190 * t204 - t210; -t178, t185, t185, t193 * t201, -t182, 0; t181, t183, t183, t191 * t201, -t179, 0; 0, t187, t187, t202 * t223, -t190, 0;];
	JR_rot = t1;
end