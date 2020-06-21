% Calculate kinetic energy for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = CloosQRC350OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 07:59:58
% EndTime: 2020-06-20 08:00:06
% DurationCPUTime: 7.98s
% Computational Cost: add. (2573->368), mult. (3404->562), div. (0->0), fcn. (3784->12), ass. (0->172)
t190 = sin(qJ(1));
t238 = pkin(2) * t190;
t195 = cos(qJ(1));
t237 = pkin(2) * t195;
t194 = cos(qJ(2));
t236 = pkin(3) * t194;
t189 = sin(qJ(2));
t235 = pkin(3) * t189;
t234 = Icges(2,4) * t190;
t233 = Icges(3,4) * t189;
t232 = Icges(3,4) * t194;
t185 = qJ(2) + qJ(3);
t181 = sin(t185);
t231 = Icges(4,4) * t181;
t182 = cos(t185);
t230 = Icges(4,4) * t182;
t188 = sin(qJ(4));
t229 = t182 * t188;
t228 = t182 * t190;
t193 = cos(qJ(4));
t227 = t182 * t193;
t226 = t182 * t195;
t225 = t188 * t195;
t224 = t190 * t188;
t223 = t190 * t193;
t222 = t193 * t195;
t221 = qJD(4) * t182;
t220 = V_base(5) * pkin(1) + V_base(1);
t173 = qJD(2) * t195 + V_base(5);
t178 = V_base(6) + qJD(1);
t143 = t235 * t190;
t217 = -t143 - t238;
t216 = t173 * t236 + t220;
t150 = qJD(3) * t195 + t173;
t215 = pkin(4) * t181 + pkin(5) * t182;
t214 = rSges(3,1) * t189 + rSges(3,2) * t194;
t213 = rSges(4,1) * t181 + rSges(4,2) * t182;
t121 = t190 * t221 + t150;
t212 = Icges(3,1) * t189 + t232;
t211 = Icges(4,1) * t181 + t230;
t210 = Icges(3,2) * t194 + t233;
t209 = Icges(4,2) * t182 + t231;
t208 = Icges(3,5) * t189 + Icges(3,6) * t194;
t207 = Icges(4,5) * t181 + Icges(4,6) * t182;
t206 = -V_base(4) * pkin(1) + t178 * t237 + V_base(2);
t156 = -qJD(4) * t181 + t178;
t139 = t181 * t224 + t222;
t94 = qJD(5) * t139 + t121;
t151 = V_base(4) + (-qJD(2) - qJD(3)) * t190;
t205 = -t237 * V_base(5) + V_base(4) * t238 + V_base(3);
t123 = qJD(5) * t229 + t156;
t122 = t195 * t221 + t151;
t204 = (Icges(4,3) * t195 + t190 * t207) * t150 + (-Icges(4,3) * t190 + t195 * t207) * t151 + (Icges(4,5) * t182 - Icges(4,6) * t181) * t178;
t174 = -qJD(2) * t190 + V_base(4);
t203 = (Icges(3,3) * t195 + t190 * t208) * t173 + (-Icges(3,3) * t190 + t195 * t208) * t174 + (Icges(3,5) * t194 - Icges(3,6) * t189) * t178;
t141 = t181 * t225 - t223;
t95 = qJD(5) * t141 + t122;
t144 = t235 * t195;
t202 = t178 * t144 - t174 * t236 + t206;
t201 = t174 * t143 - t173 * t144 + t205;
t137 = t215 * t190;
t149 = pkin(4) * t182 - pkin(5) * t181;
t200 = t150 * t149 + (-t137 + t217) * t178 + t216;
t138 = t215 * t195;
t199 = t178 * t138 - t149 * t151 + t202;
t198 = t151 * t137 - t150 * t138 + t201;
t114 = Icges(4,6) * t195 + t190 * t209;
t115 = -Icges(4,6) * t190 + t195 * t209;
t116 = Icges(4,5) * t195 + t190 * t211;
t117 = -Icges(4,5) * t190 + t195 * t211;
t146 = -Icges(4,2) * t181 + t230;
t147 = Icges(4,1) * t182 - t231;
t197 = (t115 * t182 + t117 * t181) * t151 + (t114 * t182 + t116 * t181) * t150 + (t146 * t182 + t147 * t181) * t178;
t126 = Icges(3,6) * t195 + t190 * t210;
t127 = -Icges(3,6) * t190 + t195 * t210;
t128 = Icges(3,5) * t195 + t190 * t212;
t129 = -Icges(3,5) * t190 + t195 * t212;
t160 = -Icges(3,2) * t189 + t232;
t163 = Icges(3,1) * t194 - t233;
t196 = (t127 * t194 + t129 * t189) * t174 + (t126 * t194 + t128 * t189) * t173 + (t160 * t194 + t163 * t189) * t178;
t192 = cos(qJ(5));
t191 = cos(qJ(6));
t187 = sin(qJ(5));
t186 = sin(qJ(6));
t183 = Icges(2,4) * t195;
t168 = rSges(2,1) * t195 - t190 * rSges(2,2);
t167 = rSges(3,1) * t194 - rSges(3,2) * t189;
t166 = t190 * rSges(2,1) + rSges(2,2) * t195;
t165 = Icges(2,1) * t195 - t234;
t164 = Icges(2,1) * t190 + t183;
t162 = -Icges(2,2) * t190 + t183;
t161 = Icges(2,2) * t195 + t234;
t155 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t154 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t153 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t148 = rSges(4,1) * t182 - rSges(4,2) * t181;
t142 = t181 * t222 + t224;
t140 = t181 * t223 - t225;
t136 = -t181 * t187 + t192 * t227;
t135 = -t181 * t192 - t187 * t227;
t134 = -t190 * rSges(3,3) + t195 * t214;
t133 = rSges(3,3) * t195 + t190 * t214;
t119 = -t190 * rSges(4,3) + t195 * t213;
t118 = rSges(4,3) * t195 + t190 * t213;
t110 = -rSges(5,3) * t181 + (rSges(5,1) * t193 - rSges(5,2) * t188) * t182;
t109 = -Icges(5,5) * t181 + (Icges(5,1) * t193 - Icges(5,4) * t188) * t182;
t108 = -Icges(5,6) * t181 + (Icges(5,4) * t193 - Icges(5,2) * t188) * t182;
t107 = -Icges(5,3) * t181 + (Icges(5,5) * t193 - Icges(5,6) * t188) * t182;
t106 = V_base(5) * rSges(2,3) - t166 * t178 + t220;
t105 = t168 * t178 + V_base(2) + (-pkin(1) - rSges(2,3)) * V_base(4);
t103 = t166 * V_base(4) - t168 * V_base(5) + V_base(3);
t101 = t142 * t192 + t187 * t226;
t100 = -t142 * t187 + t192 * t226;
t99 = t140 * t192 + t187 * t228;
t98 = -t140 * t187 + t192 * t228;
t97 = -t136 * t191 + t186 * t229;
t96 = t136 * t186 + t191 * t229;
t93 = qJD(6) * t135 + t123;
t92 = t142 * rSges(5,1) - t141 * rSges(5,2) + rSges(5,3) * t226;
t91 = rSges(5,1) * t140 - rSges(5,2) * t139 + rSges(5,3) * t228;
t90 = Icges(5,1) * t142 - Icges(5,4) * t141 + Icges(5,5) * t226;
t89 = Icges(5,1) * t140 - Icges(5,4) * t139 + Icges(5,5) * t228;
t88 = Icges(5,4) * t142 - Icges(5,2) * t141 + Icges(5,6) * t226;
t87 = Icges(5,4) * t140 - Icges(5,2) * t139 + Icges(5,6) * t228;
t86 = Icges(5,5) * t142 - Icges(5,6) * t141 + Icges(5,3) * t226;
t85 = Icges(5,5) * t140 - Icges(5,6) * t139 + Icges(5,3) * t228;
t84 = rSges(6,1) * t136 + rSges(6,2) * t135 + rSges(6,3) * t229;
t83 = Icges(6,1) * t136 + Icges(6,4) * t135 + Icges(6,5) * t229;
t82 = Icges(6,4) * t136 + Icges(6,2) * t135 + Icges(6,6) * t229;
t81 = Icges(6,5) * t136 + Icges(6,6) * t135 + Icges(6,3) * t229;
t80 = -t101 * t191 + t141 * t186;
t79 = t101 * t186 + t141 * t191;
t78 = t139 * t186 - t191 * t99;
t77 = t139 * t191 + t186 * t99;
t76 = t167 * t173 + (-t133 - t238) * t178 + t220;
t75 = t134 * t178 - t167 * t174 + t206;
t74 = t174 * t133 - t173 * t134 + t205;
t73 = qJD(6) * t100 + t95;
t72 = qJD(6) * t98 + t94;
t71 = rSges(6,1) * t101 + rSges(6,2) * t100 + rSges(6,3) * t141;
t70 = rSges(6,1) * t99 + rSges(6,2) * t98 + rSges(6,3) * t139;
t69 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t141;
t68 = Icges(6,1) * t99 + Icges(6,4) * t98 + Icges(6,5) * t139;
t67 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t141;
t66 = Icges(6,4) * t99 + Icges(6,2) * t98 + Icges(6,6) * t139;
t65 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t141;
t64 = Icges(6,5) * t99 + Icges(6,6) * t98 + Icges(6,3) * t139;
t63 = rSges(7,1) * t97 + rSges(7,2) * t96 + rSges(7,3) * t135;
t62 = Icges(7,1) * t97 + Icges(7,4) * t96 + Icges(7,5) * t135;
t61 = Icges(7,4) * t97 + Icges(7,2) * t96 + Icges(7,6) * t135;
t60 = Icges(7,5) * t97 + Icges(7,6) * t96 + Icges(7,3) * t135;
t59 = t148 * t150 + (-t118 + t217) * t178 + t216;
t58 = t119 * t178 - t148 * t151 + t202;
t57 = t151 * t118 - t150 * t119 + t201;
t56 = rSges(7,1) * t80 + rSges(7,2) * t79 + rSges(7,3) * t100;
t55 = rSges(7,1) * t78 + rSges(7,2) * t77 + rSges(7,3) * t98;
t54 = Icges(7,1) * t80 + Icges(7,4) * t79 + Icges(7,5) * t100;
t53 = Icges(7,1) * t78 + Icges(7,4) * t77 + Icges(7,5) * t98;
t52 = Icges(7,4) * t80 + Icges(7,2) * t79 + Icges(7,6) * t100;
t51 = Icges(7,4) * t78 + Icges(7,2) * t77 + Icges(7,6) * t98;
t50 = Icges(7,5) * t80 + Icges(7,6) * t79 + Icges(7,3) * t100;
t49 = Icges(7,5) * t78 + Icges(7,6) * t77 + Icges(7,3) * t98;
t48 = t110 * t121 - t156 * t91 + t200;
t47 = -t110 * t122 + t156 * t92 + t199;
t46 = -t121 * t92 + t122 * t91 + t198;
t45 = -t123 * t70 + t84 * t94 + t200;
t44 = t123 * t71 - t84 * t95 + t199;
t43 = t95 * t70 - t94 * t71 + t198;
t42 = -t55 * t93 + t63 * t72 + (-t98 * t123 + t135 * t94) * pkin(6) + t200;
t41 = t56 * t93 - t63 * t73 + (t100 * t123 - t135 * t95) * pkin(6) + t199;
t40 = t73 * t55 - t72 * t56 + (-t100 * t94 + t98 * t95) * pkin(6) + t198;
t1 = ((-t190 * t161 + t164 * t195 + Icges(1,4)) * V_base(5) + (-t190 * t162 + t195 * t165 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t195 * t161 + t190 * t164 + Icges(1,2)) * V_base(5) + (t162 * t195 + t190 * t165 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + m(3) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(6) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + m(5) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + m(7) * (t40 ^ 2 + t41 ^ 2 + t42 ^ 2) / 0.2e1 + V_base(5) * t178 * (Icges(2,5) * t190 + Icges(2,6) * t195) + t178 * V_base(4) * (Icges(2,5) * t195 - Icges(2,6) * t190) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((-t127 * t189 + t129 * t194) * t174 + (-t126 * t189 + t128 * t194) * t173 + (-t115 * t181 + t117 * t182) * t151 + (-t114 * t181 + t116 * t182) * t150 + (-t181 * t146 + t182 * t147 - t189 * t160 + t194 * t163 + Icges(2,3)) * t178) * t178 / 0.2e1 + t156 * ((-t107 * t156 - t85 * t121 - t86 * t122) * t181 + ((-t188 * t88 + t193 * t90) * t122 + (-t188 * t87 + t193 * t89) * t121 + (-t108 * t188 + t109 * t193) * t156) * t182) / 0.2e1 + t173 * (t196 * t190 + t203 * t195) / 0.2e1 + t174 * (-t203 * t190 + t196 * t195) / 0.2e1 + t150 * (t197 * t190 + t204 * t195) / 0.2e1 + t151 * (-t190 * t204 + t195 * t197) / 0.2e1 + m(1) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + t95 * ((t100 * t67 + t101 * t69 + t141 * t65) * t95 + (t100 * t66 + t101 * t68 + t141 * t64) * t94 + (t100 * t82 + t101 * t83 + t141 * t81) * t123) / 0.2e1 + t94 * ((t139 * t65 + t67 * t98 + t69 * t99) * t95 + (t139 * t64 + t98 * t66 + t99 * t68) * t94 + (t139 * t81 + t82 * t98 + t83 * t99) * t123) / 0.2e1 + t93 * ((t135 * t50 + t52 * t96 + t54 * t97) * t73 + (t135 * t49 + t51 * t96 + t53 * t97) * t72 + (t135 * t60 + t96 * t61 + t97 * t62) * t93) / 0.2e1 + m(2) * (t103 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + t72 * ((t50 * t98 + t52 * t77 + t54 * t78) * t73 + (t49 * t98 + t51 * t77 + t53 * t78) * t72 + (t60 * t98 + t61 * t77 + t62 * t78) * t93) / 0.2e1 + t73 * ((t100 * t50 + t52 * t79 + t54 * t80) * t73 + (t100 * t49 + t51 * t79 + t53 * t80) * t72 + (t100 * t60 + t61 * t79 + t62 * t80) * t93) / 0.2e1 + t122 * ((-t141 * t88 + t142 * t90 + t86 * t226) * t122 + (-t141 * t87 + t142 * t89 + t226 * t85) * t121 + (t107 * t226 - t141 * t108 + t142 * t109) * t156) / 0.2e1 + t121 * ((-t139 * t88 + t140 * t90 + t228 * t86) * t122 + (-t139 * t87 + t140 * t89 + t85 * t228) * t121 + (t107 * t228 - t108 * t139 + t140 * t109) * t156) / 0.2e1 + t123 * ((t135 * t67 + t136 * t69 + t229 * t65) * t95 + (t135 * t66 + t136 * t68 + t229 * t64) * t94 + (t135 * t82 + t136 * t83 + t81 * t229) * t123) / 0.2e1;
T = t1;
