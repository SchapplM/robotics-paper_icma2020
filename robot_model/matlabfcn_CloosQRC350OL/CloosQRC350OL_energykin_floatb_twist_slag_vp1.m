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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 21:55:41
% EndTime: 2020-06-23 21:55:49
% DurationCPUTime: 7.96s
% Computational Cost: add. (1571->249), mult. (2082->445), div. (0->0), fcn. (2092->12), ass. (0->157)
t144 = sin(qJ(5));
t149 = cos(qJ(5));
t142 = qJ(2) + qJ(3);
t135 = cos(t142);
t147 = sin(qJ(1));
t185 = t135 * t147;
t134 = sin(t142);
t150 = cos(qJ(4));
t178 = t147 * t150;
t145 = sin(qJ(4));
t152 = cos(qJ(1));
t182 = t145 * t152;
t91 = t134 * t178 - t182;
t59 = -t144 * t91 + t149 * t185;
t213 = t59 ^ 2 / 0.2e1;
t176 = t150 * t152;
t179 = t147 * t145;
t89 = t134 * t179 + t176;
t212 = t89 ^ 2;
t92 = t134 * t182 - t178;
t211 = t92 ^ 2;
t209 = pkin(2) * t147;
t208 = pkin(2) * t152;
t183 = t135 * t152;
t94 = t134 * t176 + t179;
t61 = -t94 * t144 + t149 * t183;
t207 = t59 * t61;
t206 = t89 * t92;
t146 = sin(qJ(2));
t205 = pkin(3) * t146;
t184 = t135 * t150;
t86 = -t134 * t144 + t149 * t184;
t204 = Icges(6,1) * t86;
t143 = sin(qJ(6));
t148 = cos(qJ(6));
t186 = t135 * t145;
t56 = t143 * t186 - t148 * t86;
t203 = Icges(7,1) * t56;
t85 = -t134 * t149 - t144 * t184;
t202 = Icges(6,2) * t85;
t55 = t143 * t86 + t148 * t186;
t201 = Icges(7,2) * t55;
t200 = Icges(7,3) * t85;
t199 = Icges(5,1) * t150;
t198 = Icges(2,4) * t147;
t197 = Icges(4,4) * t134;
t196 = Icges(4,4) * t135;
t195 = Icges(3,5) * t147;
t194 = Icges(3,5) * t152;
t151 = cos(qJ(2));
t140 = t151 ^ 2;
t193 = Icges(3,2) * t140;
t192 = Icges(5,2) * t145;
t130 = t135 ^ 2;
t191 = Icges(5,3) * t130;
t190 = Icges(5,3) * t134;
t125 = qJD(2) * t152 + V_base(5);
t189 = t125 * t151;
t126 = -qJD(2) * t147 + V_base(4);
t188 = t126 * t151;
t131 = V_base(6) + qJD(1);
t187 = t131 * t151;
t181 = t146 * t147;
t180 = t146 * t152;
t177 = t147 * t152;
t175 = qJD(4) * t135;
t174 = V_base(5) * pkin(1) + V_base(1);
t173 = Icges(6,3) * t186;
t95 = t205 * t147;
t172 = -t95 - t209;
t171 = t146 * (Icges(3,1) - Icges(3,2));
t170 = pkin(3) * t189 + t174;
t102 = qJD(3) * t152 + t125;
t169 = pkin(4) * t134 + pkin(5) * t135;
t78 = t147 * t175 + t102;
t53 = qJD(5) * t89 + t78;
t103 = V_base(4) + (-qJD(2) - qJD(3)) * t147;
t79 = t152 * t175 + t103;
t54 = qJD(5) * t92 + t79;
t168 = -t53 * t61 + t54 * t59;
t111 = -qJD(4) * t134 + t131;
t80 = qJD(5) * t186 + t111;
t167 = t53 * t85 - t59 * t80;
t166 = -t54 * t85 + t61 * t80;
t165 = rSges(4,1) * t134 + rSges(4,2) * t135;
t164 = Icges(4,1) * t134 + t196;
t163 = Icges(4,2) * t135 + t197;
t162 = Icges(4,5) * t134 + Icges(4,6) * t135;
t161 = -V_base(4) * pkin(1) + t131 * t208 + V_base(2);
t160 = -t208 * V_base(5) + V_base(4) * t209 + V_base(3);
t159 = (Icges(4,3) * t152 + t147 * t162) * t102 + (-Icges(4,3) * t147 + t152 * t162) * t103 + t131 * (Icges(4,5) * t135 - Icges(4,6) * t134);
t96 = t205 * t152;
t158 = -pkin(3) * t188 + t131 * t96 + t161;
t157 = -t125 * t96 + t126 * t95 + t160;
t101 = pkin(4) * t135 - pkin(5) * t134;
t87 = t169 * t147;
t156 = t102 * t101 + (t172 - t87) * t131 + t170;
t88 = t169 * t152;
t155 = -t103 * t101 + t131 * t88 + t158;
t154 = -t102 * t88 + t103 * t87 + t157;
t71 = Icges(4,6) * t152 + t147 * t163;
t72 = -Icges(4,6) * t147 + t152 * t163;
t73 = Icges(4,5) * t152 + t147 * t164;
t74 = -Icges(4,5) * t147 + t152 * t164;
t98 = -Icges(4,2) * t134 + t196;
t99 = Icges(4,1) * t135 - t197;
t153 = (t134 * t74 + t135 * t72) * t103 + (t134 * t73 + t135 * t71) * t102 + (t134 * t99 + t135 * t98) * t131;
t141 = t152 ^ 2;
t139 = t147 ^ 2;
t138 = t145 ^ 2;
t136 = Icges(2,4) * t152;
t120 = rSges(2,1) * t152 - t147 * rSges(2,2);
t119 = t147 * rSges(2,1) + rSges(2,2) * t152;
t118 = t177 * t193;
t117 = Icges(2,1) * t152 - t198;
t116 = Icges(2,1) * t147 + t136;
t115 = -Icges(2,2) * t147 + t136;
t114 = Icges(2,2) * t152 + t198;
t110 = rSges(3,1) * t180 - t147 * rSges(3,3);
t109 = rSges(3,1) * t181 + rSges(3,3) * t152;
t107 = Icges(3,1) * t180 - t195;
t106 = Icges(3,1) * t181 + t194;
t105 = Icges(3,5) * t180 - Icges(3,3) * t147;
t104 = Icges(3,5) * t181 + Icges(3,3) * t152;
t100 = rSges(4,1) * t135 - rSges(4,2) * t134;
t81 = t85 ^ 2;
t76 = -t147 * rSges(4,3) + t152 * t165;
t75 = rSges(4,3) * t152 + t147 * t165;
t67 = V_base(5) * rSges(2,3) - t119 * t131 + t174;
t66 = t120 * t131 + V_base(2) + (-pkin(1) - rSges(2,3)) * V_base(4);
t64 = t119 * V_base(4) - t120 * V_base(5) + V_base(3);
t62 = t144 * t183 + t94 * t149;
t60 = t144 * t185 + t149 * t91;
t58 = t61 ^ 2;
t52 = qJD(6) * t85 + t80;
t51 = rSges(3,1) * t189 + (-t109 - t209) * t131 + t174;
t50 = -rSges(3,1) * t188 + t110 * t131 + t161;
t49 = t143 * t92 - t148 * t62;
t48 = t143 * t62 + t148 * t92;
t47 = t143 * t89 - t148 * t60;
t46 = t143 * t60 + t148 * t89;
t45 = t126 * t109 - t125 * t110 + t160;
t42 = qJD(6) * t61 + t54;
t41 = qJD(6) * t59 + t53;
t39 = t100 * t102 + (t172 - t75) * t131 + t170;
t38 = -t100 * t103 + t131 * t76 + t158;
t35 = -t102 * t76 + t103 * t75 + t157;
t34 = (-t111 * t185 - t134 * t78) * rSges(5,3) + t156;
t33 = (t111 * t183 + t134 * t79) * rSges(5,3) + t155;
t31 = (t147 * t79 - t152 * t78) * t135 * rSges(5,3) + t154;
t30 = rSges(6,2) * t167 + t156;
t29 = rSges(6,2) * t166 + t155;
t26 = rSges(6,2) * t168 + t154;
t24 = (t41 * t85 - t52 * t59) * rSges(7,3) + t167 * pkin(6) + t156;
t23 = (-t42 * t85 + t52 * t61) * rSges(7,3) + t166 * pkin(6) + t155;
t22 = (-t41 * t61 + t42 * t59) * rSges(7,3) + t168 * pkin(6) + t154;
t1 = t52 * t42 * (t200 * t61 + t201 * t48 + t203 * t49) + V_base(4) * ((-t147 * t115 + t117 * t152) * V_base(4) + (-t147 * t114 + t116 * t152) * V_base(5)) / 0.2e1 + (Icges(6,1) * t86 ^ 2 + Icges(6,3) * t130 * t138 + Icges(6,2) * t81) * t80 ^ 2 / 0.2e1 + (Icges(7,1) * t56 ^ 2 + Icges(7,2) * t55 ^ 2 + Icges(7,3) * t81) * t52 ^ 2 / 0.2e1 + (Icges(7,1) * t49 ^ 2 + Icges(7,2) * t48 ^ 2 + Icges(7,3) * t58) * t42 ^ 2 / 0.2e1 + V_base(5) * ((t115 * t152 + t147 * t117) * V_base(4) + (t114 * t152 + t147 * t116) * V_base(5)) / 0.2e1 + m(2) * (t64 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(3) * (t45 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(4) * (t35 ^ 2 + t38 ^ 2 + t39 ^ 2) / 0.2e1 + m(7) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) / 0.2e1 + m(6) * (t26 ^ 2 + t29 ^ 2 + t30 ^ 2) / 0.2e1 + m(5) * (t31 ^ 2 + t33 ^ 2 + t34 ^ 2) / 0.2e1 + V_base(5) * t131 * (Icges(2,5) * t147 + Icges(2,6) * t152) + V_base(4) * t131 * (Icges(2,5) * t152 - Icges(2,6) * t147) + t126 * ((-t147 * t105 + t107 * t180 + t141 * t193) * t126 + (-t147 * t104 + t106 * t180 + t118) * t125 + (t152 * t171 - t195) * t187) / 0.2e1 + (Icges(6,1) * t62 ^ 2 + Icges(6,2) * t58 + Icges(6,3) * t211) * t54 ^ 2 / 0.2e1 + (Icges(5,1) * t94 ^ 2 + Icges(5,2) * t211 + t141 * t191) * t79 ^ 2 / 0.2e1 + (Icges(5,1) * t91 ^ 2 + Icges(5,2) * t212 + t139 * t191) * t78 ^ 2 / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + t80 * t54 * (t173 * t92 + t202 * t61 + t204 * t62) + t78 * t79 * (t94 * Icges(5,1) * t91 + Icges(5,2) * t206 + t177 * t191) + t125 * ((t105 * t152 + t107 * t181 + t118) * t126 + (t104 * t152 + t106 * t181 + t139 * t193) * t125 + (t147 * t171 + t194) * t187) / 0.2e1 + t102 * (t153 * t147 + t159 * t152) / 0.2e1 + t103 * (-t159 * t147 + t153 * t152) / 0.2e1 + ((Icges(7,1) * t47 * t49 + Icges(7,2) * t46 * t48 + Icges(7,3) * t207) * t42 + (t200 * t59 + t201 * t46 + t203 * t47) * t52 + (Icges(7,1) * t47 ^ 2 / 0.2e1 + Icges(7,2) * t46 ^ 2 / 0.2e1 + Icges(7,3) * t213) * t41) * t41 + ((Icges(6,1) * t60 * t62 + Icges(6,2) * t207 + Icges(6,3) * t206) * t54 + (t173 * t89 + t202 * t59 + t204 * t60) * t80 + (Icges(6,1) * t60 ^ 2 / 0.2e1 + Icges(6,2) * t213 + Icges(6,3) * t212 / 0.2e1) * t53) * t53 + ((Icges(5,3) * t134 ^ 2 / 0.2e1 + (Icges(5,1) * t150 ^ 2 + Icges(5,2) * t138) * t130 / 0.2e1) * t111 + ((-t147 * t190 + t192 * t89 + t199 * t91) * t78 + (-t152 * t190 + t192 * t92 + t199 * t94) * t79) * t135) * t111 + ((-t134 * t72 + t135 * t74) * t103 + (-t134 * t71 + t135 * t73) * t102 + ((-Icges(3,2) * t180 + t107) * t126 + (-Icges(3,2) * t181 + t106) * t125) * t151 + (Icges(3,2) * t146 ^ 2 + Icges(3,1) * t140 - t134 * t98 + t135 * t99 + Icges(2,3)) * t131) * t131 / 0.2e1;
T = t1;
