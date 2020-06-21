% Calculate kinetic energy for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = CloosQRC350DE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:23
% EndTime: 2020-06-19 21:37:24
% DurationCPUTime: 0.40s
% Computational Cost: add. (2783->186), mult. (3623->262), div. (0->0), fcn. (2936->12), ass. (0->99)
t1 = cos(qJ(1));
t3 = sin(qJ(1));
t5 = t1 * V_base(4) - t3 * V_base(5);
t7 = t3 * V_base(4);
t8 = t1 * V_base(5);
t9 = t7 + t8;
t11 = V_base(6) - qJD(1);
t29 = V_base(5) * pkin(1) + V_base(1);
t32 = -V_base(4) * pkin(1) + V_base(2);
t34 = t1 * t29 - t3 * t32;
t39 = t34 ^ 2;
t40 = t3 * t29;
t41 = t1 * t32;
t42 = t40 + t41;
t43 = t42 ^ 2;
t44 = V_base(3) ^ 2;
t82 = t5 * (Ifges(2,1) * t5 + Ifges(2,4) * t9 + Ifges(2,5) * t11) / 0.2e1 + t9 * (Ifges(2,4) * t5 + Ifges(2,2) * t9 + Ifges(2,6) * t11) / 0.2e1 + t11 * (Ifges(2,5) * t5 + Ifges(2,6) * t9 + Ifges(2,3) * t11) / 0.2e1 + t34 * (-t11 * mrSges(2,2) + t9 * mrSges(2,3)) + m(2) * (t39 + t43 + t44) / 0.2e1 + V_base(3) * (-t9 * mrSges(2,1) + t5 * mrSges(2,2)) + t42 * (t11 * mrSges(2,1) - t5 * mrSges(2,3)) + V_base(2) * (V_base(6) * mrSges(1,1) - V_base(4) * mrSges(1,3)) + V_base(4) * (Ifges(1,1) * V_base(4) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6)) / 0.2e1 + V_base(5) * (Ifges(1,4) * V_base(4) + Ifges(1,2) * V_base(5) + Ifges(1,6) * V_base(6)) / 0.2e1 + V_base(6) * (Ifges(1,5) * V_base(4) + Ifges(1,6) * V_base(5) + Ifges(1,3) * V_base(6)) / 0.2e1 + V_base(1) * (-V_base(6) * mrSges(1,2) + V_base(5) * mrSges(1,3));
t87 = V_base(1) ^ 2;
t88 = V_base(2) ^ 2;
t92 = cos(qJ(5));
t93 = cos(qJ(4));
t94 = cos(qJ(3));
t95 = sin(qJ(2));
t97 = cos(qJ(2));
t99 = t97 * t11 + t95 * t5;
t101 = sin(qJ(3));
t104 = -t95 * t11 + t97 * t5;
t106 = t101 * t104 + t94 * t99;
t108 = sin(qJ(4));
t109 = t7 + t8 + qJD(2) + qJD(3);
t111 = t93 * t106 - t108 * t109;
t113 = sin(qJ(5));
t114 = t101 * t99;
t115 = t94 * t104;
t116 = -t114 + t115 + qJD(4);
t118 = t92 * t111 + t113 * t116;
t120 = t113 * t111;
t121 = t92 * t116;
t122 = -t120 + t121;
t124 = t108 * t106;
t125 = t93 * t109;
t126 = t124 + t125 + qJD(5);
t144 = -t124 - t125;
t152 = -t9 * pkin(2) + V_base(3);
t154 = t97 * t152 + t95 * t34;
t155 = t94 * t154;
t156 = t97 * t34;
t157 = t95 * t152;
t158 = t7 + t8 + qJD(2);
t160 = t158 * pkin(3) + t156 - t157;
t161 = t101 * t160;
t163 = -t109 * pkin(5) + t155 + t161;
t165 = t11 * pkin(2);
t166 = t104 * pkin(3);
t168 = -t114 + t115;
t170 = -t168 * pkin(4) + t106 * pkin(5) + t165 - t166 + t40 + t41;
t172 = -t108 * t170 + t93 * t163;
t177 = t108 * t163;
t178 = t93 * t170;
t179 = -t177 - t178;
t196 = t172 ^ 2;
t197 = t179 ^ 2;
t198 = t101 * t154;
t199 = t94 * t160;
t201 = t109 * pkin(4) - t198 + t199;
t202 = t201 ^ 2;
t210 = V_base(3) * (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2)) + m(1) * (t87 + t88 + t44) / 0.2e1 + t118 * (Ifges(6,1) * t118 + Ifges(6,4) * t122 + Ifges(6,5) * t126) / 0.2e1 + t122 * (Ifges(6,4) * t118 + Ifges(6,2) * t122 + Ifges(6,6) * t126) / 0.2e1 + t126 * (Ifges(6,5) * t118 + Ifges(6,6) * t122 + Ifges(6,3) * t126) / 0.2e1 + t111 * (Ifges(5,1) * t111 + Ifges(5,4) * t144 + Ifges(5,5) * t116) / 0.2e1 + t172 * (-t116 * mrSges(5,2) + t144 * mrSges(5,3)) + t179 * (t116 * mrSges(5,1) - t111 * mrSges(5,3)) + t144 * (Ifges(5,4) * t111 + Ifges(5,2) * t144 + Ifges(5,6) * t116) / 0.2e1 + t116 * (Ifges(5,5) * t111 + Ifges(5,6) * t144 + Ifges(5,3) * t116) / 0.2e1 + m(5) * (t196 + t197 + t202) / 0.2e1 + t201 * (-t144 * mrSges(5,1) + t111 * mrSges(5,2));
t212 = t40 + t41 + t165 - t166;
t229 = t155 + t161;
t230 = t229 ^ 2;
t231 = -t198 + t199;
t232 = t231 ^ 2;
t233 = t212 ^ 2;
t255 = t156 - t157;
t266 = t40 + t41 + t165;
t271 = t154 ^ 2;
t272 = t255 ^ 2;
t273 = t266 ^ 2;
t277 = t212 * (-t168 * mrSges(4,1) + t106 * mrSges(4,2)) + t168 * (Ifges(4,4) * t106 + Ifges(4,2) * t168 + Ifges(4,6) * t109) / 0.2e1 + t109 * (Ifges(4,5) * t106 + Ifges(4,6) * t168 + Ifges(4,3) * t109) / 0.2e1 + m(4) * (t230 + t232 + t233) / 0.2e1 + t106 * (Ifges(4,1) * t106 + Ifges(4,4) * t168 + Ifges(4,5) * t109) / 0.2e1 + t229 * (-t109 * mrSges(4,2) + t168 * mrSges(4,3)) + t231 * (t109 * mrSges(4,1) - t106 * mrSges(4,3)) + t154 * (-t158 * mrSges(3,2) + t104 * mrSges(3,3)) + t255 * (t158 * mrSges(3,1) - t99 * mrSges(3,3)) + t99 * (Ifges(3,1) * t99 + Ifges(3,4) * t104 + Ifges(3,5) * t158) / 0.2e1 + t266 * (-t104 * mrSges(3,1) + t99 * mrSges(3,2)) + m(3) * (t271 + t272 + t273) / 0.2e1;
t290 = t92 * t172;
t291 = t113 * t201;
t292 = t290 + t291;
t293 = t292 ^ 2;
t296 = -t113 * t172 + t92 * t201;
t297 = t296 ^ 2;
t306 = pkin(7) * qJ(5) - qJ(6);
t307 = sin(t306);
t309 = cos(t306);
t311 = -t307 * t118 + t309 * t126;
t314 = -t309 * t118 - t307 * t126;
t318 = -pkin(7) * qJD(5) + qJD(6) - t120 + t121;
t336 = -t126 * pkin(6) + t290 + t291;
t339 = t118 * pkin(6) + t177 + t178;
t341 = -t307 * t339 - t309 * t336;
t348 = -t307 * t336 + t309 * t339;
t353 = t341 ^ 2;
t354 = t348 ^ 2;
t370 = t158 * (Ifges(3,5) * t99 + Ifges(3,6) * t104 + Ifges(3,3) * t158) / 0.2e1 + t104 * (Ifges(3,4) * t99 + Ifges(3,2) * t104 + Ifges(3,6) * t158) / 0.2e1 + m(6) * (t293 + t297 + t197) / 0.2e1 - t179 * (-t122 * mrSges(6,1) + t118 * mrSges(6,2)) + t311 * (Ifges(7,4) * t314 + Ifges(7,2) * t311 + Ifges(7,6) * t318) / 0.2e1 + t318 * (Ifges(7,5) * t314 + Ifges(7,6) * t311 + Ifges(7,3) * t318) / 0.2e1 + t314 * (Ifges(7,1) * t314 + Ifges(7,4) * t311 + Ifges(7,5) * t318) / 0.2e1 + t341 * (-t318 * mrSges(7,2) + t311 * mrSges(7,3)) + t348 * (t318 * mrSges(7,1) - t314 * mrSges(7,3)) + m(7) * (t353 + t354 + t297) / 0.2e1 + t296 * (-t311 * mrSges(7,1) + t314 * mrSges(7,2)) + t292 * (-t126 * mrSges(6,2) + t122 * mrSges(6,3)) + t296 * (t126 * mrSges(6,1) - t118 * mrSges(6,3));
t372 = t82 + t210 + t277 + t370;
T = t372;
