% Calculate kinetic energy for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
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

function T = CloosQRC350DE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:24
% EndTime: 2020-06-19 21:37:24
% DurationCPUTime: 0.18s
% Computational Cost: add. (706->117), mult. (1349->201), div. (0->0), fcn. (1004->10), ass. (0->68)
t1 = sin(qJ(3));
t2 = t1 * qJD(2);
t3 = cos(qJ(2));
t5 = t1 * t3 * qJD(1);
t6 = cos(qJ(3));
t7 = sin(qJ(2));
t9 = t6 * t7 * qJD(1);
t10 = t5 + t9;
t12 = qJD(2) + qJD(3);
t17 = t6 * qJD(2);
t22 = t1 * t7 * qJD(1) - t6 * t3 * qJD(1);
t54 = qJD(1) * pkin(2);
t55 = t7 * qJD(1);
t56 = t55 * pkin(3);
t57 = -t54 - t56;
t62 = cos(qJ(5));
t63 = cos(qJ(4));
t66 = t2 * pkin(3) - t12 * pkin(5);
t68 = sin(qJ(4));
t71 = -t10 * pkin(4) + t22 * pkin(5) - t54 - t56;
t73 = t63 * t66 - t68 * t71;
t74 = t62 * t73;
t75 = sin(qJ(5));
t78 = t17 * pkin(3) + t12 * pkin(4);
t79 = t75 * t78;
t80 = t74 + t79;
t81 = t80 ^ 2;
t84 = t62 * t78 - t75 * t73;
t85 = t84 ^ 2;
t86 = t68 * t66;
t87 = t63 * t71;
t88 = t86 + t87;
t89 = t88 ^ 2;
t95 = -t68 * t12 + t63 * t22;
t97 = t5 + t9 + qJD(4);
t99 = t62 * t95 + t75 * t97;
t101 = t75 * t95;
t102 = t62 * t97;
t103 = -t101 + t102;
t108 = pkin(7) * qJ(5) - qJ(6);
t109 = cos(t108);
t110 = t68 * t22;
t111 = t63 * t12;
t112 = t110 + t111 + qJD(5);
t114 = -t112 * pkin(6) + t74 + t79;
t116 = sin(t108);
t118 = t99 * pkin(6) + t86 + t87;
t120 = -t109 * t114 - t116 * t118;
t123 = t109 * t112 - t116 * t99;
t126 = -pkin(7) * qJD(5) + qJD(6) - t101 + t102;
t132 = t109 * t118 - t116 * t114;
t135 = -t109 * t99 - t116 * t112;
t140 = t120 ^ 2;
t141 = t132 ^ 2;
t157 = -t110 - t111;
t166 = t2 * pkin(3) * (-t12 * mrSges(4,2) + t10 * mrSges(4,3)) + t17 * pkin(3) * (t12 * mrSges(4,1) - t22 * mrSges(4,3)) + t10 * (Ifges(4,4) * t22 + Ifges(4,2) * t10 + Ifges(4,6) * t12) / 0.2e1 + t12 * (Ifges(4,5) * t22 + Ifges(4,6) * t10 + Ifges(4,3) * t12) / 0.2e1 + t22 * (Ifges(4,1) * t22 + Ifges(4,4) * t10 + Ifges(4,5) * t12) / 0.2e1 + qJD(2) * (-Ifges(3,5) * t3 * qJD(1) + Ifges(3,6) * t7 * qJD(1) + Ifges(3,3) * qJD(2)) / 0.2e1 + t57 * (-t10 * mrSges(4,1) + t22 * mrSges(4,2)) + m(6) * (t81 + t85 + t89) / 0.2e1 + t88 * (-t103 * mrSges(6,1) + t99 * mrSges(6,2)) + t120 * (-t126 * mrSges(7,2) + t123 * mrSges(7,3)) + t132 * (t126 * mrSges(7,1) - t135 * mrSges(7,3)) + m(7) * (t140 + t141 + t85) / 0.2e1 + t84 * (-t123 * mrSges(7,1) + t135 * mrSges(7,2)) + t80 * (-t112 * mrSges(6,2) + t103 * mrSges(6,3)) + t84 * (t112 * mrSges(6,1) - t99 * mrSges(6,3)) + t73 * (-t97 * mrSges(5,2) + t157 * mrSges(5,3)) - t88 * (t97 * mrSges(5,1) - t95 * mrSges(5,3));
t167 = t73 ^ 2;
t168 = t78 ^ 2;
t176 = t1 ^ 2;
t177 = qJD(2) ^ 2;
t179 = pkin(3) ^ 2;
t181 = t6 ^ 2;
t184 = t57 ^ 2;
t188 = qJD(1) ^ 2;
t246 = pkin(2) ^ 2;
t249 = t3 * qJD(1);
t270 = m(5) * (t167 + t89 + t168) / 0.2e1 + t78 * (-t157 * mrSges(5,1) + t95 * mrSges(5,2)) + m(4) * (t176 * t177 * t179 + t181 * t177 * t179 + t184) / 0.2e1 + t188 * Ifges(2,3) / 0.2e1 + t123 * (Ifges(7,4) * t135 + Ifges(7,2) * t123 + Ifges(7,6) * t126) / 0.2e1 + t126 * (Ifges(7,5) * t135 + Ifges(7,6) * t123 + Ifges(7,3) * t126) / 0.2e1 + t135 * (Ifges(7,1) * t135 + Ifges(7,4) * t123 + Ifges(7,5) * t126) / 0.2e1 + t99 * (Ifges(6,1) * t99 + Ifges(6,4) * t103 + Ifges(6,5) * t112) / 0.2e1 + t103 * (Ifges(6,4) * t99 + Ifges(6,2) * t103 + Ifges(6,6) * t112) / 0.2e1 + t112 * (Ifges(6,5) * t99 + Ifges(6,6) * t103 + Ifges(6,3) * t112) / 0.2e1 + t95 * (Ifges(5,1) * t95 + Ifges(5,4) * t157 + Ifges(5,5) * t97) / 0.2e1 + t157 * (Ifges(5,4) * t95 + Ifges(5,2) * t157 + Ifges(5,6) * t97) / 0.2e1 + t97 * (Ifges(5,5) * t95 + Ifges(5,6) * t157 + Ifges(5,3) * t97) / 0.2e1 + m(3) * t188 * t246 / 0.2e1 - t249 * (-Ifges(3,1) * t3 * qJD(1) + Ifges(3,4) * t7 * qJD(1) + Ifges(3,5) * qJD(2)) / 0.2e1 + t55 * (-Ifges(3,4) * t3 * qJD(1) + Ifges(3,2) * t7 * qJD(1) + Ifges(3,6) * qJD(2)) / 0.2e1 - t54 * (-t55 * mrSges(3,1) - t249 * mrSges(3,2));
t271 = t166 + t270;
T = t271;
