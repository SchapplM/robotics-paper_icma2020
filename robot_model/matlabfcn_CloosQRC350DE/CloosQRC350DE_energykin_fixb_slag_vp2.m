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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
% OptimizationMode: 2
% StartTime: 2020-06-23 20:59:51
% EndTime: 2020-06-23 20:59:52
% DurationCPUTime: 0.92s
% Computational Cost: add. (399->72), mult. (798->140), div. (0->0), fcn. (548->10), ass. (0->41)
t127 = sin(qJ(3));
t128 = sin(qJ(2));
t131 = cos(qJ(3));
t132 = cos(qJ(2));
t136 = t132 * qJD(1);
t113 = t131 * t128 * qJD(1) + t127 * t136;
t114 = (t127 * t128 - t131 * t132) * qJD(1);
t122 = -t128 * pkin(3) - pkin(2);
t107 = -t113 * pkin(4) + t114 * pkin(5) + t122 * qJD(1);
t124 = qJD(2) + qJD(3);
t137 = pkin(3) * qJD(2);
t116 = -t124 * pkin(5) + t127 * t137;
t126 = sin(qJ(4));
t130 = cos(qJ(4));
t101 = t130 * t107 + t126 * t116;
t139 = t101 ^ 2;
t138 = m(4) / 0.2e1;
t103 = -t126 * t107 + t130 * t116;
t117 = t124 * pkin(4) + t131 * t137;
t125 = sin(qJ(5));
t129 = cos(qJ(5));
t99 = t129 * t103 + t125 * t117;
t110 = t130 * t114 - t126 * t124;
t112 = qJD(4) + t113;
t105 = -t125 * t110 + t129 * t112;
t109 = -t126 * t114 - t130 * t124;
t123 = pkin(7) * qJ(5) - qJ(6);
t121 = cos(t123);
t120 = sin(t123);
t108 = qJD(5) - t109;
t106 = t129 * t110 + t125 * t112;
t104 = -pkin(7) * qJD(5) + qJD(6) + t105;
t98 = -t125 * t103 + t129 * t117;
t97 = t98 ^ 2;
t96 = -t121 * t106 - t120 * t108;
t95 = -t120 * t106 + t121 * t108;
t94 = -t108 * pkin(6) + t99;
t93 = t106 * pkin(6) + t101;
t92 = -t120 * t93 - t121 * t94;
t91 = -t120 * t94 + t121 * t93;
t1 = t104 ^ 2 * Ifges(7,3) / 0.2e1 + t105 ^ 2 * Ifges(6,2) / 0.2e1 + t106 ^ 2 * Ifges(6,1) / 0.2e1 + t108 ^ 2 * Ifges(6,3) / 0.2e1 + t109 ^ 2 * Ifges(5,2) / 0.2e1 + t110 ^ 2 * Ifges(5,1) / 0.2e1 + t95 ^ 2 * Ifges(7,2) / 0.2e1 + t96 ^ 2 * Ifges(7,1) / 0.2e1 + m(7) * (t91 ^ 2 + t92 ^ 2 + t97) / 0.2e1 + m(6) * (t99 ^ 2 + t139 + t97) / 0.2e1 + Ifges(4,1) * t114 ^ 2 / 0.2e1 + t112 ^ 2 * Ifges(5,3) / 0.2e1 + m(5) * (t103 ^ 2 + t117 ^ 2 + t139) / 0.2e1 + (Ifges(4,4) * t114 + Ifges(4,2) * t113 / 0.2e1) * t113 + (-t91 * t96 + t92 * t95) * mrSges(7,3) + (t101 * t110 + t103 * t109) * mrSges(5,3) + (t101 * t106 - t99 * t108) * mrSges(6,2) + (Ifges(4,5) * t114 + Ifges(4,6) * t113 + Ifges(4,3) * t124 / 0.2e1) * t124 + (-Ifges(3,5) * t136 + (t127 * (-t124 * mrSges(4,2) + t113 * mrSges(4,3)) + t131 * (t124 * mrSges(4,1) - t114 * mrSges(4,3))) * pkin(3) + ((t127 ^ 2 + t131 ^ 2) * pkin(3) ^ 2 * t138 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + (t122 * (-t113 * mrSges(4,1) + t114 * mrSges(4,2)) + (t132 ^ 2 * Ifges(3,1) / 0.2e1 + t128 ^ 2 * Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + t122 ^ 2 * t138 + (t128 * mrSges(3,1) + m(3) * pkin(2) / 0.2e1) * pkin(2)) * qJD(1)) * qJD(1);
T = t1;
