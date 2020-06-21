% Calculate kinetic energy for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = CloosQRC350OL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:00:01
% EndTime: 2020-06-20 08:00:03
% DurationCPUTime: 1.54s
% Computational Cost: add. (668->94), mult. (1311->159), div. (0->0), fcn. (1004->10), ass. (0->38)
t102 = sin(qJ(4));
t107 = cos(qJ(4));
t103 = sin(qJ(3));
t104 = sin(qJ(2));
t108 = cos(qJ(3));
t109 = cos(qJ(2));
t91 = (-t103 * t109 - t104 * t108) * qJD(1);
t92 = (-t103 * t104 + t108 * t109) * qJD(1);
t96 = (pkin(3) * t104 + pkin(2)) * qJD(1);
t85 = -pkin(4) * t91 + pkin(5) * t92 + t96;
t113 = pkin(3) * qJD(2);
t98 = qJD(2) + qJD(3);
t94 = -pkin(5) * t98 + t103 * t113;
t79 = t102 * t94 + t107 * t85;
t115 = t79 ^ 2;
t101 = sin(qJ(5));
t106 = cos(qJ(5));
t81 = -t102 * t85 + t107 * t94;
t95 = pkin(4) * t98 + t108 * t113;
t77 = t101 * t95 + t106 * t81;
t88 = -t102 * t98 + t107 * t92;
t90 = qJD(4) + t91;
t83 = -t101 * t88 + t106 * t90;
t87 = -t102 * t92 - t107 * t98;
t105 = cos(qJ(6));
t100 = sin(qJ(6));
t86 = qJD(5) - t87;
t84 = t101 * t90 + t106 * t88;
t82 = qJD(6) + t83;
t76 = -t101 * t81 + t106 * t95;
t75 = t76 ^ 2;
t74 = t100 * t86 - t105 * t84;
t73 = t100 * t84 + t105 * t86;
t72 = -pkin(6) * t86 + t77;
t71 = pkin(6) * t84 + t79;
t70 = t100 * t71 - t105 * t72;
t69 = t100 * t72 + t105 * t71;
t1 = m(4) * (t96 ^ 2 + (t103 ^ 2 + t108 ^ 2) * pkin(3) ^ 2 * qJD(2) ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t95 ^ 2 + t115) / 0.2e1 + Ifges(4,3) * t98 ^ 2 / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t75) / 0.2e1 + m(6) * (t77 ^ 2 + t115 + t75) / 0.2e1 + (-t79 * mrSges(5,1) - t81 * mrSges(5,2) + Ifges(5,3) * t90 / 0.2e1) * t90 + (t76 * mrSges(6,1) - t77 * mrSges(6,2) + Ifges(6,3) * t86 / 0.2e1) * t86 + (t69 * mrSges(7,1) - t70 * mrSges(7,2) + Ifges(7,3) * t82 / 0.2e1) * t82 + (t96 * mrSges(4,2) + Ifges(4,5) * t98 + Ifges(4,1) * t92 / 0.2e1) * t92 + (t95 * mrSges(5,2) + t79 * mrSges(5,3) + Ifges(5,5) * t90 + Ifges(5,1) * t88 / 0.2e1) * t88 + (t79 * mrSges(6,2) - t76 * mrSges(6,3) + Ifges(6,5) * t86 + Ifges(6,1) * t84 / 0.2e1) * t84 + (t76 * mrSges(7,2) - t69 * mrSges(7,3) + Ifges(7,5) * t82 + Ifges(7,1) * t74 / 0.2e1) * t74 + (-t96 * mrSges(4,1) + Ifges(4,4) * t92 + Ifges(4,6) * t98 + Ifges(4,2) * t91 / 0.2e1) * t91 + (-t95 * mrSges(5,1) + t81 * mrSges(5,3) + Ifges(5,4) * t88 + Ifges(5,6) * t90 + Ifges(5,2) * t87 / 0.2e1) * t87 + (-t79 * mrSges(6,1) + t77 * mrSges(6,3) + Ifges(6,4) * t84 + Ifges(6,6) * t86 + Ifges(6,2) * t83 / 0.2e1) * t83 + (-t76 * mrSges(7,1) + t70 * mrSges(7,3) + Ifges(7,4) * t74 + Ifges(7,6) * t82 + Ifges(7,2) * t73 / 0.2e1) * t73 + (Ifges(3,3) * qJD(2) / 0.2e1 + (Ifges(3,5) * t109 - Ifges(3,6) * t104) * qJD(1) + (t103 * (-mrSges(4,2) * t98 + mrSges(4,3) * t91) + t108 * (mrSges(4,1) * t98 - t92 * mrSges(4,3))) * pkin(3)) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * pkin(2) ^ 2 / 0.2e1 + (pkin(2) * mrSges(3,2) + Ifges(3,1) * t109 / 0.2e1) * t109 + (pkin(2) * mrSges(3,1) - Ifges(3,4) * t109 + Ifges(3,2) * t104 / 0.2e1) * t104) * qJD(1) ^ 2;
T = t1;
