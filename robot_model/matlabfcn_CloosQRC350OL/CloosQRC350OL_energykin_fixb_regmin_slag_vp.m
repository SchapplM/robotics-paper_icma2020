% Calculate minimal parameter regressor of fixed base kinetic energy for
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
% 
% Output:
% T_reg [1x36]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = CloosQRC350OL_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:18:55
% EndTime: 2020-06-20 08:18:55
% DurationCPUTime: 0.32s
% Computational Cost: add. (373->43), mult. (775->99), div. (0->0), fcn. (619->10), ass. (0->40)
t126 = sin(qJ(3));
t127 = sin(qJ(2));
t130 = cos(qJ(3));
t131 = cos(qJ(2));
t114 = (t126 * t131 + t127 * t130) * qJD(1);
t132 = qJD(1) ^ 2;
t141 = t132 / 0.2e1;
t140 = cos(qJ(5));
t139 = pkin(3) * qJD(2);
t138 = t131 * t132;
t115 = (-t126 * t127 + t130 * t131) * qJD(1);
t119 = (pkin(3) * t127 + pkin(2)) * qJD(1);
t107 = t114 * pkin(4) + t115 * pkin(5) + t119;
t121 = qJD(2) + qJD(3);
t135 = t126 * t139;
t117 = -t121 * pkin(5) + t135;
t125 = sin(qJ(4));
t129 = cos(qJ(4));
t102 = -t125 * t107 + t129 * t117;
t134 = t130 * t139;
t118 = t121 * pkin(4) + t134;
t124 = sin(qJ(5));
t137 = t140 * t102 + t124 * t118;
t136 = qJD(1) * qJD(2);
t110 = t129 * t115 - t125 * t121;
t112 = -qJD(4) + t114;
t105 = t124 * t110 + t140 * t112;
t101 = t129 * t107 + t125 * t117;
t109 = t125 * t115 + t129 * t121;
t128 = cos(qJ(6));
t123 = sin(qJ(6));
t108 = qJD(5) + t109;
t106 = t140 * t110 - t124 * t112;
t103 = -qJD(6) + t105;
t99 = -t124 * t102 + t140 * t118;
t98 = -t128 * t106 + t123 * t108;
t97 = t123 * t106 + t128 * t108;
t96 = -t108 * pkin(6) + t137;
t95 = t106 * pkin(6) + t101;
t1 = [t141, t131 ^ 2 * t141, -t127 * t138, t131 * t136, -t127 * t136, qJD(2) ^ 2 / 0.2e1, t132 * pkin(2) * t127, pkin(2) * t138, t115 ^ 2 / 0.2e1, -t115 * t114, t115 * t121, -t114 * t121, t121 ^ 2 / 0.2e1, t119 * t114 + t121 * t134, t119 * t115 - t121 * t135, t110 ^ 2 / 0.2e1, -t110 * t109, -t110 * t112, t109 * t112, t112 ^ 2 / 0.2e1, t101 * t112 + t118 * t109, t102 * t112 + t118 * t110, t106 ^ 2 / 0.2e1, -t106 * t105, t106 * t108, -t105 * t108, t108 ^ 2 / 0.2e1, t101 * t105 + t99 * t108, t101 * t106 - t137 * t108, t98 ^ 2 / 0.2e1, t98 * t97, -t98 * t103, -t97 * t103, t103 ^ 2 / 0.2e1, -(t123 * t96 + t128 * t95) * t103 - t99 * t97, (t123 * t95 - t128 * t96) * t103 + t99 * t98;];
T_reg = t1;
