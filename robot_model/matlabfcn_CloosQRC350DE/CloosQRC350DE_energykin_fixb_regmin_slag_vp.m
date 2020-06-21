% Calculate minimal parameter regressor of fixed base kinetic energy for
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
% 
% Output:
% T_reg [1x36]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = CloosQRC350DE_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:33
% EndTime: 2020-06-19 21:39:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (396->43), mult. (798->108), div. (0->0), fcn. (619->10), ass. (0->95)
unknown=NaN(1,36);
t1 = qJD(1) ^ 2;
t3 = cos(qJ(2));
t4 = t3 ^ 2;
t8 = sin(qJ(2));
t12 = t8 * qJD(1);
t14 = qJD(2) ^ 2;
t16 = t1 * pkin(2);
t19 = cos(qJ(3));
t22 = sin(qJ(3));
t25 = -t19 * t3 * qJD(1) + t22 * t8 * qJD(1);
t26 = t25 ^ 2;
t29 = t22 * t3 * qJD(1);
t31 = t19 * t8 * qJD(1);
t32 = t29 + t31;
t34 = qJD(2) + qJD(3);
t37 = t34 ^ 2;
t39 = t19 * qJD(2);
t42 = qJD(1) * pkin(2);
t43 = t12 * pkin(3);
t44 = -t42 - t43;
t47 = t22 * qJD(2);
t52 = cos(qJ(4));
t54 = sin(qJ(4));
t56 = t52 * t25 - t54 * t34;
t57 = t56 ^ 2;
t59 = t54 * t25;
t60 = t52 * t34;
t61 = -t59 - t60;
t63 = t29 + t31 + qJD(4);
t66 = t63 ^ 2;
t70 = t47 * pkin(3) - t34 * pkin(5);
t71 = t54 * t70;
t74 = -t32 * pkin(4) + t25 * pkin(5) - t42 - t43;
t75 = t52 * t74;
t76 = -t71 - t75;
t80 = t39 * pkin(3) + t34 * pkin(4);
t85 = t52 * t70 - t54 * t74;
t89 = cos(qJ(5));
t91 = sin(qJ(5));
t93 = t89 * t56 + t91 * t63;
t94 = t93 ^ 2;
t96 = t91 * t56;
t97 = t89 * t63;
t98 = -t96 + t97;
t100 = t59 + t60 + qJD(5);
t103 = t100 ^ 2;
t107 = t89 * t80 - t91 * t85;
t111 = t89 * t85;
t112 = t91 * t80;
t118 = pkin(7) * qJ(5) - qJ(6);
t119 = cos(t118);
t121 = sin(t118);
t123 = -t121 * t100 - t119 * t93;
t124 = t123 ^ 2;
t128 = t119 * t100 - t121 * t93;
t131 = -pkin(7) * qJD(5) + qJD(6) - t96 + t97;
t134 = t131 ^ 2;
t137 = -t100 * pkin(6) + t111 + t112;
t140 = t93 * pkin(6) + t71 + t75;
unknown(1,1) = t1 / 0.2e1;
unknown(1,2) = t4 * t1 / 0.2e1;
unknown(1,3) = -t3 * t1 * t8;
unknown(1,4) = -t3 * qJD(1) * qJD(2);
unknown(1,5) = t12 * qJD(2);
unknown(1,6) = t14 / 0.2e1;
unknown(1,7) = t8 * t16;
unknown(1,8) = t16 * t3;
unknown(1,9) = t26 / 0.2e1;
unknown(1,10) = t25 * t32;
unknown(1,11) = t25 * t34;
unknown(1,12) = t32 * t34;
unknown(1,13) = t37 / 0.2e1;
unknown(1,14) = t39 * pkin(3) * t34 - t44 * t32;
unknown(1,15) = -t47 * pkin(3) * t34 + t44 * t25;
unknown(1,16) = t57 / 0.2e1;
unknown(1,17) = t56 * t61;
unknown(1,18) = t56 * t63;
unknown(1,19) = t61 * t63;
unknown(1,20) = t66 / 0.2e1;
unknown(1,21) = -t80 * t61 + t76 * t63;
unknown(1,22) = t80 * t56 - t85 * t63;
unknown(1,23) = t94 / 0.2e1;
unknown(1,24) = t93 * t98;
unknown(1,25) = t93 * t100;
unknown(1,26) = t98 * t100;
unknown(1,27) = t103 / 0.2e1;
unknown(1,28) = t107 * t100 + t76 * t98;
unknown(1,29) = -(t111 + t112) * t100 - t76 * t93;
unknown(1,30) = t124 / 0.2e1;
unknown(1,31) = t123 * t128;
unknown(1,32) = t123 * t131;
unknown(1,33) = t128 * t131;
unknown(1,34) = t134 / 0.2e1;
unknown(1,35) = (t119 * t140 - t121 * t137) * t131 - t107 * t128;
unknown(1,36) = -(-t119 * t137 - t121 * t140) * t131 + t107 * t123;
T_reg = unknown;
