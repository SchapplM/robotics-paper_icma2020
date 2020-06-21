% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(6*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = CloosQRC350DE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:33
% EndTime: 2020-06-19 21:39:33
% DurationCPUTime: 0.28s
% Computational Cost: add. (711->55), mult. (1461->159), div. (0->0), fcn. (1075->10), ass. (0->139)
unknown=NaN(1,60);
t1 = qJD(1) ^ 2;
t3 = cos(qJ(2));
t4 = t3 ^ 2;
t8 = sin(qJ(2));
t12 = t8 ^ 2;
t15 = t8 * qJD(1);
t17 = qJD(2) ^ 2;
t19 = t1 * pkin(2);
t22 = pkin(2) ^ 2;
t25 = cos(qJ(3));
t28 = sin(qJ(3));
t31 = -t25 * t3 * qJD(1) + t28 * t8 * qJD(1);
t32 = t31 ^ 2;
t35 = t28 * t3 * qJD(1);
t37 = t25 * t8 * qJD(1);
t38 = t35 + t37;
t40 = qJD(2) + qJD(3);
t42 = t38 ^ 2;
t45 = t40 ^ 2;
t47 = t25 * qJD(2);
t50 = qJD(1) * pkin(2);
t51 = t15 * pkin(3);
t52 = -t50 - t51;
t55 = t28 * qJD(2);
t65 = t28 ^ 2;
t67 = pkin(3) ^ 2;
t69 = t25 ^ 2;
t72 = t52 ^ 2;
t74 = cos(qJ(4));
t76 = sin(qJ(4));
t78 = t74 * t31 - t76 * t40;
t79 = t78 ^ 2;
t81 = t76 * t31;
t82 = t74 * t40;
t83 = -t81 - t82;
t85 = t35 + t37 + qJD(4);
t87 = t83 ^ 2;
t90 = t85 ^ 2;
t94 = t55 * pkin(3) - t40 * pkin(5);
t95 = t76 * t94;
t98 = -t38 * pkin(4) + t31 * pkin(5) - t50 - t51;
t99 = t74 * t98;
t100 = -t95 - t99;
t104 = t47 * pkin(3) + t40 * pkin(4);
t109 = t74 * t94 - t76 * t98;
t116 = t109 ^ 2;
t117 = t100 ^ 2;
t118 = t104 ^ 2;
t120 = cos(qJ(5));
t122 = sin(qJ(5));
t124 = t120 * t78 + t122 * t85;
t125 = t124 ^ 2;
t127 = t122 * t78;
t128 = t120 * t85;
t129 = -t127 + t128;
t131 = t81 + t82 + qJD(5);
t133 = t129 ^ 2;
t136 = t131 ^ 2;
t140 = t120 * t104 - t122 * t109;
t144 = t120 * t109;
t145 = t122 * t104;
t146 = t144 + t145;
t153 = t146 ^ 2;
t154 = t140 ^ 2;
t157 = pkin(7) * qJ(5) - qJ(6);
t158 = cos(t157);
t160 = sin(t157);
t162 = -t158 * t124 - t160 * t131;
t163 = t162 ^ 2;
t167 = -t160 * t124 + t158 * t131;
t170 = -pkin(7) * qJD(5) + qJD(6) - t127 + t128;
t172 = t167 ^ 2;
t175 = t170 ^ 2;
t178 = -t131 * pkin(6) + t144 + t145;
t181 = t124 * pkin(6) + t95 + t99;
t183 = t158 * t181 - t160 * t178;
t189 = -t158 * t178 - t160 * t181;
t196 = t189 ^ 2;
t197 = t183 ^ 2;
unknown(1,1) = 0;
unknown(1,2) = 0;
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = 0;
unknown(1,6) = (t1 / 0.2e1);
unknown(1,7) = 0;
unknown(1,8) = 0;
unknown(1,9) = 0;
unknown(1,10) = 0;
unknown(1,11) = (t4 * t1 / 0.2e1);
unknown(1,12) = -(t3 * t1 * t8);
unknown(1,13) = -(t3 * qJD(1) * qJD(2));
unknown(1,14) = (t12 * t1 / 0.2e1);
unknown(1,15) = (t15 * qJD(2));
unknown(1,16) = (t17 / 0.2e1);
unknown(1,17) = (t19 * t8);
unknown(1,18) = (t19 * t3);
unknown(1,19) = 0;
unknown(1,20) = (t1 * t22 / 0.2e1);
unknown(1,21) = (t32 / 0.2e1);
unknown(1,22) = (t31 * t38);
unknown(1,23) = (t31 * t40);
unknown(1,24) = (t42 / 0.2e1);
unknown(1,25) = (t38 * t40);
unknown(1,26) = (t45 / 0.2e1);
unknown(1,27) = (t47 * pkin(3) * t40 - t52 * t38);
unknown(1,28) = (-t55 * pkin(3) * t40 + t52 * t31);
unknown(1,29) = (-t47 * pkin(3) * t31 + t55 * pkin(3) * t38);
unknown(1,30) = (t65 * t17 * t67 / 0.2e1 + t69 * t17 * t67 / 0.2e1 + t72 / 0.2e1);
unknown(1,31) = (t79 / 0.2e1);
unknown(1,32) = (t78 * t83);
unknown(1,33) = (t78 * t85);
unknown(1,34) = (t87 / 0.2e1);
unknown(1,35) = (t83 * t85);
unknown(1,36) = (t90 / 0.2e1);
unknown(1,37) = (t100 * t85 - t104 * t83);
unknown(1,38) = (t104 * t78 - t109 * t85);
unknown(1,39) = (-t100 * t78 + t109 * t83);
unknown(1,40) = (t116 / 0.2e1 + t117 / 0.2e1 + t118 / 0.2e1);
unknown(1,41) = (t125 / 0.2e1);
unknown(1,42) = (t124 * t129);
unknown(1,43) = (t124 * t131);
unknown(1,44) = (t133 / 0.2e1);
unknown(1,45) = (t129 * t131);
unknown(1,46) = (t136 / 0.2e1);
unknown(1,47) = (t100 * t129 + t140 * t131);
unknown(1,48) = (-t100 * t124 - t146 * t131);
unknown(1,49) = (-t140 * t124 + t146 * t129);
unknown(1,50) = (t153 / 0.2e1 + t154 / 0.2e1 + t117 / 0.2e1);
unknown(1,51) = (t163 / 0.2e1);
unknown(1,52) = (t162 * t167);
unknown(1,53) = (t162 * t170);
unknown(1,54) = (t172 / 0.2e1);
unknown(1,55) = (t167 * t170);
unknown(1,56) = (t175 / 0.2e1);
unknown(1,57) = (-t140 * t167 + t183 * t170);
unknown(1,58) = (t140 * t162 - t189 * t170);
unknown(1,59) = (-t183 * t162 + t189 * t167);
unknown(1,60) = (t196 / 0.2e1 + t197 / 0.2e1 + t154 / 0.2e1);
T_reg = unknown;
