% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(6*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = CloosQRC350OL_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:18:54
% EndTime: 2020-06-20 08:18:55
% DurationCPUTime: 0.47s
% Computational Cost: add. (674->53), mult. (1424->141), div. (0->0), fcn. (1075->10), ass. (0->45)
t47 = qJD(1) ^ 2;
t57 = t47 / 0.2e1;
t56 = cos(qJ(5));
t55 = sin(qJ(6));
t45 = cos(qJ(2));
t54 = t45 * t47;
t40 = sin(qJ(3));
t41 = sin(qJ(2));
t44 = cos(qJ(3));
t28 = (-t40 * t45 - t41 * t44) * qJD(1);
t30 = (-t40 * t41 + t44 * t45) * qJD(1);
t34 = (pkin(3) * t41 + pkin(2)) * qJD(1);
t20 = -t28 * pkin(4) + t30 * pkin(5) + t34;
t36 = qJD(2) + qJD(3);
t53 = pkin(3) * qJD(2);
t50 = t40 * t53;
t32 = -t36 * pkin(5) + t50;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t14 = -t39 * t20 + t43 * t32;
t49 = t44 * t53;
t33 = t36 * pkin(4) + t49;
t38 = sin(qJ(5));
t10 = t56 * t14 + t38 * t33;
t12 = t43 * t20 + t39 * t32;
t52 = t12 ^ 2 / 0.2e1;
t51 = qJD(1) * qJD(2);
t25 = t43 * t30 - t39 * t36;
t27 = qJD(4) + t28;
t17 = t38 * t25 - t56 * t27;
t23 = t39 * t30 + t43 * t36;
t46 = qJD(2) ^ 2;
t42 = cos(qJ(6));
t22 = qJD(5) + t23;
t19 = t56 * t25 + t38 * t27;
t15 = -qJD(6) + t17;
t9 = -t38 * t14 + t56 * t33;
t8 = t9 ^ 2 / 0.2e1;
t6 = t42 * t19 - t55 * t22;
t5 = t55 * t19 + t42 * t22;
t4 = -t22 * pkin(6) + t10;
t3 = t19 * pkin(6) + t12;
t2 = t55 * t3 - t42 * t4;
t1 = t42 * t3 + t55 * t4;
t7 = [0, 0, 0, 0, 0, t57, 0, 0, 0, 0, t45 ^ 2 * t57, -t41 * t54, t45 * t51, t41 ^ 2 * t57, -t41 * t51, t46 / 0.2e1, t47 * pkin(2) * t41, pkin(2) * t54, 0, pkin(2) ^ 2 * t57, t30 ^ 2 / 0.2e1, t30 * t28, t30 * t36, t28 ^ 2 / 0.2e1, t28 * t36, t36 ^ 2 / 0.2e1, -t34 * t28 + t36 * t49, t34 * t30 - t36 * t50, (t28 * t40 - t30 * t44) * t53, t34 ^ 2 / 0.2e1 + (t40 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1) * pkin(3) ^ 2 * t46, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t27, t23 ^ 2 / 0.2e1, -t23 * t27, t27 ^ 2 / 0.2e1, -t12 * t27 + t33 * t23, -t14 * t27 + t33 * t25, t12 * t25 - t14 * t23, t14 ^ 2 / 0.2e1 + t52 + t33 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t22, t17 ^ 2 / 0.2e1, -t17 * t22, t22 ^ 2 / 0.2e1, t12 * t17 + t9 * t22, -t10 * t22 + t12 * t19, -t10 * t17 - t9 * t19, t10 ^ 2 / 0.2e1 + t8 + t52, t6 ^ 2 / 0.2e1, -t6 * t5, t6 * t15, t5 ^ 2 / 0.2e1, -t5 * t15, t15 ^ 2 / 0.2e1, -t1 * t15 - t9 * t5, t2 * t15 - t9 * t6, t1 * t6 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8;];
T_reg = t7;
