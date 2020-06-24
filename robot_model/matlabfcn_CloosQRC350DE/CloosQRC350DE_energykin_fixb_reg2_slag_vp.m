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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
% OptimizationMode: 2
% StartTime: 2020-06-23 21:06:21
% EndTime: 2020-06-23 21:06:21
% DurationCPUTime: 0.30s
% Computational Cost: add. (363->44), mult. (769->118), div. (0->0), fcn. (542->10), ass. (0->46)
t44 = qJD(1) ^ 2;
t52 = t44 / 0.2e1;
t37 = sin(qJ(3));
t38 = sin(qJ(2));
t41 = cos(qJ(3));
t42 = cos(qJ(2));
t23 = (t37 * t42 + t38 * t41) * qJD(1);
t50 = qJD(1) * t42;
t24 = -t37 * t38 * qJD(1) + t41 * t50;
t32 = -t38 * pkin(3) - pkin(2);
t49 = t32 * qJD(1);
t17 = -t23 * pkin(4) - t24 * pkin(5) + t49;
t34 = qJD(2) + qJD(3);
t51 = pkin(3) * qJD(2);
t47 = t37 * t51;
t27 = -t34 * pkin(5) + t47;
t36 = sin(qJ(4));
t40 = cos(qJ(4));
t13 = -t36 * t17 + t40 * t27;
t46 = t41 * t51;
t28 = t34 * pkin(4) + t46;
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t9 = t39 * t13 + t35 * t28;
t11 = t40 * t17 + t36 * t27;
t48 = t11 ^ 2 / 0.2e1;
t20 = -t40 * t24 - t36 * t34;
t22 = qJD(4) + t23;
t15 = -t35 * t20 + t39 * t22;
t19 = t36 * t24 - t40 * t34;
t43 = qJD(2) ^ 2;
t33 = pkin(7) * qJ(5) - qJ(6);
t31 = cos(t33);
t30 = sin(t33);
t18 = qJD(5) - t19;
t16 = t39 * t20 + t35 * t22;
t14 = -pkin(7) * qJD(5) + qJD(6) + t15;
t8 = -t35 * t13 + t39 * t28;
t7 = t8 ^ 2 / 0.2e1;
t6 = -t31 * t16 - t30 * t18;
t5 = -t30 * t16 + t31 * t18;
t4 = -t18 * pkin(6) + t9;
t3 = t16 * pkin(6) + t11;
t2 = -t30 * t3 - t31 * t4;
t1 = t31 * t3 - t30 * t4;
t10 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t42 ^ 2 * t52, 0, -qJD(2) * t50, t38 ^ 2 * t52, 0, t43 / 0.2e1, t44 * pkin(2) * t38, 0, 0, pkin(2) ^ 2 * t52, t24 ^ 2 / 0.2e1, -t24 * t23, -t24 * t34, t23 ^ 2 / 0.2e1, t23 * t34, t34 ^ 2 / 0.2e1, -t23 * t49 + t34 * t46, -t24 * t49 - t34 * t47, (t23 * t37 + t24 * t41) * t51, t32 ^ 2 * t52 + (t37 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1) * pkin(3) ^ 2 * t43, t20 ^ 2 / 0.2e1, 0, 0, t19 ^ 2 / 0.2e1, 0, t22 ^ 2 / 0.2e1, 0, 0, t11 * t20 + t13 * t19, t13 ^ 2 / 0.2e1 + t48 + t28 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, 0, 0, t15 ^ 2 / 0.2e1, 0, t18 ^ 2 / 0.2e1, 0, t11 * t16 - t9 * t18, 0, t9 ^ 2 / 0.2e1 + t7 + t48, t6 ^ 2 / 0.2e1, 0, 0, t5 ^ 2 / 0.2e1, 0, t14 ^ 2 / 0.2e1, 0, 0, -t1 * t6 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7;];
T_reg = t10;
