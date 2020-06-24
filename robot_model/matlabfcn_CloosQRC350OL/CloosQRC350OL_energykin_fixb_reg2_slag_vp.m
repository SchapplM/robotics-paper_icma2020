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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 22:02:20
% EndTime: 2020-06-23 22:02:20
% DurationCPUTime: 0.35s
% Computational Cost: add. (345->39), mult. (772->116), div. (0->0), fcn. (542->10), ass. (0->43)
t44 = qJD(1) ^ 2;
t50 = t44 / 0.2e1;
t36 = sin(qJ(3));
t37 = sin(qJ(2));
t41 = cos(qJ(3));
t42 = cos(qJ(2));
t23 = (-t36 * t42 - t37 * t41) * qJD(1);
t25 = (-t36 * t37 + t41 * t42) * qJD(1);
t29 = (pkin(3) * t37 + pkin(2)) * qJD(1);
t17 = -t23 * pkin(4) + t25 * pkin(5) + t29;
t31 = qJD(2) + qJD(3);
t49 = pkin(3) * qJD(2);
t47 = t36 * t49;
t27 = -t31 * pkin(5) + t47;
t35 = sin(qJ(4));
t40 = cos(qJ(4));
t13 = -t35 * t17 + t40 * t27;
t46 = t41 * t49;
t28 = t31 * pkin(4) + t46;
t34 = sin(qJ(5));
t39 = cos(qJ(5));
t9 = t39 * t13 + t34 * t28;
t11 = t40 * t17 + t35 * t27;
t48 = t11 ^ 2 / 0.2e1;
t20 = t40 * t25 - t35 * t31;
t22 = qJD(4) + t23;
t15 = -t34 * t20 + t39 * t22;
t19 = -t35 * t25 - t40 * t31;
t43 = qJD(2) ^ 2;
t38 = cos(qJ(6));
t33 = sin(qJ(6));
t18 = qJD(5) - t19;
t16 = t39 * t20 + t34 * t22;
t14 = qJD(6) + t15;
t8 = -t34 * t13 + t39 * t28;
t7 = t8 ^ 2 / 0.2e1;
t6 = -t38 * t16 + t33 * t18;
t5 = t33 * t16 + t38 * t18;
t4 = -t18 * pkin(6) + t9;
t3 = t16 * pkin(6) + t11;
t2 = t33 * t3 - t38 * t4;
t1 = t38 * t3 + t33 * t4;
t10 = [0, 0, 0, 0, 0, t50, 0, 0, 0, 0, t42 ^ 2 * t50, 0, t42 * qJD(1) * qJD(2), t37 ^ 2 * t50, 0, t43 / 0.2e1, t44 * pkin(2) * t37, 0, 0, pkin(2) ^ 2 * t50, t25 ^ 2 / 0.2e1, t25 * t23, t25 * t31, t23 ^ 2 / 0.2e1, t23 * t31, t31 ^ 2 / 0.2e1, -t29 * t23 + t31 * t46, t29 * t25 - t31 * t47, (t23 * t36 - t25 * t41) * t49, t29 ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1) * pkin(3) ^ 2 * t43, t20 ^ 2 / 0.2e1, 0, 0, t19 ^ 2 / 0.2e1, 0, t22 ^ 2 / 0.2e1, 0, 0, t11 * t20 + t13 * t19, t13 ^ 2 / 0.2e1 + t48 + t28 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, 0, 0, t15 ^ 2 / 0.2e1, 0, t18 ^ 2 / 0.2e1, 0, t11 * t16 - t9 * t18, 0, t9 ^ 2 / 0.2e1 + t7 + t48, t6 ^ 2 / 0.2e1, 0, 0, t5 ^ 2 / 0.2e1, 0, t14 ^ 2 / 0.2e1, 0, 0, -t1 * t6 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7;];
T_reg = t10;
