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
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:29
% EndTime: 2020-06-23 21:07:29
% DurationCPUTime: 0.13s
% Computational Cost: add. (83->25), mult. (194->68), div. (0->0), fcn. (140->10), ass. (0->27)
t39 = qJD(1) ^ 2;
t44 = t39 / 0.2e1;
t43 = pkin(3) * qJD(2);
t34 = sin(qJ(2));
t42 = (-t34 * pkin(3) - pkin(2)) * qJD(1);
t33 = sin(qJ(3));
t41 = t33 * t43;
t37 = cos(qJ(3));
t40 = t37 * t43;
t38 = cos(qJ(2));
t25 = (t33 * t38 + t34 * t37) * qJD(1);
t36 = cos(qJ(4));
t35 = cos(qJ(5));
t32 = sin(qJ(4));
t31 = sin(qJ(5));
t30 = qJD(2) + qJD(3);
t29 = pkin(7) * qJ(5) - qJ(6);
t27 = -t30 * pkin(5) + t41;
t26 = (t33 * t34 - t37 * t38) * qJD(1);
t24 = qJD(4) + t25;
t23 = t36 * t26 - t32 * t30;
t22 = t32 * t26 + t36 * t30 + qJD(5);
t21 = -t25 * pkin(4) + t26 * pkin(5) + t42;
t20 = t35 * t23 + t31 * t24;
t19 = -pkin(7) * qJD(5) - t31 * t23 + t35 * t24 + qJD(6);
t18 = -cos(t29) * t20 - sin(t29) * t22;
t1 = [t44, t38 ^ 2 * t44, -t38 * qJD(1) * qJD(2), qJD(2) ^ 2 / 0.2e1, t39 * pkin(2) * t34, t26 ^ 2 / 0.2e1, t26 * t25, t26 * t30, t25 * t30, t30 ^ 2 / 0.2e1, -t25 * t42 + t30 * t40, t26 * t42 - t30 * t41, t23 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -(t35 * (-t32 * t21 + t36 * t27) + t31 * (t30 * pkin(4) + t40)) * t22 + (t36 * t21 + t32 * t27) * t20, t18 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1;];
T_reg = t1;
