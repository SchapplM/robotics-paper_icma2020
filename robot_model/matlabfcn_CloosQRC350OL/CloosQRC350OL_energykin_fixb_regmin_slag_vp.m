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
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 22:03:17
% EndTime: 2020-06-23 22:03:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (79->24), mult. (195->66), div. (0->0), fcn. (140->10), ass. (0->26)
t34 = sin(qJ(3));
t35 = sin(qJ(2));
t38 = cos(qJ(3));
t39 = cos(qJ(2));
t25 = (t34 * t39 + t35 * t38) * qJD(1);
t40 = qJD(1) ^ 2;
t45 = t40 / 0.2e1;
t28 = (pkin(3) * t35 + pkin(2)) * qJD(1);
t44 = pkin(3) * qJD(2);
t43 = t34 * t44;
t42 = t38 * t44;
t37 = cos(qJ(4));
t36 = cos(qJ(5));
t33 = sin(qJ(4));
t32 = sin(qJ(5));
t30 = qJD(2) + qJD(3);
t27 = -t30 * pkin(5) + t43;
t26 = (-t34 * t35 + t38 * t39) * qJD(1);
t24 = qJD(4) - t25;
t23 = t37 * t26 - t33 * t30;
t22 = t33 * t26 + t37 * t30 + qJD(5);
t21 = t25 * pkin(4) + t26 * pkin(5) + t28;
t20 = t36 * t23 + t32 * t24;
t19 = -t32 * t23 + t36 * t24 + qJD(6);
t18 = -cos(qJ(6)) * t20 + sin(qJ(6)) * t22;
t1 = [t45, t39 ^ 2 * t45, t39 * qJD(1) * qJD(2), qJD(2) ^ 2 / 0.2e1, t40 * pkin(2) * t35, t26 ^ 2 / 0.2e1, -t26 * t25, t26 * t30, -t25 * t30, t30 ^ 2 / 0.2e1, t28 * t25 + t30 * t42, t28 * t26 - t30 * t43, t23 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -(t36 * (-t33 * t21 + t37 * t27) + t32 * (t30 * pkin(4) + t42)) * t22 + (t37 * t21 + t33 * t27) * t20, t18 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1;];
T_reg = t1;
