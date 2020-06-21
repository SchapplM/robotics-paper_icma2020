% Calculate kinetic energy for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
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

function T = CloosQRC350OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 07:59:58
% EndTime: 2020-06-20 08:00:01
% DurationCPUTime: 3.34s
% Computational Cost: add. (2745->146), mult. (3585->218), div. (0->0), fcn. (2936->12), ass. (0->57)
t60 = sin(qJ(1));
t66 = cos(qJ(1));
t48 = t60 * V_base(5) + t66 * V_base(4);
t54 = V_base(6) + qJD(1);
t59 = sin(qJ(2));
t65 = cos(qJ(2));
t39 = t48 * t59 + t54 * t65;
t40 = t48 * t65 - t54 * t59;
t58 = sin(qJ(3));
t64 = cos(qJ(3));
t31 = -t39 * t58 + t64 * t40;
t32 = t39 * t64 + t40 * t58;
t50 = V_base(5) * pkin(1) + V_base(1);
t51 = -V_base(4) * pkin(1) + V_base(2);
t42 = -t50 * t60 + t66 * t51;
t38 = t54 * pkin(2) + t42;
t33 = -pkin(3) * t40 + t38;
t14 = -pkin(4) * t31 + pkin(5) * t32 + t33;
t43 = t50 * t66 + t51 * t60;
t47 = -t60 * V_base(4) + t66 * V_base(5);
t44 = -pkin(2) * t47 + V_base(3);
t36 = t65 * t43 - t44 * t59;
t46 = qJD(2) + t47;
t29 = pkin(3) * t46 + t36;
t35 = t43 * t59 + t44 * t65;
t22 = t58 * t29 + t64 * t35;
t45 = qJD(3) + t46;
t20 = -pkin(5) * t45 + t22;
t57 = sin(qJ(4));
t63 = cos(qJ(4));
t9 = t14 * t63 + t20 * t57;
t68 = t9 ^ 2;
t11 = -t14 * t57 + t20 * t63;
t21 = t64 * t29 - t35 * t58;
t19 = pkin(4) * t45 + t21;
t56 = sin(qJ(5));
t62 = cos(qJ(5));
t7 = t62 * t11 + t56 * t19;
t25 = t32 * t63 - t45 * t57;
t30 = qJD(4) + t31;
t16 = -t25 * t56 + t62 * t30;
t24 = -t32 * t57 - t45 * t63;
t67 = V_base(3) ^ 2;
t61 = cos(qJ(6));
t55 = sin(qJ(6));
t23 = qJD(5) - t24;
t17 = t25 * t62 + t30 * t56;
t15 = qJD(6) + t16;
t13 = -t17 * t61 + t23 * t55;
t12 = t17 * t55 + t23 * t61;
t6 = -t11 * t56 + t19 * t62;
t5 = t6 ^ 2;
t4 = pkin(6) * t17 + t9;
t3 = -pkin(6) * t23 + t7;
t2 = -t3 * t61 + t4 * t55;
t1 = t3 * t55 + t4 * t61;
t8 = (-t9 * mrSges(6,1) + t7 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t23 + Ifges(6,2) * t16 / 0.2e1) * t16 + (-t6 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t13 + Ifges(7,6) * t15 + Ifges(7,2) * t12 / 0.2e1) * t12 + (t6 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t15 + Ifges(7,1) * t13 / 0.2e1) * t13 + (t42 * mrSges(2,1) - t43 * mrSges(2,2) + Ifges(2,3) * t54 / 0.2e1) * t54 + (-V_base(3) * mrSges(2,1) + t43 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t54 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t38 * mrSges(3,1) + t35 * mrSges(3,3) + Ifges(3,6) * t46 + Ifges(3,2) * t40 / 0.2e1) * t40 + (-t33 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,4) * t32 + Ifges(4,6) * t45 + Ifges(4,2) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t6 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,3) * t23 / 0.2e1) * t23 + (t36 * mrSges(3,1) - t35 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t67) / 0.2e1 + m(2) * (t42 ^ 2 + t43 ^ 2 + t67) / 0.2e1 + (-t9 * mrSges(5,1) - t11 * mrSges(5,2) + Ifges(5,3) * t30 / 0.2e1) * t30 + (t19 * mrSges(5,2) + t9 * mrSges(5,3) + Ifges(5,5) * t30 + Ifges(5,1) * t25 / 0.2e1) * t25 + (t33 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,5) * t45 + Ifges(4,1) * t32 / 0.2e1) * t32 + m(4) * (t21 ^ 2 + t22 ^ 2 + t33 ^ 2) / 0.2e1 + m(3) * (t35 ^ 2 + t36 ^ 2 + t38 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5) / 0.2e1 + (-t19 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t30 + Ifges(5,2) * t24 / 0.2e1) * t24 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t15 / 0.2e1) * t15 + (t38 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,4) * t40 + Ifges(3,5) * t46 + Ifges(3,1) * t39 / 0.2e1) * t39 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + m(5) * (t11 ^ 2 + t19 ^ 2 + t68) / 0.2e1 + m(6) * (t7 ^ 2 + t5 + t68) / 0.2e1 + (t9 * mrSges(6,2) - t6 * mrSges(6,3) + Ifges(6,5) * t23 + Ifges(6,1) * t17 / 0.2e1) * t17 + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * t45 / 0.2e1) * t45 + (V_base(3) * mrSges(2,2) - t42 * mrSges(2,3) + Ifges(2,5) * t54 + Ifges(2,1) * t48 / 0.2e1) * t48;
T = t8;
