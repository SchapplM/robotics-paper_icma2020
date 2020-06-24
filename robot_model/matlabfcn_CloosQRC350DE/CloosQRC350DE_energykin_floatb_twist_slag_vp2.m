% Calculate kinetic energy for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = CloosQRC350DE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 20:59:50
% EndTime: 2020-06-23 20:59:51
% DurationCPUTime: 1.46s
% Computational Cost: add. (1689->126), mult. (2404->203), div. (0->0), fcn. (2101->12), ass. (0->63)
t77 = sin(qJ(1));
t82 = cos(qJ(1));
t51 = -t77 * V_base(5) + t82 * V_base(4);
t90 = -t77 * V_base(2) + t82 * V_base(1);
t72 = V_base(6) - qJD(1);
t76 = sin(qJ(2));
t81 = cos(qJ(2));
t37 = t76 * t51 + t81 * t72;
t38 = t81 * t51 - t76 * t72;
t75 = sin(qJ(3));
t80 = cos(qJ(3));
t28 = -t75 * t37 + t80 * t38;
t29 = t80 * t37 + t75 * t38;
t30 = V_base(1) * t77 + V_base(2) * t82 - t51 * (t81 * pkin(3) + pkin(1)) + t72 * (t76 * pkin(3) + pkin(2));
t14 = -t28 * pkin(4) + t29 * pkin(5) + t30;
t48 = t77 * V_base(4) + t82 * V_base(5);
t49 = t81 * t75 + t80 * t76;
t50 = -t76 * t75 + t81 * t80;
t55 = t81 * pkin(1) + t76 * pkin(2) + pkin(3);
t56 = -t76 * pkin(1) + t81 * pkin(2);
t88 = pkin(3) * qJD(2);
t21 = -t48 * (-t55 * t75 + t80 * t56) + t50 * V_base(3) + t75 * t88 + t90 * t49;
t47 = qJD(2) + t48;
t44 = qJD(3) + t47;
t19 = -t44 * pkin(5) + t21;
t74 = sin(qJ(4));
t79 = cos(qJ(4));
t9 = t79 * t14 + t74 * t19;
t89 = t9 ^ 2;
t11 = -t74 * t14 + t79 * t19;
t22 = t48 * (t55 * t80 + t75 * t56) - t49 * V_base(3) + t80 * t88 + t90 * t50;
t20 = t44 * pkin(4) + t22;
t73 = sin(qJ(5));
t78 = cos(qJ(5));
t7 = t78 * t11 + t73 * t20;
t59 = V_base(5) * pkin(1) + V_base(1);
t60 = V_base(4) * pkin(1) - V_base(2);
t39 = t77 * t59 - t82 * t60;
t25 = t79 * t29 - t74 * t44;
t27 = qJD(4) + t28;
t16 = -t73 * t25 + t78 * t27;
t24 = -t74 * t29 - t79 * t44;
t83 = V_base(3) ^ 2;
t71 = pkin(7) * qJ(5) - qJ(6);
t62 = cos(t71);
t61 = sin(t71);
t43 = -t48 * pkin(2) + V_base(3);
t40 = t82 * t59 + t77 * t60;
t35 = t72 * pkin(2) + t39;
t32 = t81 * t40 - t76 * t43;
t31 = t76 * t40 + t81 * t43;
t23 = qJD(5) - t24;
t17 = t78 * t25 + t73 * t27;
t15 = -pkin(7) * qJD(5) + qJD(6) + t16;
t13 = -t62 * t17 - t61 * t23;
t12 = -t61 * t17 + t62 * t23;
t6 = -t73 * t11 + t78 * t20;
t5 = t6 ^ 2;
t4 = t17 * pkin(6) + t9;
t3 = -t23 * pkin(6) + t7;
t2 = -t62 * t3 - t61 * t4;
t1 = -t61 * t3 + t62 * t4;
t8 = (t22 * mrSges(4,1) - t21 * mrSges(4,2) + Ifges(4,3) * t44 / 0.2e1) * t44 + (-t35 * mrSges(3,1) + t31 * mrSges(3,3) + Ifges(3,2) * t38 / 0.2e1) * t38 + (-t32 * mrSges(3,3) + Ifges(3,5) * t47 + Ifges(3,1) * t37 / 0.2e1) * t37 + (V_base(3) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,5) * t72 + Ifges(2,1) * t51 / 0.2e1) * t51 + (-t30 * mrSges(4,1) + t21 * mrSges(4,3) + Ifges(4,4) * t29 + Ifges(4,6) * t44 + Ifges(4,2) * t28 / 0.2e1) * t28 + t27 ^ 2 * Ifges(5,3) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t30 ^ 2) / 0.2e1 + m(3) * (t31 ^ 2 + t32 ^ 2 + t35 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(2,1) + t40 * mrSges(2,3) + Ifges(2,4) * t51 + Ifges(2,6) * t72 + Ifges(2,2) * t48 / 0.2e1) * t48 + (t32 * mrSges(3,1) + Ifges(3,3) * t47 / 0.2e1) * t47 + (t30 * mrSges(4,2) - t22 * mrSges(4,3) + Ifges(4,5) * t44 + Ifges(4,1) * t29 / 0.2e1) * t29 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5) / 0.2e1 + t12 ^ 2 * Ifges(7,2) / 0.2e1 + t13 ^ 2 * Ifges(7,1) / 0.2e1 + (t39 * mrSges(2,1) - t40 * mrSges(2,2) + Ifges(2,3) * t72 / 0.2e1) * t72 + (t9 * t17 - t7 * t23) * mrSges(6,2) + t15 ^ 2 * Ifges(7,3) / 0.2e1 + t16 ^ 2 * Ifges(6,2) / 0.2e1 + t17 ^ 2 * Ifges(6,1) / 0.2e1 + t23 ^ 2 * Ifges(6,3) / 0.2e1 + t24 ^ 2 * Ifges(5,2) / 0.2e1 + t25 ^ 2 * Ifges(5,1) / 0.2e1 + (t11 * t24 + t9 * t25) * mrSges(5,3) + m(2) * (t39 ^ 2 + t40 ^ 2 + t83) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t83) / 0.2e1 + m(6) * (t7 ^ 2 + t5 + t89) / 0.2e1 + m(5) * (t11 ^ 2 + t20 ^ 2 + t89) / 0.2e1 + (-t1 * t13 + t2 * t12) * mrSges(7,3);
T = t8;
