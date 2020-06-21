% Calculate Gravitation load on the joints for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350OL_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:00:26
% EndTime: 2020-06-20 08:00:27
% DurationCPUTime: 0.44s
% Computational Cost: add. (164->52), mult. (196->82), div. (0->0), fcn. (194->10), ass. (0->27)
t35 = pkin(6) + rSges(7,3);
t16 = sin(qJ(6));
t20 = cos(qJ(6));
t34 = -rSges(7,1) * t20 + rSges(7,2) * t16;
t15 = qJ(2) + qJ(3);
t13 = sin(t15);
t14 = cos(t15);
t21 = cos(qJ(5));
t17 = sin(qJ(5));
t22 = cos(qJ(4));
t28 = t17 * t22;
t7 = t13 * t28 - t14 * t21;
t26 = t21 * t22;
t8 = -t13 * t26 - t14 * t17;
t33 = t8 * rSges(6,1) + t7 * rSges(6,2);
t19 = sin(qJ(2));
t32 = t19 * pkin(3);
t18 = sin(qJ(4));
t29 = t16 * t18;
t27 = t18 * t20;
t25 = (-t13 * t27 + t8 * t16) * rSges(7,2) + (-t13 * t29 - t8 * t20) * rSges(7,1) + t35 * t7;
t24 = -t14 * pkin(5) - t32;
t23 = (m(4) * rSges(4,1) - m(5) * (-rSges(5,1) * t22 - pkin(4)) - m(6) * (-rSges(6,3) * t18 - pkin(4)) + m(7) * pkin(4)) * t13;
t11 = t13 * t18 * rSges(5,2);
t10 = -t13 * t17 + t14 * t26;
t9 = -t13 * t21 - t14 * t28;
t1 = [0, (-m(3) * (-t19 * rSges(3,1) - cos(qJ(2)) * rSges(3,2)) - m(4) * (-t14 * rSges(4,2) - t32) - m(5) * (-t14 * rSges(5,3) + t11 + t24) - m(6) * (t24 + t33) - m(7) * (t24 + t25) + t23) * g(3), (-m(5) * t11 - m(6) * t33 - m(7) * t25 + (m(4) * rSges(4,2) + m(5) * rSges(5,3) + (m(5) + m(6) + m(7)) * pkin(5)) * t14 + t23) * g(3), ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * (t16 * rSges(7,1) + t20 * rSges(7,2))) * t22 + (m(5) * rSges(5,1) - m(6) * (-rSges(6,1) * t21 + rSges(6,2) * t17) - m(7) * (t35 * t17 - t34 * t21)) * t18) * g(3) * t14, (-m(6) * (t9 * rSges(6,1) - t10 * rSges(6,2)) - m(7) * (-t35 * t10 + t34 * t9)) * g(3), -m(7) * g(3) * ((t10 * t16 + t14 * t27) * rSges(7,1) + (t10 * t20 - t14 * t29) * rSges(7,2))];
taug = t1(:);
