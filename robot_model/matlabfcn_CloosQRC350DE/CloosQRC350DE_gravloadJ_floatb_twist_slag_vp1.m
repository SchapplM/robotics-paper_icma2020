% Calculate Gravitation load on the joints for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:02:09
% EndTime: 2020-06-23 21:02:09
% DurationCPUTime: 0.21s
% Computational Cost: add. (72->22), mult. (101->35), div. (0->0), fcn. (42->8), ass. (0->17)
t12 = cos(qJ(3));
t10 = cos(qJ(5));
t26 = m(6) + m(7);
t6 = (pkin(6) + rSges(7,3)) * m(7) + m(6) * rSges(6,2);
t25 = -t6 * t10 - m(4) * rSges(4,2) - t26 * pkin(5) - (pkin(5) + rSges(5,3)) * m(5);
t11 = cos(qJ(4));
t21 = m(5) + t26;
t7 = sin(qJ(5));
t24 = t6 * t7;
t4 = t21 * pkin(4) + m(4) * rSges(4,1) - t11 * t24;
t8 = sin(qJ(3));
t28 = t4 * t12 + t25 * t8;
t13 = cos(qJ(2));
t27 = t13 * (-t25 * t12 + t4 * t8);
t22 = t10 * t11;
t9 = sin(qJ(2));
t1 = [0, (t27 + (m(3) * rSges(3,1) + (m(4) + t21) * pkin(3) + t28) * t9) * g(3), g(3) * (t9 * t28 + t27), -g(3) * sin(qJ(4)) * (t13 * t12 - t8 * t9) * t24, g(3) * t6 * (-(t12 * t7 + t8 * t22) * t9 + t13 * (t12 * t22 - t7 * t8)), 0];
taug = t1(:);
