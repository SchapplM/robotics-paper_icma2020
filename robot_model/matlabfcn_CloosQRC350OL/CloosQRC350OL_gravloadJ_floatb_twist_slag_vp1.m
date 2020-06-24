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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 21:56:05
% EndTime: 2020-06-23 21:56:05
% DurationCPUTime: 0.22s
% Computational Cost: add. (79->23), mult. (94->31), div. (0->0), fcn. (77->7), ass. (0->20)
t24 = pkin(6) + rSges(7,3);
t13 = cos(qJ(4));
t8 = qJ(2) + qJ(3);
t7 = cos(t8);
t19 = t7 * cos(qJ(5));
t6 = sin(t8);
t9 = sin(qJ(5));
t21 = t6 * t9;
t4 = t13 * t21 - t19;
t23 = t24 * t4;
t22 = pkin(4) * t6;
t11 = sin(qJ(2));
t20 = pkin(3) * t11;
t18 = -pkin(5) * t7 - t22;
t17 = -rSges(4,1) * t6 - rSges(4,2) * t7;
t16 = -t22 + (-pkin(5) - rSges(5,3)) * t7;
t15 = t18 - t20;
t14 = (-m(6) * rSges(6,2) - m(7) * t24) * g(3);
t2 = t4 * rSges(6,2);
t1 = [0, (m(3) * t11 * rSges(3,1) - m(4) * (t17 - t20) - m(5) * (t16 - t20) - m(6) * (t15 + t2) - m(7) * (t15 + t23)) * g(3), (-m(4) * t17 - m(5) * t16 - m(6) * (t18 + t2) - m(7) * (t18 + t23)) * g(3), t9 * t7 * sin(qJ(4)) * t14, (-t13 * t19 + t21) * t14, 0];
taug = t1(:);
