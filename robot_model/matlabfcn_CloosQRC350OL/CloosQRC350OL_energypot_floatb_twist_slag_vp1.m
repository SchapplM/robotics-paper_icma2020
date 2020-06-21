% Calculate potential energy for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = CloosQRC350OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 07:59:58
% EndTime: 2020-06-20 07:59:58
% DurationCPUTime: 0.26s
% Computational Cost: add. (72->32), mult. (65->35), div. (0->0), fcn. (59->10), ass. (0->19)
t10 = sin(qJ(4));
t7 = qJ(2) + qJ(3);
t5 = cos(t7);
t20 = t10 * t5;
t13 = cos(qJ(4));
t19 = t13 * t5;
t18 = pkin(1) + r_base(3);
t14 = cos(qJ(2));
t17 = t14 * pkin(3) + t18;
t16 = t5 * pkin(4) + t17;
t4 = sin(t7);
t15 = -t4 * pkin(5) + t16;
t12 = cos(qJ(5));
t11 = cos(qJ(6));
t9 = sin(qJ(5));
t8 = sin(qJ(6));
t2 = t12 * t19 - t4 * t9;
t1 = -t12 * t4 - t9 * t19;
t3 = (-m(1) * (r_base(3) + rSges(1,3)) - m(2) * (rSges(2,3) + t18) - m(3) * (t14 * rSges(3,1) - sin(qJ(2)) * rSges(3,2) + t18) - m(4) * (rSges(4,1) * t5 - rSges(4,2) * t4 + t17) - m(5) * ((rSges(5,1) * t13 - rSges(5,2) * t10) * t5 + (-pkin(5) - rSges(5,3)) * t4 + t16) - m(6) * (rSges(6,1) * t2 + rSges(6,2) * t1 + rSges(6,3) * t20 + t15) - m(7) * ((-t11 * t2 + t8 * t20) * rSges(7,1) + (t11 * t20 + t2 * t8) * rSges(7,2) + (pkin(6) + rSges(7,3)) * t1 + t15)) * g(3);
U = t3;
