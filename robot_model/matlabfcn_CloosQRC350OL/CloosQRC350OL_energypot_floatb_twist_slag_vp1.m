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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 21:55:40
% EndTime: 2020-06-23 21:55:41
% DurationCPUTime: 0.15s
% Computational Cost: add. (48->21), mult. (40->20), div. (0->0), fcn. (29->6), ass. (0->10)
t11 = pkin(1) + r_base(3);
t7 = cos(qJ(2));
t10 = t7 * pkin(3) + t11;
t6 = qJ(2) + qJ(3);
t4 = cos(t6);
t9 = t4 * pkin(4) + t10;
t3 = sin(t6);
t8 = -pkin(5) * t3 + t9;
t1 = -t4 * cos(qJ(4)) * sin(qJ(5)) - t3 * cos(qJ(5));
t2 = (-m(1) * r_base(3) - m(2) * (rSges(2,3) + t11) - m(3) * (rSges(3,1) * t7 + t11) - m(4) * (rSges(4,1) * t4 - rSges(4,2) * t3 + t10) - m(5) * ((-pkin(5) - rSges(5,3)) * t3 + t9) - m(6) * (rSges(6,2) * t1 + t8) - m(7) * ((pkin(6) + rSges(7,3)) * t1 + t8)) * g(3);
U = t2;
