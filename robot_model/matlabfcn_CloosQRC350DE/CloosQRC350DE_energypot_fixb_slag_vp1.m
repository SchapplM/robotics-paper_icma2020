% Calculate potential energy for
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = CloosQRC350DE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:01:18
% EndTime: 2020-06-23 21:01:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (42->21), mult. (49->21), div. (0->0), fcn. (13->7), ass. (0->9)
t28 = -m(7) - m(6);
t30 = m(5) - t28;
t31 = m(4) + t30;
t19 = (pkin(6) + rSges(7,3)) * m(7) + m(6) * rSges(6,2);
t29 = t28 * pkin(5) - (pkin(5) + rSges(5,3)) * m(5) - m(4) * rSges(4,2) - t19 * cos(qJ(5));
t23 = cos(qJ(3));
t20 = sin(qJ(3));
t18 = -t19 * sin(qJ(5)) * cos(qJ(4)) + m(4) * rSges(4,1) + t30 * pkin(4);
t1 = -g(3) * ((t31 * pkin(3) + m(3) * rSges(3,1) + t18 * t23 + t29 * t20) * cos(qJ(2)) + m(2) * rSges(2,3) + (-t18 * t20 + t29 * t23) * sin(qJ(2)) + (m(2) + m(3) + t31) * pkin(1));
U = t1;
