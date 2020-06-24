% Calculate potential energy for
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = CloosQRC350OL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:55:41
% EndTime: 2020-06-23 21:55:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (41->19), mult. (38->19), div. (0->0), fcn. (29->6), ass. (0->9)
t18 = cos(qJ(2));
t21 = t18 * pkin(3) + pkin(1);
t17 = qJ(2) + qJ(3);
t15 = cos(t17);
t20 = t15 * pkin(4) + t21;
t14 = sin(t17);
t19 = -t14 * pkin(5) + t20;
t12 = -t15 * cos(qJ(4)) * sin(qJ(5)) - t14 * cos(qJ(5));
t1 = (-m(2) * (pkin(1) + rSges(2,3)) - m(3) * (t18 * rSges(3,1) + pkin(1)) - m(4) * (t15 * rSges(4,1) - t14 * rSges(4,2) + t21) - m(5) * ((-pkin(5) - rSges(5,3)) * t14 + t20) - m(6) * (t12 * rSges(6,2) + t19) - m(7) * ((pkin(6) + rSges(7,3)) * t12 + t19)) * g(3);
U = t1;
