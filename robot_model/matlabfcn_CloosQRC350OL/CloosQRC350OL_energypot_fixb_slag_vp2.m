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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = CloosQRC350OL_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:55:41
% EndTime: 2020-06-23 21:55:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (41->21), mult. (38->17), div. (0->0), fcn. (29->6), ass. (0->7)
t19 = qJ(2) + qJ(3);
t17 = cos(t19);
t20 = cos(qJ(2));
t18 = t20 * pkin(3);
t21 = t17 * pkin(4) + t18;
t16 = sin(t19);
t1 = (-mrSges(2,3) - t20 * mrSges(3,1) - m(4) * t18 - t17 * mrSges(4,1) - m(7) * t21 + (m(7) * pkin(5) + mrSges(4,2) + mrSges(5,3)) * t16 + (-m(5) - m(6)) * (-t16 * pkin(5) + pkin(1) + t21) + (-m(7) * pkin(6) - mrSges(6,2) - mrSges(7,3)) * (-t17 * cos(qJ(4)) * sin(qJ(5)) - t16 * cos(qJ(5))) + (-m(2) - m(3) - m(4) - m(7)) * pkin(1)) * g(3);
U = t1;
