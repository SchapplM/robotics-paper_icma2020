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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = CloosQRC350DE_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 20:58:57
% EndTime: 2020-06-23 20:58:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (42->19), mult. (38->15), div. (0->0), fcn. (13->7), ass. (0->8)
t20 = m(5) + m(6) + m(7);
t22 = m(4) + t20;
t14 = pkin(6) * m(7) + mrSges(6,2) + mrSges(7,3);
t21 = t20 * pkin(5) + t14 * cos(qJ(5)) + mrSges(4,2) + mrSges(5,3);
t18 = cos(qJ(3));
t15 = sin(qJ(3));
t13 = -sin(qJ(5)) * t14 * cos(qJ(4)) + mrSges(4,1) + t20 * pkin(4);
t1 = -((t22 * pkin(3) + t13 * t18 - t21 * t15 + mrSges(3,1)) * cos(qJ(2)) + mrSges(2,3) + (-t13 * t15 - t21 * t18) * sin(qJ(2)) + (m(2) + m(3) + t22) * pkin(1)) * g(3);
U = t1;
