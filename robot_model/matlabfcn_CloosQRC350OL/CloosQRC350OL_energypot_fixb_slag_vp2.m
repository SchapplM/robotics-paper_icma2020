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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 07:59:58
% EndTime: 2020-06-20 07:59:58
% DurationCPUTime: 0.13s
% Computational Cost: add. (65->32), mult. (64->27), div. (0->0), fcn. (59->10), ass. (0->13)
t28 = qJ(2) + qJ(3);
t26 = cos(t28);
t34 = cos(qJ(4));
t37 = t26 * t34;
t35 = cos(qJ(2));
t27 = t35 * pkin(3);
t36 = t26 * pkin(4) + t27;
t33 = cos(qJ(5));
t32 = cos(qJ(6));
t30 = sin(qJ(5));
t29 = sin(qJ(6));
t25 = sin(t28);
t1 = (-mrSges(1,3) - mrSges(2,3) - t35 * mrSges(3,1) + sin(qJ(2)) * mrSges(3,2) - m(4) * t27 - m(7) * t36 + (m(7) * pkin(5) + mrSges(4,2) + mrSges(5,3)) * t25 + (-m(5) - m(6)) * (-t25 * pkin(5) + pkin(1) + t36) + (-m(7) * pkin(6) - mrSges(6,2) - mrSges(7,3)) * (-t25 * t33 - t30 * t37) + (-m(2) - m(3) - m(4) - m(7)) * pkin(1) + (t32 * mrSges(7,1) - t29 * mrSges(7,2) - mrSges(6,1)) * (-t25 * t30 + t33 * t37) + (-mrSges(5,1) * t34 - mrSges(4,1) + (-mrSges(7,1) * t29 - mrSges(7,2) * t32 + mrSges(5,2) - mrSges(6,3)) * sin(qJ(4))) * t26) * g(3);
U = t1;
