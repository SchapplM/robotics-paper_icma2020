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

function U = CloosQRC350OL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:55:40
% EndTime: 2020-06-23 21:55:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (48->21), mult. (40->18), div. (0->0), fcn. (29->6), ass. (0->8)
t12 = pkin(1) + r_base(3);
t9 = cos(qJ(2));
t11 = t9 * pkin(3) + t12;
t8 = qJ(2) + qJ(3);
t5 = cos(t8);
t10 = t5 * pkin(4) + t11;
t4 = sin(t8);
t1 = (-m(1) * r_base(3) - mrSges(2,3) - t9 * mrSges(3,1) - m(4) * t11 - t5 * mrSges(4,1) - m(7) * t10 + (-m(2) - m(3)) * t12 + (m(7) * pkin(5) + mrSges(4,2) + mrSges(5,3)) * t4 + (-m(7) * pkin(6) - mrSges(6,2) - mrSges(7,3)) * (-t5 * cos(qJ(4)) * sin(qJ(5)) - t4 * cos(qJ(5))) + (-m(5) - m(6)) * (-pkin(5) * t4 + t10)) * g(3);
U = t1;
