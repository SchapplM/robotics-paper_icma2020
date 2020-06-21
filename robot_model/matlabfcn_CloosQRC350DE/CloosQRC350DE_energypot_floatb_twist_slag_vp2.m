% Calculate potential energy for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = CloosQRC350DE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:23
% EndTime: 2020-06-19 21:37:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (76->39), mult. (69->42), div. (0->0), fcn. (59->10), ass. (0->22)
t4 = pkin(1) + r_base(3);
t9 = cos(qJ(2));
t11 = sin(qJ(2));
t15 = t9 * pkin(3);
t18 = qJ(2) + qJ(3);
t19 = cos(t18);
t21 = sin(t18);
t25 = t19 * pkin(4);
t26 = t21 * pkin(5);
t27 = t25 - t26 + t15 + pkin(1) + r_base(3);
t29 = cos(qJ(4));
t30 = t19 * t29;
t32 = sin(qJ(4));
t33 = t19 * t32;
t39 = cos(qJ(5));
t41 = sin(qJ(5));
t43 = -t21 * t41 + t30 * t39;
t47 = -t21 * t39 - t30 * t41;
t56 = pkin(7) * qJ(5) - qJ(6);
t57 = cos(t56);
t59 = sin(t56);
t70 = -g(3) * (m(1) * r_base(3) + mrSges(1,3)) - g(3) * (m(2) * t4 + mrSges(2,3)) - g(3) * (m(3) * t4 + t9 * mrSges(3,1) - t11 * mrSges(3,2)) - g(3) * (m(4) * (t15 + pkin(1) + r_base(3)) + t19 * mrSges(4,1) - t21 * mrSges(4,2)) - g(3) * (m(5) * t27 + t30 * mrSges(5,1) - t33 * mrSges(5,2) - t21 * mrSges(5,3)) - g(3) * (m(6) * t27 + t43 * mrSges(6,1) + t47 * mrSges(6,2) + t33 * mrSges(6,3)) - g(3) * (m(7) * (t47 * pkin(6) + pkin(1) + t15 + t25 - t26 + r_base(3)) + (-t33 * t59 - t43 * t57) * mrSges(7,1) + (t33 * t57 - t43 * t59) * mrSges(7,2) + t47 * mrSges(7,3));
U = t70;
