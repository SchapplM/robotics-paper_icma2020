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
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
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
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:23
% EndTime: 2020-06-19 21:37:24
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->34), mult. (68->41), div. (0->0), fcn. (59->10), ass. (0->21)
t6 = cos(qJ(2));
t8 = sin(qJ(2));
t12 = t6 * pkin(3);
t15 = qJ(2) + qJ(3);
t16 = cos(t15);
t18 = sin(t15);
t22 = t16 * pkin(4);
t23 = t18 * pkin(5);
t24 = t22 - t23 + t12 + pkin(1);
t26 = cos(qJ(4));
t27 = t16 * t26;
t29 = sin(qJ(4));
t30 = t16 * t29;
t36 = cos(qJ(5));
t38 = sin(qJ(5));
t40 = -t18 * t38 + t27 * t36;
t44 = -t18 * t36 - t27 * t38;
t53 = pkin(7) * qJ(5) - qJ(6);
t54 = cos(t53);
t56 = sin(t53);
t67 = -g(3) * mrSges(1,3) - g(3) * (m(2) * pkin(1) + mrSges(2,3)) - g(3) * (m(3) * pkin(1) + t6 * mrSges(3,1) - t8 * mrSges(3,2)) - g(3) * (m(4) * (t12 + pkin(1)) + t16 * mrSges(4,1) - t18 * mrSges(4,2)) - g(3) * (m(5) * t24 + t27 * mrSges(5,1) - t30 * mrSges(5,2) - t18 * mrSges(5,3)) - g(3) * (m(6) * t24 + t40 * mrSges(6,1) + t44 * mrSges(6,2) + t30 * mrSges(6,3)) - g(3) * (m(7) * (t44 * pkin(6) + pkin(1) + t12 + t22 - t23) + (-t30 * t56 - t40 * t54) * mrSges(7,1) + (t30 * t54 - t40 * t56) * mrSges(7,2) + t44 * mrSges(7,3));
U = t67;
