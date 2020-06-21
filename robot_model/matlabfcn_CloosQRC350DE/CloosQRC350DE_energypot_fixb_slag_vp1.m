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
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
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
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:23
% EndTime: 2020-06-19 21:37:24
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->37), mult. (69->42), div. (0->0), fcn. (59->10), ass. (0->20)
t7 = cos(qJ(2));
t9 = sin(qJ(2));
t14 = t7 * pkin(3);
t15 = qJ(2) + qJ(3);
t16 = cos(t15);
t18 = sin(t15);
t23 = t16 * pkin(4);
t24 = t18 * pkin(5);
t25 = cos(qJ(4));
t26 = t16 * t25;
t28 = sin(qJ(4));
t29 = t16 * t28;
t35 = cos(qJ(5));
t37 = sin(qJ(5));
t39 = -t18 * t37 + t26 * t35;
t43 = -t18 * t35 - t26 * t37;
t51 = pkin(7) * qJ(5) - qJ(6);
t52 = cos(t51);
t54 = sin(t51);
t65 = -m(1) * g(3) * rSges(1,3) - m(2) * g(3) * (pkin(1) + rSges(2,3)) - m(3) * g(3) * (t7 * rSges(3,1) - t9 * rSges(3,2) + pkin(1)) - m(4) * g(3) * (t16 * rSges(4,1) - t18 * rSges(4,2) + pkin(1) + t14) - m(5) * g(3) * (t26 * rSges(5,1) - t29 * rSges(5,2) - t18 * rSges(5,3) + pkin(1) + t14 + t23 - t24) - m(6) * g(3) * (t39 * rSges(6,1) + t43 * rSges(6,2) + t29 * rSges(6,3) + pkin(1) + t14 + t23 - t24) - m(7) * g(3) * (t43 * pkin(6) + t23 - t24 + t14 + pkin(1) + (-t29 * t54 - t39 * t52) * rSges(7,1) + (t29 * t52 - t39 * t54) * rSges(7,2) + t43 * rSges(7,3));
U = t65;
