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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 07:59:58
% EndTime: 2020-06-20 07:59:58
% DurationCPUTime: 0.12s
% Computational Cost: add. (65->30), mult. (65->35), div. (0->0), fcn. (59->10), ass. (0->18)
t34 = cos(qJ(2));
t39 = t34 * pkin(3) + pkin(1);
t27 = qJ(2) + qJ(3);
t25 = cos(t27);
t30 = sin(qJ(4));
t38 = t25 * t30;
t33 = cos(qJ(4));
t37 = t25 * t33;
t36 = t25 * pkin(4) + t39;
t24 = sin(t27);
t35 = -t24 * pkin(5) + t36;
t32 = cos(qJ(5));
t31 = cos(qJ(6));
t29 = sin(qJ(5));
t28 = sin(qJ(6));
t22 = -t24 * t29 + t32 * t37;
t21 = -t24 * t32 - t29 * t37;
t1 = (-m(1) * rSges(1,3) - m(2) * (pkin(1) + rSges(2,3)) - m(3) * (pkin(1) + t34 * rSges(3,1) - sin(qJ(2)) * rSges(3,2)) - m(4) * (t25 * rSges(4,1) - t24 * rSges(4,2) + t39) - m(5) * ((rSges(5,1) * t33 - rSges(5,2) * t30) * t25 + (-pkin(5) - rSges(5,3)) * t24 + t36) - m(6) * (t22 * rSges(6,1) + t21 * rSges(6,2) + rSges(6,3) * t38 + t35) - m(7) * ((-t22 * t31 + t28 * t38) * rSges(7,1) + (t22 * t28 + t31 * t38) * rSges(7,2) + (pkin(6) + rSges(7,3)) * t21 + t35)) * g(3);
U = t1;
