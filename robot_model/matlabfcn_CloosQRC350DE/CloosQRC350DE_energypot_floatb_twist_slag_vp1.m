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

function U = CloosQRC350DE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:23
% EndTime: 2020-06-19 21:37:23
% DurationCPUTime: 0.14s
% Computational Cost: add. (76->44), mult. (69->42), div. (0->0), fcn. (59->10), ass. (0->20)
t8 = cos(qJ(2));
t10 = sin(qJ(2));
t15 = t8 * pkin(3);
t16 = qJ(2) + qJ(3);
t17 = cos(t16);
t19 = sin(t16);
t24 = t17 * pkin(4);
t25 = t19 * pkin(5);
t26 = cos(qJ(4));
t27 = t17 * t26;
t29 = sin(qJ(4));
t30 = t17 * t29;
t36 = cos(qJ(5));
t38 = sin(qJ(5));
t40 = -t19 * t38 + t27 * t36;
t44 = -t19 * t36 - t27 * t38;
t52 = pkin(7) * qJ(5) - qJ(6);
t53 = cos(t52);
t55 = sin(t52);
t66 = -m(1) * g(3) * (r_base(3) + rSges(1,3)) - m(2) * g(3) * (pkin(1) + r_base(3) + rSges(2,3)) - m(3) * g(3) * (rSges(3,1) * t8 - rSges(3,2) * t10 + pkin(1) + r_base(3)) - m(4) * g(3) * (rSges(4,1) * t17 - rSges(4,2) * t19 + pkin(1) + t15 + r_base(3)) - m(5) * g(3) * (rSges(5,1) * t27 - rSges(5,2) * t30 - rSges(5,3) * t19 + pkin(1) + t15 + t24 - t25 + r_base(3)) - m(6) * g(3) * (rSges(6,1) * t40 + rSges(6,2) * t44 + rSges(6,3) * t30 + pkin(1) + t15 + t24 - t25 + r_base(3)) - m(7) * g(3) * (t44 * pkin(6) + t24 - t25 + t15 + pkin(1) + r_base(3) + (-t30 * t55 - t40 * t53) * rSges(7,1) + (t30 * t53 - t40 * t55) * rSges(7,2) + t44 * rSges(7,3));
U = t66;
