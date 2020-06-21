% Calculate Gravitation load on the joints for
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
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:31
% EndTime: 2020-06-19 21:37:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (198->72), mult. (234->94), div. (0->0), fcn. (212->10), ass. (0->46)
unknown=NaN(6,1);
t2 = sin(qJ(2));
t4 = cos(qJ(2));
t8 = m(4) * g(3);
t9 = t2 * pkin(3);
t10 = qJ(2) + qJ(3);
t11 = sin(t10);
t12 = t11 * rSges(4,1);
t13 = cos(t10);
t14 = t13 * rSges(4,2);
t17 = m(5) * g(3);
t18 = t11 * pkin(4);
t19 = t13 * pkin(5);
t20 = cos(qJ(4));
t21 = t11 * t20;
t22 = t21 * rSges(5,1);
t23 = sin(qJ(4));
t24 = t11 * t23;
t25 = t24 * rSges(5,2);
t26 = t13 * rSges(5,3);
t29 = m(6) * g(3);
t30 = cos(qJ(5));
t32 = sin(qJ(5));
t34 = -t13 * t32 - t21 * t30;
t35 = t34 * rSges(6,1);
t38 = -t13 * t30 + t21 * t32;
t39 = t38 * rSges(6,2);
t40 = t24 * rSges(6,3);
t43 = m(7) * g(3);
t44 = t38 * pkin(6);
t46 = pkin(7) * qJ(5) - qJ(6);
t47 = cos(t46);
t49 = sin(t46);
t52 = (t24 * t49 - t34 * t47) * rSges(7,1);
t56 = (-t24 * t47 - t34 * t49) * rSges(7,2);
t57 = t38 * rSges(7,3);
t70 = t13 * t23;
t72 = t13 * t20;
t102 = -t11 * t30 - t72 * t32;
t106 = t11 * t32 - t72 * t30;
t112 = -t106 * pkin(7);
unknown(1) = 0;
unknown(2) = (-m(3) * g(3) * (-t2 * rSges(3,1) - t4 * rSges(3,2)) - t8 * (-t9 - t12 - t14) - t17 * (-t18 - t19 - t9 - t22 + t25 - t26) - t29 * (-t18 - t19 - t9 + t35 + t39 - t40) - t43 * (t44 - t18 - t19 - t9 + t52 + t56 + t57));
unknown(3) = (-t8 * (-t12 - t14) - t17 * (-t18 - t19 - t22 + t25 - t26) - t29 * (-t18 - t19 + t35 + t39 - t40) - t43 * (t44 - t18 - t19 + t52 + t56 + t57));
unknown(4) = (-t17 * (-t70 * rSges(5,1) - t72 * rSges(5,2)) - t29 * (-t70 * t30 * rSges(6,1) + t70 * t32 * rSges(6,2) + t72 * rSges(6,3)) - t43 * (t70 * t32 * pkin(6) + (t70 * t30 * t47 - t72 * t49) * rSges(7,1) + (t70 * t30 * t49 + t72 * t47) * rSges(7,2) + t70 * t32 * rSges(7,3)));
unknown(5) = (-t29 * (t102 * rSges(6,1) + t106 * rSges(6,2)) - t43 * (t106 * pkin(6) + (-t70 * pkin(7) * t47 - t102 * t47 + t112 * t49) * rSges(7,1) + (-t70 * pkin(7) * t49 - t102 * t49 - t112 * t47) * rSges(7,2) + t106 * rSges(7,3)));
unknown(6) = -(t43 * ((t106 * t49 + t70 * t47) * rSges(7,1) + (-t106 * t47 + t70 * t49) * rSges(7,2)));
taug = unknown(:);
