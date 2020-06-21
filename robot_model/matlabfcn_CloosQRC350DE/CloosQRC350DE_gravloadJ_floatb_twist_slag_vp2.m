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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350DE_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:37:31
% EndTime: 2020-06-19 21:37:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (198->69), mult. (228->100), div. (0->0), fcn. (212->10), ass. (0->44)
unknown=NaN(6,1);
t1 = cos(qJ(2));
t3 = sin(qJ(2));
t9 = qJ(2) + qJ(3);
t10 = sin(t9);
t11 = t10 * mrSges(4,1);
t12 = cos(t9);
t13 = t12 * mrSges(4,2);
t16 = t10 * pkin(4);
t17 = t12 * pkin(5);
t18 = t3 * pkin(3);
t19 = -t16 - t17 - t18;
t21 = cos(qJ(4));
t22 = t10 * t21;
t23 = t22 * mrSges(5,1);
t24 = sin(qJ(4));
t25 = t10 * t24;
t26 = t25 * mrSges(5,2);
t27 = t12 * mrSges(5,3);
t31 = cos(qJ(5));
t33 = sin(qJ(5));
t35 = -t12 * t33 - t22 * t31;
t36 = t35 * mrSges(6,1);
t39 = -t12 * t31 + t22 * t33;
t40 = t39 * mrSges(6,2);
t41 = t25 * mrSges(6,3);
t44 = t39 * pkin(6);
t48 = pkin(7) * qJ(5) - qJ(6);
t49 = cos(t48);
t51 = sin(t48);
t54 = (t25 * t51 - t35 * t49) * mrSges(7,1);
t58 = (-t25 * t49 - t35 * t51) * mrSges(7,2);
t59 = t39 * mrSges(7,3);
t65 = -t16 - t17;
t77 = t12 * t24;
t79 = t12 * t21;
t111 = -t10 * t31 - t79 * t33;
t115 = t10 * t33 - t79 * t31;
t122 = -t115 * pkin(7);
unknown(1) = 0;
unknown(2) = (-g(3) * (-t3 * mrSges(3,1) - t1 * mrSges(3,2)) - g(3) * (-m(4) * t3 * pkin(3) - t11 - t13) - g(3) * (m(5) * t19 - t23 + t26 - t27) - g(3) * (m(6) * t19 + t36 + t40 - t41) - g(3) * (m(7) * (t44 - t16 - t17 - t18) + t54 + t58 + t59));
unknown(3) = (-g(3) * (-t11 - t13) - g(3) * (m(5) * t65 - t23 + t26 - t27) - g(3) * (m(6) * t65 + t36 + t40 - t41) - g(3) * (m(7) * (t44 - t16 - t17) + t54 + t58 + t59));
unknown(4) = (-g(3) * (-t77 * mrSges(5,1) - t79 * mrSges(5,2)) - g(3) * (-t77 * t31 * mrSges(6,1) + t77 * t33 * mrSges(6,2) + t79 * mrSges(6,3)) - g(3) * (m(7) * t12 * t24 * t33 * pkin(6) + (t77 * t31 * t49 - t79 * t51) * mrSges(7,1) + (t77 * t31 * t51 + t79 * t49) * mrSges(7,2) + t77 * t33 * mrSges(7,3)));
unknown(5) = (-g(3) * (t111 * mrSges(6,1) + t115 * mrSges(6,2)) - g(3) * (m(7) * t115 * pkin(6) + (-t77 * pkin(7) * t49 - t111 * t49 + t122 * t51) * mrSges(7,1) + (-t77 * pkin(7) * t51 - t111 * t51 - t122 * t49) * mrSges(7,2) + t115 * mrSges(7,3)));
unknown(6) = -(g(3) * ((t115 * t51 + t77 * t49) * mrSges(7,1) + (-t115 * t49 + t77 * t51) * mrSges(7,2)));
taug = unknown(:);
