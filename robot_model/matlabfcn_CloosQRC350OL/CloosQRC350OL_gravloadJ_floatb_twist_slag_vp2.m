% Calculate Gravitation load on the joints for
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
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350OL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:00:28
% EndTime: 2020-06-20 08:00:29
% DurationCPUTime: 0.34s
% Computational Cost: add. (164->39), mult. (190->53), div. (0->0), fcn. (194->10), ass. (0->22)
t36 = -mrSges(5,2) + mrSges(6,3);
t35 = -pkin(6) * m(7) - mrSges(6,2) - mrSges(7,3);
t18 = sin(qJ(6));
t22 = cos(qJ(6));
t26 = mrSges(7,1) * t22 - mrSges(7,2) * t18 - mrSges(6,1);
t33 = -m(5) - m(6);
t20 = sin(qJ(4));
t32 = t18 * t20;
t19 = sin(qJ(5));
t24 = cos(qJ(4));
t31 = t19 * t24;
t30 = t20 * t22;
t23 = cos(qJ(5));
t29 = t23 * t24;
t17 = qJ(2) + qJ(3);
t15 = sin(t17);
t16 = cos(t17);
t27 = -t15 * pkin(4) - t16 * pkin(5);
t25 = (m(7) * pkin(5) + mrSges(4,2) + mrSges(5,3)) * t16 + t26 * (-t15 * t29 - t16 * t19) + (m(7) * pkin(4) + mrSges(5,1) * t24 + t32 * mrSges(7,1) + t30 * mrSges(7,2) + t36 * t20 + mrSges(4,1)) * t15 + t35 * (t15 * t31 - t16 * t23);
t21 = sin(qJ(2));
t10 = -t15 * t19 + t16 * t29;
t1 = [0, (cos(qJ(2)) * mrSges(3,2) + (mrSges(3,1) + (m(4) + m(7)) * pkin(3)) * t21 + t33 * (-t21 * pkin(3) + t27) + t25) * g(3), (t33 * t27 + t25) * g(3), ((-t18 * mrSges(7,1) - t22 * mrSges(7,2) - t36) * t24 + (t19 * t35 - t26 * t23 + mrSges(5,1)) * t20) * g(3) * t16, (t26 * (-t15 * t23 - t16 * t31) - t35 * t10) * g(3), -g(3) * ((t10 * t18 + t16 * t30) * mrSges(7,1) + (t10 * t22 - t16 * t32) * mrSges(7,2))];
taug = t1(:);
