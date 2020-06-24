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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
% OptimizationMode: 2
% StartTime: 2020-06-23 21:02:09
% EndTime: 2020-06-23 21:02:09
% DurationCPUTime: 0.18s
% Computational Cost: add. (72->20), mult. (72->30), div. (0->0), fcn. (42->8), ass. (0->16)
t15 = m(7) + m(5) + m(6);
t11 = cos(qJ(3));
t10 = cos(qJ(4));
t5 = pkin(6) * m(7) + mrSges(6,2) + mrSges(7,3);
t6 = sin(qJ(5));
t17 = t5 * t6;
t3 = t15 * pkin(4) - t10 * t17 + mrSges(4,1);
t9 = cos(qJ(5));
t4 = t15 * pkin(5) + t5 * t9 + mrSges(4,2) + mrSges(5,3);
t7 = sin(qJ(3));
t18 = t3 * t11 - t4 * t7;
t16 = t10 * t9;
t12 = cos(qJ(2));
t8 = sin(qJ(2));
t1 = (t4 * t11 + t3 * t7) * t12;
t2 = [0, g(3) * (t1 + (mrSges(3,1) + (m(4) + t15) * pkin(3) + t18) * t8), g(3) * (t8 * t18 + t1), -g(3) * sin(qJ(4)) * (t12 * t11 - t8 * t7) * t17, ((t11 * t16 - t6 * t7) * t12 - t8 * (t6 * t11 + t7 * t16)) * g(3) * t5, 0];
taug = t2(:);
