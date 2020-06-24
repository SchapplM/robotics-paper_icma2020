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
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
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
% StartTime: 2020-06-23 21:56:05
% EndTime: 2020-06-23 21:56:05
% DurationCPUTime: 0.19s
% Computational Cost: add. (79->18), mult. (90->24), div. (0->0), fcn. (77->7), ass. (0->15)
t23 = -pkin(6) * m(7) - mrSges(6,2) - mrSges(7,3);
t10 = qJ(2) + qJ(3);
t8 = sin(t10);
t22 = pkin(4) * t8;
t21 = -m(5) - m(6);
t11 = sin(qJ(5));
t20 = t8 * t11;
t9 = cos(t10);
t19 = t9 * cos(qJ(5));
t18 = -pkin(5) * t9 - t22;
t17 = t23 * g(3);
t15 = cos(qJ(4));
t16 = t8 * mrSges(4,1) + m(7) * t22 + (pkin(5) * m(7) + mrSges(4,2) + mrSges(5,3)) * t9 + t23 * (t15 * t20 - t19);
t13 = sin(qJ(2));
t1 = [0, (t21 * (-pkin(3) * t13 + t18) + (mrSges(3,1) + (m(4) + m(7)) * pkin(3)) * t13 + t16) * g(3), (t21 * t18 + t16) * g(3), t9 * sin(qJ(4)) * t11 * t17, (-t15 * t19 + t20) * t17, 0];
taug = t1(:);
