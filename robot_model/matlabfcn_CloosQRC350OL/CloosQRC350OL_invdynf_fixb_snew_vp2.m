% Calculate vector of cutting forces with Newton-Euler
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
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
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = CloosQRC350OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'CloosQRC350OL_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:56:10
% EndTime: 2020-06-23 21:56:19
% DurationCPUTime: 6.56s
% Computational Cost: add. (38344->135), mult. (81253->185), div. (0->0), fcn. (57700->12), ass. (0->86)
t79 = sin(qJ(3));
t80 = sin(qJ(2));
t85 = cos(qJ(3));
t86 = cos(qJ(2));
t61 = (-t79 * t86 - t80 * t85) * qJD(1);
t106 = qJD(1) * qJD(2);
t104 = t86 * t106;
t65 = -qJDD(1) * t80 - t104;
t103 = t80 * t106;
t66 = qJDD(1) * t86 - t103;
t41 = qJD(3) * t61 + t65 * t79 + t66 * t85;
t62 = (-t79 * t80 + t85 * t86) * qJD(1);
t73 = qJD(2) + qJD(3);
t78 = sin(qJ(4));
t84 = cos(qJ(4));
t49 = -t62 * t78 - t73 * t84;
t72 = qJDD(2) + qJDD(3);
t30 = qJD(4) * t49 + t41 * t84 - t72 * t78;
t50 = t62 * t84 - t73 * t78;
t60 = qJD(4) + t61;
t77 = sin(qJ(5));
t83 = cos(qJ(5));
t37 = -t50 * t77 + t83 * t60;
t40 = -qJD(3) * t62 + t85 * t65 - t66 * t79;
t22 = t83 * t30 + t77 * (qJDD(4) + t40) + t37 * qJD(5);
t46 = qJD(5) - t49;
t101 = t37 * t46 + t22;
t93 = -t78 * t41 - t84 * t72;
t29 = qJD(4) * t50 + qJDD(5) - t93;
t38 = t50 * t83 + t60 * t77;
t100 = t37 * t38 - t29;
t75 = qJDD(1) * pkin(2);
t97 = t75 + (t104 - t65) * pkin(3);
t23 = (t61 * t73 + t41) * pkin(5) + (t62 * t73 - t40) * pkin(4) + t97;
t88 = qJD(1) ^ 2;
t112 = t86 * t88;
t96 = -t80 * t112 + qJDD(2);
t98 = -pkin(2) * t112 + t80 * g(3);
t57 = t96 * pkin(3) + t98;
t74 = t80 ^ 2;
t113 = t74 * t88;
t92 = -pkin(2) * t80 * t88 - g(3) * t86;
t59 = (-qJD(2) ^ 2 - t113) * pkin(3) + t92;
t109 = t79 * t57 + t85 * t59;
t44 = -pkin(4) * t61 + pkin(5) * t62;
t70 = t73 ^ 2;
t28 = -pkin(4) * t70 - pkin(5) * t72 + t44 * t61 + t109;
t21 = -t23 * t78 + t28 * t84;
t99 = t85 * t57 - t59 * t79;
t27 = pkin(4) * t72 - pkin(5) * t70 - t44 * t62 + t99;
t111 = t83 * t21 + t77 * t27;
t16 = t100 * pkin(6) + t111;
t20 = t23 * t84 + t28 * t78;
t17 = t101 * pkin(6) + t20;
t76 = sin(qJ(6));
t82 = cos(qJ(6));
t33 = t38 * t76 + t46 * t82;
t12 = m(7) * (t16 * t76 + t17 * t82) + (t82 * t22 - t76 * t29 + t33 * t37) * mrSges(7,3);
t34 = -t38 * t82 + t46 * t76;
t13 = m(7) * (-t16 * t82 + t17 * t76) + (t76 * t22 + t82 * t29 + t34 * t37) * mrSges(7,3);
t114 = -t101 * mrSges(6,2) - t12 * t82 - t13 * t76;
t110 = -t38 ^ 2 - t46 ^ 2;
t108 = qJD(1) * t86;
t43 = -mrSges(4,1) * t61 + mrSges(4,2) * t62;
t55 = mrSges(4,1) * t73 - mrSges(4,3) * t62;
t11 = m(6) * t111 + t100 * mrSges(6,2) + t76 * t12 - t82 * t13;
t102 = -t21 * t77 + t83 * t27;
t95 = m(7) * (t110 * pkin(6) + t102) + (-t33 ^ 2 - t34 ^ 2) * mrSges(7,3);
t15 = m(6) * t102 + t110 * mrSges(6,2) + t95;
t7 = m(5) * t21 + t83 * t11 - t77 * t15 + ((-qJD(4) + t60) * t50 + t93) * mrSges(5,3);
t9 = (-m(5) - m(6)) * t20 + (t49 * t60 - t30) * mrSges(5,3) + t114;
t6 = m(4) * t109 - t72 * mrSges(4,2) + t40 * mrSges(4,3) + t61 * t43 - t73 * t55 + t84 * t7 - t78 * t9;
t54 = -mrSges(4,2) * t73 + mrSges(4,3) * t61;
t91 = t77 * t11 + t83 * t15 + m(5) * t27 + (-t49 ^ 2 - t50 ^ 2) * mrSges(5,3);
t8 = m(4) * t99 + t72 * mrSges(4,1) - t41 * mrSges(4,3) - t62 * t43 + t73 * t54 + t91;
t3 = m(3) * t98 + t79 * t6 + t85 * t8 + (-t66 - t103) * mrSges(3,3) + t96 * mrSges(3,1);
t67 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t108;
t4 = m(3) * t92 - mrSges(3,1) * t113 + t65 * mrSges(3,3) - qJD(2) * t67 + t85 * t6 - t79 * t8;
t105 = -t3 * t80 + t86 * t4;
t90 = m(4) * t97 - mrSges(4,1) * t40 + t41 * mrSges(4,2) - t54 * t61 + t62 * t55 - t7 * t78 - t84 * t9;
t89 = m(3) * t75 - mrSges(3,1) * t65 + t67 * t108 + t90;
t87 = cos(qJ(1));
t81 = sin(qJ(1));
t5 = qJDD(1) * mrSges(2,1) + (-mrSges(3,3) * t74 - mrSges(2,2)) * t88 + t89;
t1 = -mrSges(2,1) * t88 - qJDD(1) * mrSges(2,2) + t3 * t86 + t4 * t80;
t2 = [t1 * t87 - t5 * t81, t1, t4, t6, t7, t11, t13; t1 * t81 + t5 * t87, t5, t3, t8, t9, t15, t12; (-m(1) - m(2)) * g(3) + t105, -m(2) * g(3) + t105, -mrSges(3,3) * t113 + t89, t90, t91, m(6) * t20 - t114, t95;];
f_new = t2;
