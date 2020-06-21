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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 08:00:35
% EndTime: 2020-06-20 08:00:53
% DurationCPUTime: 9.74s
% Computational Cost: add. (57110->174), mult. (115827->236), div. (0->0), fcn. (88050->12), ass. (0->97)
t100 = cos(qJ(5));
t101 = cos(qJ(4));
t103 = cos(qJ(2));
t118 = qJD(1) * qJD(2);
t114 = t103 * t118;
t97 = sin(qJ(2));
t82 = -t97 * qJDD(1) - t114;
t92 = qJDD(1) * pkin(2);
t113 = t92 + (t114 - t82) * pkin(3);
t102 = cos(qJ(3));
t96 = sin(qJ(3));
t76 = (t102 * t103 - t96 * t97) * qJD(1);
t83 = t103 * qJDD(1) - t97 * t118;
t55 = -t76 * qJD(3) + t102 * t82 - t96 * t83;
t75 = (-t102 * t97 - t103 * t96) * qJD(1);
t56 = t75 * qJD(3) + t102 * t83 + t96 * t82;
t91 = qJD(2) + qJD(3);
t30 = (t75 * t91 + t56) * pkin(5) + (t76 * t91 - t55) * pkin(4) + t113;
t105 = qJD(1) ^ 2;
t120 = t103 * t105;
t112 = -pkin(2) * t120 + t97 * g(3);
t71 = (-t97 * t120 + qJDD(2)) * pkin(3) + t112;
t111 = -t97 * t105 * pkin(2) - t103 * g(3);
t73 = (-t105 * t97 ^ 2 - qJD(2) ^ 2) * pkin(3) + t111;
t122 = t102 * t73 + t96 * t71;
t61 = -t75 * pkin(4) + t76 * pkin(5);
t88 = t91 ^ 2;
t90 = qJDD(2) + qJDD(3);
t37 = -t88 * pkin(4) - t90 * pkin(5) + t75 * t61 + t122;
t95 = sin(qJ(4));
t24 = t101 * t37 - t95 * t30;
t115 = t102 * t71 - t96 * t73;
t36 = t90 * pkin(4) - t88 * pkin(5) - t76 * t61 + t115;
t94 = sin(qJ(5));
t123 = t100 * t24 + t94 * t36;
t121 = qJD(1) * t97;
t119 = qJD(1) * t103;
t60 = -t75 * mrSges(4,1) + t76 * mrSges(4,2);
t69 = t91 * mrSges(4,1) - t76 * mrSges(4,3);
t64 = t101 * t76 - t95 * t91;
t41 = -t64 * qJD(4) - t101 * t90 - t95 * t56;
t40 = qJDD(5) - t41;
t74 = qJD(4) + t75;
t51 = t100 * t74 - t94 * t64;
t52 = t100 * t64 + t94 * t74;
t16 = (t51 * t52 - t40) * pkin(6) + t123;
t23 = t101 * t30 + t95 * t37;
t63 = -t101 * t91 - t95 * t76;
t42 = t63 * qJD(4) + t101 * t56 - t95 * t90;
t54 = qJDD(4) + t55;
t28 = t51 * qJD(5) + t100 * t42 + t94 * t54;
t62 = qJD(5) - t63;
t17 = (t51 * t62 + t28) * pkin(6) + t23;
t93 = sin(qJ(6));
t99 = cos(qJ(6));
t43 = t93 * t52 + t99 * t62;
t21 = t43 * qJD(6) - t99 * t28 + t93 * t40;
t27 = -t52 * qJD(5) + t100 * t54 - t94 * t42;
t26 = qJDD(6) + t27;
t44 = -t99 * t52 + t93 * t62;
t29 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t50 = qJD(6) + t51;
t31 = -t50 * mrSges(7,2) + t43 * mrSges(7,3);
t14 = m(7) * (t93 * t16 + t99 * t17) - t21 * mrSges(7,3) + t26 * mrSges(7,1) - t44 * t29 + t50 * t31;
t20 = -t44 * qJD(6) + t93 * t28 + t99 * t40;
t32 = t50 * mrSges(7,1) - t44 * mrSges(7,3);
t15 = m(7) * (-t99 * t16 + t93 * t17) + t20 * mrSges(7,3) - t26 * mrSges(7,2) + t43 * t29 - t50 * t32;
t38 = -t51 * mrSges(6,1) + t52 * mrSges(6,2);
t46 = t62 * mrSges(6,1) - t52 * mrSges(6,3);
t11 = m(6) * t123 - t40 * mrSges(6,2) + t27 * mrSges(6,3) + t93 * t14 - t99 * t15 + t51 * t38 - t62 * t46;
t116 = t100 * t36 - t94 * t24;
t110 = -t20 * mrSges(7,1) - t43 * t31 + m(7) * ((-t52 ^ 2 - t62 ^ 2) * pkin(6) + t116) + t21 * mrSges(7,2) + t44 * t32;
t45 = -t62 * mrSges(6,2) + t51 * mrSges(6,3);
t13 = m(6) * t116 + t40 * mrSges(6,1) - t28 * mrSges(6,3) - t52 * t38 + t62 * t45 + t110;
t48 = -t63 * mrSges(5,1) + t64 * mrSges(5,2);
t58 = t74 * mrSges(5,1) - t64 * mrSges(5,3);
t8 = m(5) * t24 - t54 * mrSges(5,2) + t41 * mrSges(5,3) + t100 * t11 - t94 * t13 + t63 * t48 - t74 * t58;
t107 = t27 * mrSges(6,1) - t28 * mrSges(6,2) - t99 * t14 - t93 * t15 + t51 * t45 - t52 * t46;
t57 = -t74 * mrSges(5,2) + t63 * mrSges(5,3);
t9 = t54 * mrSges(5,1) - t42 * mrSges(5,3) - t64 * t48 + t74 * t57 + (-m(5) - m(6)) * t23 + t107;
t6 = m(4) * t122 - t90 * mrSges(4,2) + t55 * mrSges(4,3) + t101 * t8 + t75 * t60 - t91 * t69 - t95 * t9;
t109 = m(5) * t36 - t41 * mrSges(5,1) + t42 * mrSges(5,2) + t100 * t13 + t94 * t11 - t63 * t57 + t64 * t58;
t68 = -t91 * mrSges(4,2) + t75 * mrSges(4,3);
t7 = m(4) * t115 + t90 * mrSges(4,1) - t56 * mrSges(4,3) - t76 * t60 + t91 * t68 + t109;
t81 = (mrSges(3,1) * t97 + mrSges(3,2) * t103) * qJD(1);
t85 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t119;
t3 = m(3) * t111 - qJDD(2) * mrSges(3,2) + t82 * mrSges(3,3) - qJD(2) * t85 + t102 * t6 - t81 * t121 - t96 * t7;
t84 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t121;
t4 = m(3) * t112 + qJDD(2) * mrSges(3,1) - t83 * mrSges(3,3) + qJD(2) * t84 + t102 * t7 - t81 * t119 + t96 * t6;
t117 = t103 * t3 - t97 * t4;
t108 = m(4) * t113 - t55 * mrSges(4,1) + t56 * mrSges(4,2) - t101 * t9 - t75 * t68 + t76 * t69 - t95 * t8;
t106 = m(3) * t92 - t82 * mrSges(3,1) + t83 * mrSges(3,2) + t85 * t119 + t84 * t121 + t108;
t104 = cos(qJ(1));
t98 = sin(qJ(1));
t5 = qJDD(1) * mrSges(2,1) - t105 * mrSges(2,2) + t106;
t1 = -t105 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t103 * t4 + t97 * t3;
t2 = [t104 * t1 - t98 * t5, t1, t3, t6, t8, t11, t15; t98 * t1 + t104 * t5, t5, t4, t7, t9, t13, t14; (-m(1) - m(2)) * g(3) + t117, -m(2) * g(3) + t117, t106, t108, t109, m(6) * t23 - t107, t110;];
f_new = t2;
