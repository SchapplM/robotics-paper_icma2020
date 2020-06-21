% Return the minimum parameter vector for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
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
% MPV [36x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = CloosQRC350DE_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_convert_par2_MPV_fixb: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350DE_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350DE_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350DE_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
unknown=NaN(36,1);
t1 = (pkin(4) ^ 2);
t2 = (m(5) + m(6) + m(7));
t3 = (t1 * t2);
t4 = pkin(3) ^ 2;
t5 = m(4) + m(5) + m(6) + m(7);
t6 = t4 * t5;
t7 = pkin(2) ^ 2;
t17 = (pkin(5) ^ 2);
t20 = 2 * pkin(5) * mrSges(5,3);
t36 = pkin(6) ^ 2;
t37 = (t36 * m(7));
t39 = 2 * pkin(6) * mrSges(7,3);
unknown(1,1) = Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t3 + t6 + t7 * (m(3) + m(4) + m(5) + m(6) + m(7));
unknown(2,1) = Ifges(3,1) - Ifges(3,2) - t6;
unknown(3,1) = Ifges(3,4);
unknown(4,1) = -pkin(3) * mrSges(4,3) + Ifges(3,5);
unknown(5,1) = Ifges(3,6);
unknown(6,1) = Ifges(3,3) + t6;
unknown(7,1) = pkin(3) * t5 + mrSges(3,1);
unknown(8,1) = mrSges(3,2);
unknown(9,1) = t17 * t2 + Ifges(4,1) - Ifges(4,2) + Ifges(5,2) + t20 - t3;
unknown(10,1) = -pkin(4) * pkin(5) * t2 - pkin(4) * mrSges(5,3) + Ifges(4,4);
unknown(11,1) = Ifges(4,5);
unknown(12,1) = Ifges(4,6);
unknown(13,1) = Ifges(4,3) + Ifges(5,2) + t20 + (t1 + t17) * t2;
unknown(14,1) = pkin(4) * t2 + mrSges(4,1);
unknown(15,1) = pkin(5) * t2 + mrSges(4,2) + mrSges(5,3);
unknown(16,1) = Ifges(5,1) + Ifges(6,2) - Ifges(5,2);
unknown(17,1) = Ifges(5,4);
unknown(18,1) = Ifges(5,5);
unknown(19,1) = Ifges(5,6);
unknown(20,1) = Ifges(5,3) + Ifges(6,2);
unknown(21,1) = mrSges(5,1);
unknown(22,1) = mrSges(5,2) - mrSges(6,3);
unknown(23,1) = t37 + t39 + Ifges(6,1) - Ifges(6,2) + Ifges(7,2);
unknown(24,1) = Ifges(6,4);
unknown(25,1) = Ifges(6,5);
unknown(26,1) = Ifges(6,6);
unknown(27,1) = t37 + t39 + Ifges(7,2) + Ifges(6,3);
unknown(28,1) = mrSges(6,1);
unknown(29,1) = pkin(6) * m(7) + mrSges(6,2) + mrSges(7,3);
unknown(30,1) = Ifges(7,1) - Ifges(7,2);
unknown(31,1) = Ifges(7,4);
unknown(32,1) = Ifges(7,5);
unknown(33,1) = Ifges(7,6);
unknown(34,1) = Ifges(7,3);
unknown(35,1) = mrSges(7,1);
unknown(36,1) = mrSges(7,2);
MPV = unknown;
