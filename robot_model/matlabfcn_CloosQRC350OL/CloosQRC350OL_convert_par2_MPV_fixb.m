% Return the minimum parameter vector for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
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
% MPV [19x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = CloosQRC350OL_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'CloosQRC350OL_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'CloosQRC350OL_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'CloosQRC350OL_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t17 = (pkin(6) * m(7));
t5 = (m(5) + m(6) + m(7));
t16 = (t5 * pkin(5));
t4 = m(4) + t5;
t1 = t4 * pkin(3) ^ 2;
t15 = -Ifges(3,2) - t1;
t3 = (t5 * pkin(4) ^ 2);
t14 = (-Ifges(4,2) - t3);
t13 = Ifges(5,2) + (2 * mrSges(5,3) + t16) * pkin(5);
t12 = Ifges(7,2) + (2 * mrSges(7,3) + t17) * pkin(6);
t11 = -mrSges(5,3) - t16;
t2 = [(m(3) + t4) * pkin(2) ^ 2 + Ifges(2,3) - t14 - t15; Ifges(3,1) + t15; -pkin(3) * mrSges(4,3) + Ifges(3,5); t1 + Ifges(3,3); t4 * pkin(3) + mrSges(3,1); Ifges(4,1) + t13 + t14; t11 * pkin(4) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + t3 + t13; t5 * pkin(4) + mrSges(4,1); mrSges(4,2) - t11; Ifges(5,1) - Ifges(5,2) + Ifges(6,2); Ifges(5,3) + Ifges(6,2); Ifges(6,1) - Ifges(6,2) + t12; Ifges(6,3) + t12; mrSges(6,2) + mrSges(7,3) + t17; Ifges(7,1) - Ifges(7,2); Ifges(7,3);];
MPV = t2;
