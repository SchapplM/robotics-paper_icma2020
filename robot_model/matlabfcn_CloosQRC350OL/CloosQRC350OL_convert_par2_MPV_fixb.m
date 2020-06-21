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
% MPV [36x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
t94 = (pkin(6) * m(7));
t83 = (m(5) + m(6) + m(7));
t82 = (m(4) + t83);
t81 = pkin(3) ^ 2 * t82;
t93 = (-Ifges(3,2) - t81);
t92 = 2 * pkin(5) * mrSges(5,3) + Ifges(5,2);
t91 = Ifges(7,2) + (2 * mrSges(7,3) + t94) * pkin(6);
t90 = -pkin(5) * t83 - mrSges(5,3);
t88 = (pkin(4) ^ 2);
t87 = pkin(5) ^ 2;
t1 = [Ifges(2,3) + Ifges(4,2) + t88 * t83 + pkin(2) ^ 2 * (m(3) + t82) - t93; Ifges(3,1) + t93; Ifges(3,4); -pkin(3) * mrSges(4,3) + Ifges(3,5); Ifges(3,6); Ifges(3,3) + t81; pkin(3) * t82 + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) + (t87 - t88) * t83 + t92; t90 * pkin(4) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t87 + t88) * t83 + t92; pkin(4) * t83 + mrSges(4,1); mrSges(4,2) - t90; Ifges(5,1) + Ifges(6,2) - Ifges(5,2); Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + Ifges(6,2); mrSges(5,1); mrSges(5,2) - mrSges(6,3); Ifges(6,1) - Ifges(6,2) + t91; Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + t91; mrSges(6,1); mrSges(6,2) + mrSges(7,3) + t94; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV = t1;
