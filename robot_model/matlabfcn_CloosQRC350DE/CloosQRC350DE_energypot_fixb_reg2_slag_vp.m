% Calculate inertial parameters regressor of potential energy for
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
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = CloosQRC350DE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:33
% EndTime: 2020-06-19 21:39:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (65->30), mult. (61->29), div. (0->0), fcn. (59->10), ass. (0->85)
unknown=NaN(1,60);
t1 = (g(3) * pkin(1));
t2 = cos(qJ(2));
t4 = sin(qJ(2));
t6 = qJ(2) + qJ(3);
t7 = cos(t6);
t8 = g(3) * t7;
t9 = sin(t6);
t10 = g(3) * t9;
t11 = t2 * pkin(3);
t14 = cos(qJ(4));
t16 = sin(qJ(4));
t17 = t8 * t16;
t18 = t7 * pkin(4);
t19 = t9 * pkin(5);
t21 = g(3) * (t18 - t19 + t11 + pkin(1));
t22 = t7 * t14;
t23 = cos(qJ(5));
t25 = sin(qJ(5));
t27 = t22 * t23 - t9 * t25;
t31 = -t22 * t25 - t9 * t23;
t32 = g(3) * t31;
t34 = pkin(7) * qJ(5) - qJ(6);
t35 = cos(t34);
t37 = t7 * t16;
t38 = sin(t34);
unknown(1,1) = 0;
unknown(1,2) = 0;
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = 0;
unknown(1,6) = 0;
unknown(1,7) = 0;
unknown(1,8) = 0;
unknown(1,9) = -g(3);
unknown(1,10) = -t1;
unknown(1,11) = 0;
unknown(1,12) = 0;
unknown(1,13) = 0;
unknown(1,14) = 0;
unknown(1,15) = 0;
unknown(1,16) = 0;
unknown(1,17) = -(g(3) * t2);
unknown(1,18) = (g(3) * t4);
unknown(1,19) = 0;
unknown(1,20) = -t1;
unknown(1,21) = 0;
unknown(1,22) = 0;
unknown(1,23) = 0;
unknown(1,24) = 0;
unknown(1,25) = 0;
unknown(1,26) = 0;
unknown(1,27) = -t8;
unknown(1,28) = t10;
unknown(1,29) = 0;
unknown(1,30) = -(g(3) * (t11 + pkin(1)));
unknown(1,31) = 0;
unknown(1,32) = 0;
unknown(1,33) = 0;
unknown(1,34) = 0;
unknown(1,35) = 0;
unknown(1,36) = 0;
unknown(1,37) = -(t8 * t14);
unknown(1,38) = t17;
unknown(1,39) = t10;
unknown(1,40) = -t21;
unknown(1,41) = 0;
unknown(1,42) = 0;
unknown(1,43) = 0;
unknown(1,44) = 0;
unknown(1,45) = 0;
unknown(1,46) = 0;
unknown(1,47) = -(g(3) * t27);
unknown(1,48) = -t32;
unknown(1,49) = -t17;
unknown(1,50) = -t21;
unknown(1,51) = 0;
unknown(1,52) = 0;
unknown(1,53) = 0;
unknown(1,54) = 0;
unknown(1,55) = 0;
unknown(1,56) = 0;
unknown(1,57) = -(g(3) * (-t27 * t35 - t37 * t38));
unknown(1,58) = -(g(3) * (-t27 * t38 + t37 * t35));
unknown(1,59) = -t32;
unknown(1,60) = -(g(3) * (t31 * pkin(6) + pkin(1) + t11 + t18 - t19));
U_reg = unknown;
