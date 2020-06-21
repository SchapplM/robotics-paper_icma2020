% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x36]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = CloosQRC350DE_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:33
% EndTime: 2020-06-19 21:39:33
% DurationCPUTime: 0.07s
% Computational Cost: add. (31->13), mult. (34->21), div. (0->0), fcn. (36->10), ass. (0->52)
unknown=NaN(1,36);
t1 = cos(qJ(2));
t3 = sin(qJ(2));
t5 = qJ(2) + qJ(3);
t6 = cos(t5);
t7 = g(3) * t6;
t8 = sin(t5);
t10 = cos(qJ(4));
t12 = sin(qJ(4));
t14 = t6 * t10;
t15 = cos(qJ(5));
t17 = sin(qJ(5));
t19 = t14 * t15 - t8 * t17;
t26 = pkin(7) * qJ(5) - qJ(6);
t27 = cos(t26);
t29 = t6 * t12;
t30 = sin(t26);
unknown(1,1) = 0;
unknown(1,2) = 0;
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = 0;
unknown(1,6) = 0;
unknown(1,7) = -(g(3) * t1);
unknown(1,8) = (g(3) * t3);
unknown(1,9) = 0;
unknown(1,10) = 0;
unknown(1,11) = 0;
unknown(1,12) = 0;
unknown(1,13) = 0;
unknown(1,14) = -t7;
unknown(1,15) = (g(3) * t8);
unknown(1,16) = 0;
unknown(1,17) = 0;
unknown(1,18) = 0;
unknown(1,19) = 0;
unknown(1,20) = 0;
unknown(1,21) = -(t7 * t10);
unknown(1,22) = (t7 * t12);
unknown(1,23) = 0;
unknown(1,24) = 0;
unknown(1,25) = 0;
unknown(1,26) = 0;
unknown(1,27) = 0;
unknown(1,28) = -(g(3) * t19);
unknown(1,29) = -(g(3) * (-t14 * t17 - t8 * t15));
unknown(1,30) = 0;
unknown(1,31) = 0;
unknown(1,32) = 0;
unknown(1,33) = 0;
unknown(1,34) = 0;
unknown(1,35) = -(g(3) * (-t19 * t27 - t29 * t30));
unknown(1,36) = -(g(3) * (-t19 * t30 + t29 * t27));
U_reg = unknown;
