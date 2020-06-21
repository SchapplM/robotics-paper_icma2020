% Calculate minimal parameter regressor of potential energy for
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
% 
% Output:
% U_reg [1x36]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = CloosQRC350OL_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:18:55
% EndTime: 2020-06-20 08:18:56
% DurationCPUTime: 0.06s
% Computational Cost: add. (27->12), mult. (30->20), div. (0->0), fcn. (36->10), ass. (0->11)
t145 = qJ(2) + qJ(3);
t144 = cos(t145);
t153 = t144 * sin(qJ(4));
t152 = t144 * cos(qJ(4));
t150 = cos(qJ(5));
t149 = cos(qJ(6));
t147 = sin(qJ(5));
t146 = sin(qJ(6));
t143 = sin(t145);
t142 = -t143 * t147 + t150 * t152;
t1 = [0, 0, 0, 0, 0, 0, -g(3) * cos(qJ(2)), g(3) * sin(qJ(2)), 0, 0, 0, 0, 0, -g(3) * t144, g(3) * t143, 0, 0, 0, 0, 0, -g(3) * t152, g(3) * t153, 0, 0, 0, 0, 0, -g(3) * t142, -g(3) * (-t143 * t150 - t147 * t152), 0, 0, 0, 0, 0, -g(3) * (-t142 * t149 + t146 * t153), -g(3) * (t142 * t146 + t149 * t153);];
U_reg = t1;
