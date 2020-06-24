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
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
% OptimizationMode: 2
% StartTime: 2020-06-23 21:06:21
% EndTime: 2020-06-23 21:06:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (34->21), mult. (60->30), div. (0->0), fcn. (63->7), ass. (0->15)
t68 = pkin(1) * g(3);
t67 = sin(qJ(5)) * cos(qJ(4));
t66 = cos(qJ(2));
t65 = cos(qJ(3));
t63 = cos(qJ(5));
t62 = sin(qJ(2));
t61 = sin(qJ(3));
t59 = pkin(6) * t63 + pkin(5);
t58 = pkin(6) * t67 - pkin(4);
t57 = t66 * t61 + t65 * t62;
t56 = t62 * t61 - t66 * t65;
t55 = t57 * g(3);
t54 = ((-pkin(4) * t65 + pkin(5) * t61 - pkin(3)) * t66 - pkin(1) + (pkin(4) * t61 + pkin(5) * t65) * t62) * g(3);
t53 = (-t56 * t67 + t63 * t57) * g(3);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, -g(3), -t68, 0, 0, 0, 0, 0, 0, -t66 * g(3), 0, 0, -t68, 0, 0, 0, 0, 0, 0, t56 * g(3), t55, 0, -(t66 * pkin(3) + pkin(1)) * g(3), 0, 0, 0, 0, 0, 0, 0, 0, t55, t54, 0, 0, 0, 0, 0, 0, 0, t53, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(3) * ((t58 * t65 + t59 * t61 - pkin(3)) * t66 - pkin(1) + (-t61 * t58 + t59 * t65) * t62);];
U_reg = t1;
