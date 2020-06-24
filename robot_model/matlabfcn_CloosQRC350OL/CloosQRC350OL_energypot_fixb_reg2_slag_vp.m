% Calculate inertial parameters regressor of potential energy for
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
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = CloosQRC350OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:02:20
% EndTime: 2020-06-23 22:02:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->17), mult. (32->15), div. (0->0), fcn. (29->6), ass. (0->12)
t63 = g(3) * pkin(1);
t58 = qJ(2) + qJ(3);
t55 = sin(t58);
t56 = cos(t58);
t51 = -t56 * cos(qJ(4)) * sin(qJ(5)) - t55 * cos(qJ(5));
t62 = g(3) * t51;
t59 = cos(qJ(2));
t60 = t59 * pkin(3) + pkin(1);
t52 = t56 * pkin(4) - t55 * pkin(5) + t60;
t61 = g(3) * t52;
t54 = g(3) * t55;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, -g(3), -t63, 0, 0, 0, 0, 0, 0, -g(3) * t59, 0, 0, -t63, 0, 0, 0, 0, 0, 0, -g(3) * t56, t54, 0, -g(3) * t60, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t61, 0, 0, 0, 0, 0, 0, 0, -t62, 0, -t61, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -g(3) * (t51 * pkin(6) + t52);];
U_reg = t1;
