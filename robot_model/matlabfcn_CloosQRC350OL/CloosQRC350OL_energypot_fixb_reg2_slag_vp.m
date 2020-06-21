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
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
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
% StartTime: 2020-06-20 08:18:55
% EndTime: 2020-06-20 08:18:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (61->25), mult. (57->28), div. (0->0), fcn. (59->10), ass. (0->20)
t80 = g(3) * pkin(1);
t66 = qJ(2) + qJ(3);
t63 = sin(t66);
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t64 = cos(t66);
t75 = t64 * cos(qJ(4));
t58 = -t63 * t71 - t68 * t75;
t79 = g(3) * t58;
t73 = cos(qJ(2));
t77 = t73 * pkin(3) + pkin(1);
t60 = t64 * pkin(4) - t63 * pkin(5) + t77;
t78 = g(3) * t60;
t76 = t64 * sin(qJ(4));
t74 = g(3) * t76;
t70 = cos(qJ(6));
t67 = sin(qJ(6));
t62 = g(3) * t63;
t59 = -t63 * t68 + t71 * t75;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, -g(3), -t80, 0, 0, 0, 0, 0, 0, -g(3) * t73, g(3) * sin(qJ(2)), 0, -t80, 0, 0, 0, 0, 0, 0, -g(3) * t64, t62, 0, -g(3) * t77, 0, 0, 0, 0, 0, 0, -g(3) * t75, t74, t62, -t78, 0, 0, 0, 0, 0, 0, -g(3) * t59, -t79, -t74, -t78, 0, 0, 0, 0, 0, 0, -g(3) * (-t59 * t70 + t67 * t76), -g(3) * (t59 * t67 + t70 * t76), -t79, -g(3) * (t58 * pkin(6) + t60);];
U_reg = t1;
