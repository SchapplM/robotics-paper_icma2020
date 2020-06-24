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
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
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
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:29
% EndTime: 2020-06-23 21:07:29
% DurationCPUTime: 0.04s
% Computational Cost: add. (6->4), mult. (15->11), div. (0->0), fcn. (20->7), ass. (0->7)
t50 = cos(qJ(2));
t49 = cos(qJ(3));
t48 = sin(qJ(2));
t47 = sin(qJ(3));
t46 = t50 * t47 + t49 * t48;
t45 = t48 * t47 - t50 * t49;
t1 = [0, 0, 0, 0, -t50 * g(3), 0, 0, 0, 0, 0, t45 * g(3), t46 * g(3), 0, 0, 0, 0, (-cos(qJ(4)) * t45 * sin(qJ(5)) + cos(qJ(5)) * t46) * g(3), 0, 0;];
U_reg = t1;
