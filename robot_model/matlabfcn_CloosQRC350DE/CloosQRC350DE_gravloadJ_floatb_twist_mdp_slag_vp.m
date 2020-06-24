% Calculate Gravitation load on the joints for
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350DE_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:36
% EndTime: 2020-06-23 21:07:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (19->7), mult. (51->19), div. (0->0), fcn. (56->8), ass. (0->13)
t32 = cos(qJ(2));
t31 = cos(qJ(3));
t23 = sin(qJ(3));
t24 = sin(qJ(2));
t20 = t32 * t23 + t31 * t24;
t22 = sin(qJ(5));
t30 = t22 * t20;
t19 = t24 * t23 - t32 * t31;
t29 = cos(qJ(5)) * t19;
t28 = MDP(17) * g(3);
t26 = cos(qJ(4));
t27 = (-t26 * t30 - t29) * t28 + (t20 * MDP(11) - t19 * MDP(12)) * g(3);
t1 = [0; t24 * g(3) * MDP(5) + t27; t27; sin(qJ(4)) * t19 * t22 * t28; (-t26 * t29 - t30) * t28; 0;];
taug = t1;
