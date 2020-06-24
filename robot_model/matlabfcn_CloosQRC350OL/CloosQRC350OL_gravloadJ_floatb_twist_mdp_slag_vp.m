% Calculate Gravitation load on the joints for
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350OL_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:19
% EndTime: 2020-06-23 22:03:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (21->9), mult. (29->16), div. (0->0), fcn. (23->7), ass. (0->11)
t18 = qJ(2) + qJ(3);
t16 = sin(t18);
t19 = sin(qJ(5));
t26 = t16 * t19;
t17 = cos(t18);
t25 = t17 * cos(qJ(5));
t24 = (MDP(11) * t16 + MDP(12) * t17) * g(3);
t23 = MDP(17) * g(3);
t21 = cos(qJ(4));
t22 = MDP(17) * (t21 * t26 - t25);
t1 = [0; (-t22 + MDP(5) * sin(qJ(2))) * g(3) + t24; -g(3) * t22 + t24; -t17 * sin(qJ(4)) * t19 * t23; -(-t21 * t25 + t26) * t23; 0;];
taug = t1;
