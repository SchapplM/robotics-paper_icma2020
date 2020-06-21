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
% MDP [36x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350OL_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(36,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [36 1]), ...
  'CloosQRC350OL_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [36x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:18:58
% EndTime: 2020-06-20 08:18:59
% DurationCPUTime: 0.15s
% Computational Cost: add. (97->28), mult. (142->50), div. (0->0), fcn. (134->10), ass. (0->20)
t44 = qJ(2) + qJ(3);
t43 = cos(t44);
t58 = g(3) * t43;
t45 = sin(qJ(6));
t47 = sin(qJ(4));
t57 = t45 * t47;
t46 = sin(qJ(5));
t50 = cos(qJ(4));
t56 = t46 * t50;
t48 = cos(qJ(6));
t55 = t47 * t48;
t49 = cos(qJ(5));
t54 = t49 * t50;
t42 = sin(t44);
t53 = MDP(15) * t58 + (MDP(21) * t50 + MDP(14)) * g(3) * t42;
t52 = MDP(35) * t48 - MDP(36) * t45 - MDP(28);
t36 = -t42 * t54 - t43 * t46;
t51 = -MDP(22) * t42 * t47 - MDP(28) * t36 - MDP(29) * (t42 * t56 - t43 * t49) - MDP(35) * (-t36 * t48 - t42 * t57) - MDP(36) * (t36 * t45 - t42 * t55);
t38 = -t42 * t46 + t43 * t54;
t1 = [0; (MDP(7) * sin(qJ(2)) + MDP(8) * cos(qJ(2)) + t51) * g(3) + t53; t51 * g(3) + t53; ((-t45 * MDP(35) - t48 * MDP(36) + MDP(22)) * t50 + (-t46 * MDP(29) - t52 * t49 + MDP(21)) * t47) * t58; (MDP(29) * t38 + t52 * (-t42 * t49 - t43 * t56)) * g(3); (-(t38 * t45 + t43 * t55) * MDP(35) - (t38 * t48 - t43 * t57) * MDP(36)) * g(3);];
taug = t1;
