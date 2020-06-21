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
% MDP [36x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = CloosQRC350DE_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(7,1),zeros(36,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [36 1]), ...
  'CloosQRC350DE_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [36x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:35
% EndTime: 2020-06-19 21:39:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (131->42), mult. (180->82), div. (0->0), fcn. (152->10), ass. (0->36)
unknown=NaN(6,1);
t1 = sin(qJ(2));
t4 = cos(qJ(2));
t7 = qJ(2) + qJ(3);
t8 = sin(t7);
t9 = g(3) * t8;
t10 = t9 * MDP(14);
t11 = cos(t7);
t12 = g(3) * t11;
t13 = t12 * MDP(15);
t14 = cos(qJ(4));
t16 = t9 * t14 * MDP(21);
t17 = sin(qJ(4));
t19 = t9 * t17 * MDP(22);
t20 = t8 * t14;
t21 = cos(qJ(5));
t23 = sin(qJ(5));
t25 = -t11 * t23 - t20 * t21;
t27 = g(3) * t25 * MDP(28);
t32 = g(3) * (-t11 * t21 + t20 * t23) * MDP(29);
t34 = pkin(7) * qJ(5) - qJ(6);
t35 = cos(t34);
t37 = t8 * t17;
t38 = sin(t34);
t42 = g(3) * (-t25 * t35 + t37 * t38) * MDP(35);
t47 = g(3) * (-t25 * t38 - t37 * t35) * MDP(36);
t60 = t11 * t17;
t63 = t11 * t14;
t77 = -t8 * t21 - t63 * t23;
t82 = -t63 * t21 + t8 * t23;
t86 = -t82 * pkin(7);
unknown(1,1) = 0;
unknown(2,1) = (g(3) * t1 * MDP(7) + g(3) * t4 * MDP(8) + t10 + t13 + t16 - t19 - t27 - t32 - t42 - t47);
unknown(3,1) = (t10 + t13 + t16 - t19 - t27 - t32 - t42 - t47);
unknown(4,1) = (t12 * t17 * MDP(21) + t12 * t14 * MDP(22) + t12 * t17 * t21 * MDP(28) - t12 * t17 * t23 * MDP(29) - g(3) * (t60 * t21 * t35 - t63 * t38) * MDP(35) - g(3) * (t60 * t21 * t38 + t63 * t35) * MDP(36));
unknown(5,1) = (-g(3) * t77 * MDP(28) - g(3) * t82 * MDP(29) - g(3) * (-t60 * pkin(7) * t35 - t77 * t35 + t86 * t38) * MDP(35) - g(3) * (-t60 * pkin(7) * t38 - t86 * t35 - t77 * t38) * MDP(36));
unknown(6,1) = (-g(3) * (t60 * t35 + t82 * t38) * MDP(35) - g(3) * (-t82 * t35 + t60 * t38) * MDP(36));
taug = unknown;
