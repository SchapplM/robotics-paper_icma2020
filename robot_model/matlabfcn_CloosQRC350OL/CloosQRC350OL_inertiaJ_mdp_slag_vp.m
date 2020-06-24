% Calculate joint inertia matrix for
% CloosQRC350OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350OL_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350OL_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 22:03:20
% EndTime: 2020-06-23 22:03:21
% DurationCPUTime: 0.41s
% Computational Cost: add. (213->76), mult. (507->129), div. (0->0), fcn. (534->10), ass. (0->48)
t90 = sin(qJ(2));
t76 = t90 * pkin(3) + pkin(2);
t114 = 0.2e1 * t76;
t113 = 2 * MDP(17);
t89 = sin(qJ(3));
t94 = cos(qJ(3));
t95 = cos(qJ(2));
t68 = -t89 * t90 + t94 * t95;
t88 = sin(qJ(4));
t112 = t68 * t88;
t93 = cos(qJ(4));
t111 = t68 * t93;
t83 = t88 ^ 2;
t92 = cos(qJ(5));
t110 = t83 * t92;
t109 = t92 * t93;
t86 = sin(qJ(6));
t106 = MDP(18) * t86;
t91 = cos(qJ(6));
t62 = t91 * t92 * t88 + t86 * t93;
t79 = t93 * MDP(16);
t108 = t62 * t106 + t79;
t66 = t89 * t95 + t94 * t90;
t107 = MDP(17) * (t66 * pkin(4) + t68 * pkin(5) + t76);
t105 = MDP(18) * t91;
t104 = t66 * MDP(14);
t103 = t92 * MDP(15);
t102 = t92 * MDP(19);
t75 = t89 * pkin(3) - pkin(5);
t101 = t75 * t110;
t100 = t62 * t105;
t87 = sin(qJ(5));
t72 = t87 * t88 * MDP(19);
t56 = t68 * t109 - t87 * t66;
t53 = t86 * t112 - t91 * t56;
t55 = -t87 * t111 - t92 * t66;
t99 = t53 * t62 * MDP(18) + t68 * MDP(8) - t66 * MDP(9) + t79 * t112 + t55 * t72;
t82 = t87 ^ 2;
t84 = t92 ^ 2;
t85 = t93 ^ 2;
t98 = t62 ^ 2 * MDP(18) + t85 * MDP(16) + MDP(10) + (MDP(15) * t84 + MDP(19) * t82 + MDP(13)) * t83;
t97 = (t94 * MDP(11) - t89 * MDP(12)) * pkin(3);
t96 = -MDP(13) * t111 - t56 * t103;
t74 = pkin(5) * t110;
t70 = t92 * t72;
t69 = -t87 * pkin(4) + pkin(5) * t109;
t58 = t75 * t109 + t87 * (t94 * pkin(3) + pkin(4));
t1 = [0.2e1 * pkin(2) * t90 * MDP(5) + t68 * MDP(12) * t114 + t56 ^ 2 * MDP(15) + t53 ^ 2 * MDP(18) + t55 ^ 2 * MDP(19) + t95 ^ 2 * MDP(2) + MDP(1) + 0.2e1 * (t68 * t110 + t56 * t93) * t107 + (MDP(11) * t114 - 0.2e1 * t68 * MDP(7) + t104) * t66 + (t85 * MDP(13) + t83 * MDP(16) + MDP(6)) * t68 ^ 2; t95 * MDP(3) + ((t56 * t75 - t58 * t68) * MDP(17) + t96) * t88 + t99; MDP(4) + (-t58 * t93 - t101) * t113 + 0.2e1 * t97 + t98; ((-pkin(5) * t56 + t68 * t69) * MDP(17) + t96) * t88 + t99; (-t101 + t74 + (-t58 + t69) * t93) * MDP(17) + t97 + t98; (t69 * t93 + t74) * t113 + t98; t55 * t102 - t104 + (MDP(15) * t56 - t53 * t105 + t93 * t107) * t87; t70 + (-t100 + (t75 * MDP(17) - t103) * t88) * t87; t70 + (-t100 + (-pkin(5) * MDP(17) - t103) * t88) * t87; t84 * MDP(19) + MDP(14) + (MDP(18) * t91 ^ 2 + MDP(15)) * t82; t53 * t106 + (MDP(16) * t68 + t92 * t107) * t88; -t58 * MDP(17) + t108; t69 * MDP(17) + t108; -t87 * t86 * t105; t86 ^ 2 * MDP(18) + MDP(16); t55 * MDP(19); t72; t72; t102; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
