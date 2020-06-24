% Calculate joint inertia matrix for
% CloosQRC350DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = CloosQRC350DE_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350DE_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350DE_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:38
% EndTime: 2020-06-23 21:07:39
% DurationCPUTime: 0.38s
% Computational Cost: add. (259->81), mult. (565->136), div. (0->0), fcn. (550->10), ass. (0->50)
t88 = sin(qJ(3));
t89 = sin(qJ(2));
t92 = cos(qJ(3));
t93 = cos(qJ(2));
t65 = t89 * t88 - t93 * t92;
t113 = 0.2e1 * t65;
t112 = 2 * MDP(17);
t87 = sin(qJ(4));
t111 = t65 * t87;
t91 = cos(qJ(4));
t110 = t65 * t91;
t83 = t87 ^ 2;
t90 = cos(qJ(5));
t109 = t83 * t90;
t108 = t90 * t91;
t66 = t93 * t88 + t92 * t89;
t75 = -t89 * pkin(3) - pkin(2);
t107 = MDP(17) * (-t66 * pkin(4) + t65 * pkin(5) + t75);
t77 = pkin(7) * qJ(5) - qJ(6);
t73 = sin(t77);
t106 = MDP(18) * t73;
t74 = cos(t77);
t105 = MDP(18) * t74;
t86 = sin(qJ(5));
t56 = -t86 * t110 + t90 * t66;
t104 = t56 * MDP(19);
t103 = t66 * MDP(14);
t102 = t90 * MDP(15);
t101 = t90 * MDP(19);
t80 = t91 * MDP(16);
t76 = pkin(3) * t88 - pkin(5);
t100 = t76 * t109;
t59 = t74 * t90 * t87 - t73 * t91;
t99 = t59 * t105;
t70 = t86 * t87 * MDP(19);
t55 = t65 * t108 + t86 * t66;
t53 = -t73 * t111 - t74 * t55;
t98 = t53 * t59 * MDP(18) + t65 * MDP(8) + t66 * MDP(9) + t80 * t111 + t56 * t70;
t82 = t86 ^ 2;
t84 = t90 ^ 2;
t85 = t91 ^ 2;
t97 = t59 ^ 2 * MDP(18) + t85 * MDP(16) + MDP(10) + (MDP(15) * t84 + MDP(19) * t82 + MDP(13)) * t83;
t96 = (t92 * MDP(11) - t88 * MDP(12)) * pkin(3);
t95 = -MDP(13) * t110 - t55 * t102;
t94 = -pkin(7) * t70 - t59 * t106 + t80;
t72 = pkin(5) * t109;
t68 = t90 * t70;
t67 = -t86 * pkin(4) + pkin(5) * t108;
t60 = t76 * t108 + t86 * (pkin(3) * t92 + pkin(4));
t1 = [0.2e1 * pkin(2) * t89 * MDP(5) + t75 * MDP(12) * t113 + t55 ^ 2 * MDP(15) + t53 ^ 2 * MDP(18) + t56 ^ 2 * MDP(19) + t93 ^ 2 * MDP(2) + MDP(1) + 0.2e1 * (t65 * t109 + t55 * t91) * t107 + (-0.2e1 * t75 * MDP(11) + MDP(7) * t113 + t103) * t66 + (t85 * MDP(13) + t83 * MDP(16) + MDP(6)) * t65 ^ 2; -t93 * MDP(3) + ((t55 * t76 - t60 * t65) * MDP(17) + t95) * t87 + t98; MDP(4) + (-t60 * t91 - t100) * t112 + 0.2e1 * t96 + t97; ((-pkin(5) * t55 + t65 * t67) * MDP(17) + t95) * t87 + t98; (-t100 + t72 + (-t60 + t67) * t91) * MDP(17) + t96 + t97; (t67 * t91 + t72) * t112 + t97; t56 * t101 + t103 + (MDP(15) * t55 - t53 * t105 + t91 * t107) * t86; t68 + (-t99 + (t76 * MDP(17) - t102) * t87) * t86; t68 + (-t99 + (-pkin(5) * MDP(17) - t102) * t87) * t86; t84 * MDP(19) + MDP(14) + (MDP(18) * t74 ^ 2 + MDP(15)) * t82; -t53 * t106 - pkin(7) * t104 + (MDP(16) * t65 + t90 * t107) * t87; -t60 * MDP(17) + t94; t67 * MDP(17) + t94; t86 * t73 * t105 - pkin(7) * t101; t73 ^ 2 * MDP(18) + pkin(7) ^ 2 * MDP(19) + MDP(16); t104; t70; t70; t101; -pkin(7) * MDP(19); MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
