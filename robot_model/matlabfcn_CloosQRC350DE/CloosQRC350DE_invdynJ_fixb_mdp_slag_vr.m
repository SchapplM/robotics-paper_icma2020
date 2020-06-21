% Calculate vector of inverse dynamics joint torques for
% CloosQRC350DE
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [139x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see CloosQRC350DE_invdynJ_fixb_regmin2vec.m
% MDP [36x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 21:40
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = CloosQRC350DE_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(139,1), zeros(36,1)}
assert(isreal(MDP) && all(size(MDP) == [36 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vr: MDP has to be [36x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 1
% StartTime: 2020-06-19 21:39:43
% EndTime: 2020-06-19 21:39:43
% DurationCPUTime: 0.27s
% Computational Cost: add. (133->133), mult. (139->139), div. (0->0), fcn. (139->139), ass. (0->153)
unknown=NaN(6,1);
t1 = RV(134);
t3 = RV(122);
t5 = RV(128);
t7 = RV(110);
t9 = RV(116);
t11 = RV(98);
t13 = RV(104);
t15 = RV(88);
t17 = RV(93);
t19 = RV(78);
t21 = RV(83);
t23 = RV(68);
t25 = RV(73);
t27 = RV(59);
t29 = RV(63);
t31 = RV(47);
t33 = RV(51);
t35 = t1 * MDP(36) + t3 * MDP(34) + t5 * MDP(35) + t7 * MDP(32) + t9 * MDP(33) + t11 * MDP(30) + t13 * MDP(31) + t15 * MDP(28) + t17 * MDP(29) + t19 * MDP(26) + t21 * MDP(27) + t23 * MDP(24) + t25 * MDP(25) + t27 * MDP(22) + t29 * MDP(23) + t31 * MDP(19) + t33 * MDP(20);
t36 = RV(55);
t38 = RV(35);
t40 = RV(39);
t42 = RV(43);
t44 = RV(24);
t46 = RV(29);
t48 = RV(32);
t50 = RV(11);
t52 = RV(13);
t54 = RV(15);
t56 = RV(18);
t58 = RV(21);
t60 = RV(1);
t62 = RV(2);
t64 = RV(4);
t66 = RV(6);
t68 = RV(8);
t70 = t36 * MDP(21) + t38 * MDP(16) + t40 * MDP(17) + t42 * MDP(18) + t44 * MDP(12) + t46 * MDP(14) + t48 * MDP(15) + t50 * MDP(7) + t52 * MDP(8) + t54 * MDP(9) + t56 * MDP(10) + t58 * MDP(11) + t60 * MDP(1) + t62 * MDP(2) + t64 * MDP(3) + t66 * MDP(4) + t68 * MDP(5);
t72 = RV(135);
t74 = RV(123);
t76 = RV(129);
t78 = RV(111);
t80 = RV(117);
t82 = RV(99);
t84 = RV(105);
t86 = RV(89);
t88 = RV(94);
t90 = RV(74);
t92 = RV(79);
t94 = RV(84);
t96 = RV(69);
t98 = RV(56);
t100 = RV(60);
t102 = RV(64);
t104 = RV(48);
t106 = t72 * MDP(36) + t74 * MDP(34) + t76 * MDP(35) + t78 * MDP(32) + t80 * MDP(33) + t82 * MDP(30) + t84 * MDP(31) + t86 * MDP(28) + t88 * MDP(29) + t90 * MDP(25) + t92 * MDP(26) + t94 * MDP(27) + t96 * MDP(24) + t98 * MDP(21) + t100 * MDP(22) + t102 * MDP(23) + t104 * MDP(19);
t107 = RV(52);
t109 = RV(36);
t111 = RV(40);
t113 = RV(44);
t115 = RV(22);
t117 = RV(25);
t119 = RV(27);
t121 = RV(30);
t123 = RV(33);
t125 = RV(9);
t127 = RV(10);
t129 = RV(12);
t131 = RV(14);
t133 = RV(16);
t135 = RV(19);
t137 = RV(3);
t139 = RV(5);
t141 = RV(7);
t143 = t107 * MDP(20) + t109 * MDP(16) + t111 * MDP(17) + t113 * MDP(18) + t115 * MDP(11) + t117 * MDP(12) + t119 * MDP(13) + t121 * MDP(14) + t123 * MDP(15) + t125 * MDP(5) + t127 * MDP(6) + t129 * MDP(7) + t131 * MDP(8) + t133 * MDP(9) + t135 * MDP(10) + t137 * MDP(2) + t139 * MDP(3) + t141 * MDP(4);
t145 = RV(136);
t147 = RV(124);
t149 = RV(130);
t151 = RV(112);
t153 = RV(118);
t155 = RV(100);
t157 = RV(106);
t159 = RV(85);
t161 = RV(90);
t163 = RV(95);
t165 = RV(75);
t167 = RV(80);
t169 = RV(70);
t171 = RV(57);
t173 = MDP(21) * t171 + MDP(24) * t169 + MDP(25) * t165 + MDP(26) * t167 + MDP(27) * t159 + MDP(28) * t161 + MDP(29) * t163 + MDP(30) * t155 + MDP(31) * t157 + MDP(32) * t151 + MDP(33) * t153 + MDP(34) * t147 + MDP(35) * t149 + MDP(36) * t145;
t174 = RV(61);
t176 = RV(65);
t178 = RV(45);
t180 = RV(49);
t182 = RV(53);
t184 = RV(34);
t186 = RV(37);
t188 = RV(41);
t190 = RV(23);
t192 = RV(26);
t194 = RV(28);
t196 = RV(31);
t198 = RV(17);
t200 = RV(20);
t202 = MDP(10) * t200 + MDP(11) * t190 + MDP(12) * t192 + MDP(13) * t194 + MDP(14) * t196 + MDP(15) * t184 + MDP(16) * t186 + MDP(17) * t188 + MDP(18) * t178 + MDP(19) * t180 + MDP(20) * t182 + MDP(22) * t174 + MDP(23) * t176 + MDP(9) * t198;
t204 = RV(38);
t206 = RV(42);
t208 = RV(46);
t210 = RV(50);
t212 = RV(54);
t214 = RV(58);
t216 = RV(62);
t218 = RV(66);
t220 = RV(71);
t222 = RV(76);
t225 = RV(81);
t227 = RV(86);
t229 = RV(91);
t231 = RV(96);
t233 = RV(101);
t235 = RV(107);
t237 = RV(113);
t239 = RV(119);
t241 = RV(125);
t243 = RV(131);
t245 = RV(137);
t247 = MDP(26) * t225 + MDP(27) * t227 + MDP(28) * t229 + MDP(29) * t231 + MDP(30) * t233 + MDP(31) * t235 + MDP(32) * t237 + MDP(33) * t239 + MDP(34) * t241 + MDP(35) * t243 + MDP(36) * t245;
t249 = RV(67);
t251 = RV(72);
t253 = RV(77);
t255 = RV(82);
t257 = RV(87);
t259 = RV(92);
t261 = RV(97);
t263 = RV(102);
t265 = RV(108);
t267 = RV(114);
t269 = RV(120);
t271 = RV(126);
t273 = RV(132);
t275 = RV(138);
t277 = MDP(23) * t249 + MDP(24) * t251 + MDP(25) * t253 + MDP(26) * t255 + MDP(27) * t257 + MDP(28) * t259 + MDP(29) * t261 + MDP(30) * t263 + MDP(31) * t265 + MDP(32) * t267 + MDP(33) * t269 + MDP(34) * t271 + MDP(35) * t273 + MDP(36) * t275;
t278 = RV(103);
t280 = RV(109);
t282 = RV(115);
t284 = RV(121);
t286 = RV(127);
t288 = RV(133);
t290 = RV(139);
unknown(1,1) = t35 + t70;
unknown(2,1) = t106 + t143;
unknown(3,1) = t173 + t202;
unknown(4,1) = MDP(16) * t204 + MDP(17) * t206 + MDP(18) * t208 + MDP(19) * t210 + MDP(20) * t212 + MDP(21) * t214 + MDP(22) * t216 + MDP(23) * t218 + MDP(24) * t220 + MDP(25) * t222 + t247;
unknown(5,1) = t277;
unknown(6,1) = MDP(30) * t278 + MDP(31) * t280 + MDP(32) * t282 + MDP(33) * t284 + MDP(34) * t286 + MDP(35) * t288 + MDP(36) * t290;
tauJ = unknown;
