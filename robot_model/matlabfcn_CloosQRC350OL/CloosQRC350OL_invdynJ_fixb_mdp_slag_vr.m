% Calculate vector of inverse dynamics joint torques for
% CloosQRC350OL
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [139x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see CloosQRC350OL_invdynJ_fixb_regmin2vec.m
% MDP [36x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = CloosQRC350OL_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(139,1), zeros(36,1)}
assert(isreal(MDP) && all(size(MDP) == [36 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vr: MDP has to be [36x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-20 08:20:07
% EndTime: 2020-06-20 08:20:07
% DurationCPUTime: 0.13s
% Computational Cost: add. (133->133), mult. (139->139), div. (0->0), fcn. (139->139), ass. (0->1)
t1 = [RV(73) * MDP(25) + RV(78) * MDP(26) + RV(83) * MDP(27) + RV(128) * MDP(35) + RV(134) * MDP(36) + RV(47) * MDP(19) + RV(51) * MDP(20) + RV(55) * MDP(21) + RV(59) * MDP(22) + RV(88) * MDP(28) + RV(93) * MDP(29) + RV(63) * MDP(23) + RV(68) * MDP(24) + RV(116) * MDP(33) + RV(122) * MDP(34) + RV(1) * MDP(1) + RV(2) * MDP(2) + RV(4) * MDP(3) + RV(6) * MDP(4) + RV(8) * MDP(5) + RV(11) * MDP(7) + RV(13) * MDP(8) + RV(15) * MDP(9) + RV(18) * MDP(10) + RV(35) * MDP(16) + RV(39) * MDP(17) + RV(43) * MDP(18) + RV(21) * MDP(11) + RV(24) * MDP(12) + RV(29) * MDP(14) + RV(32) * MDP(15) + RV(98) * MDP(30) + RV(104) * MDP(31) + RV(110) * MDP(32); RV(74) * MDP(25) + RV(79) * MDP(26) + RV(84) * MDP(27) + RV(123) * MDP(34) + RV(129) * MDP(35) + RV(135) * MDP(36) + RV(48) * MDP(19) + RV(52) * MDP(20) + RV(56) * MDP(21) + RV(89) * MDP(28) + RV(94) * MDP(29) + RV(60) * MDP(22) + RV(64) * MDP(23) + RV(69) * MDP(24) + RV(111) * MDP(32) + RV(117) * MDP(33) + RV(3) * MDP(2) + RV(5) * MDP(3) + RV(7) * MDP(4) + RV(9) * MDP(5) + RV(10) * MDP(6) + RV(12) * MDP(7) + RV(14) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10) + RV(36) * MDP(16) + RV(40) * MDP(17) + RV(44) * MDP(18) + RV(22) * MDP(11) + RV(25) * MDP(12) + RV(27) * MDP(13) + RV(30) * MDP(14) + RV(33) * MDP(15) + RV(99) * MDP(30) + RV(105) * MDP(31); RV(75) * MDP(25) + RV(80) * MDP(26) + RV(124) * MDP(34) + RV(130) * MDP(35) + RV(49) * MDP(19) + RV(53) * MDP(20) + RV(57) * MDP(21) + RV(85) * MDP(27) + RV(90) * MDP(28) + RV(95) * MDP(29) + RV(61) * MDP(22) + RV(65) * MDP(23) + RV(70) * MDP(24) + RV(112) * MDP(32) + RV(118) * MDP(33) + RV(136) * MDP(36) + RV(17) * MDP(9) + RV(20) * MDP(10) + RV(34) * MDP(15) + RV(37) * MDP(16) + RV(41) * MDP(17) + RV(45) * MDP(18) + RV(23) * MDP(11) + RV(26) * MDP(12) + RV(28) * MDP(13) + RV(31) * MDP(14) + RV(100) * MDP(30) + RV(106) * MDP(31); RV(38) * MDP(16) + RV(42) * MDP(17) + RV(46) * MDP(18) + RV(50) * MDP(19) + RV(54) * MDP(20) + RV(58) * MDP(21) + RV(62) * MDP(22) + RV(66) * MDP(23) + RV(71) * MDP(24) + RV(76) * MDP(25) + RV(81) * MDP(26) + RV(86) * MDP(27) + RV(91) * MDP(28) + RV(96) * MDP(29) + RV(101) * MDP(30) + RV(107) * MDP(31) + RV(113) * MDP(32) + RV(119) * MDP(33) + RV(125) * MDP(34) + RV(131) * MDP(35) + RV(137) * MDP(36); RV(67) * MDP(23) + RV(72) * MDP(24) + RV(77) * MDP(25) + RV(82) * MDP(26) + RV(87) * MDP(27) + RV(92) * MDP(28) + RV(97) * MDP(29) + RV(102) * MDP(30) + RV(108) * MDP(31) + RV(114) * MDP(32) + RV(120) * MDP(33) + RV(126) * MDP(34) + RV(132) * MDP(35) + RV(138) * MDP(36); RV(103) * MDP(30) + RV(109) * MDP(31) + RV(115) * MDP(32) + RV(121) * MDP(33) + RV(127) * MDP(34) + RV(133) * MDP(35) + RV(139) * MDP(36);];
tauJ = t1;
