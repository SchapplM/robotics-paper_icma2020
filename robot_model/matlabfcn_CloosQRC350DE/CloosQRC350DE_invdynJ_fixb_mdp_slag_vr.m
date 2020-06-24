% Calculate vector of inverse dynamics joint torques for
% CloosQRC350DE
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [63x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see CloosQRC350DE_invdynJ_fixb_regmin2vec.m
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = CloosQRC350DE_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(63,1), zeros(19,1)}
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vr: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-23 21:07:42
% EndTime: 2020-06-23 21:07:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (57->57), mult. (63->63), div. (0->0), fcn. (63->63), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(4) * MDP(3) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(12) * MDP(7) + RV(15) * MDP(8) + RV(18) * MDP(9) + RV(23) * MDP(11) + RV(26) * MDP(12) + RV(29) * MDP(13) + RV(33) * MDP(14) + RV(37) * MDP(15) + RV(42) * MDP(16) + RV(47) * MDP(17) + RV(52) * MDP(18) + RV(58) * MDP(19); RV(3) * MDP(2) + RV(5) * MDP(3) + RV(6) * MDP(4) + RV(8) * MDP(5) + RV(10) * MDP(6) + RV(13) * MDP(7) + RV(16) * MDP(8) + RV(19) * MDP(9) + RV(21) * MDP(10) + RV(24) * MDP(11) + RV(27) * MDP(12) + RV(30) * MDP(13) + RV(34) * MDP(14) + RV(38) * MDP(15) + RV(43) * MDP(16) + RV(48) * MDP(17) + RV(53) * MDP(18) + RV(59) * MDP(19); RV(11) * MDP(6) + RV(14) * MDP(7) + RV(17) * MDP(8) + RV(20) * MDP(9) + RV(22) * MDP(10) + RV(25) * MDP(11) + RV(28) * MDP(12) + RV(31) * MDP(13) + RV(35) * MDP(14) + RV(39) * MDP(15) + RV(44) * MDP(16) + RV(49) * MDP(17) + RV(54) * MDP(18) + RV(60) * MDP(19); RV(32) * MDP(13) + RV(36) * MDP(14) + RV(40) * MDP(15) + RV(45) * MDP(16) + RV(50) * MDP(17) + RV(55) * MDP(18) + RV(61) * MDP(19); RV(41) * MDP(15) + RV(46) * MDP(16) + RV(51) * MDP(17) + RV(56) * MDP(18) + RV(62) * MDP(19); RV(57) * MDP(18) + RV(63) * MDP(19);];
tauJ = t1;
