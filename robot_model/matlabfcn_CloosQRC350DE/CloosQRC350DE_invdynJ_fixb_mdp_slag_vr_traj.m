% Inverse Dynamik für komplette Trajektorie für
% CloosQRC350DE
%
% Eingabe:
% RV_Traj [NTx63]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see CloosQRC350DE_invdynJ_fixb_regmin2vec.m
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350DE_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx6]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = CloosQRC350DE_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,63]), zeros(19,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 63), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx63] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350DE_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [19x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 6);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = CloosQRC350DE_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
