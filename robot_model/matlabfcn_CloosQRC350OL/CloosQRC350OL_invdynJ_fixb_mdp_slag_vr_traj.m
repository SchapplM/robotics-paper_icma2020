% Inverse Dynamik für komplette Trajektorie für
% CloosQRC350OL
%
% Eingabe:
% RV_Traj [NTx139]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see CloosQRC350OL_invdynJ_fixb_regmin2vec.m
% MDP [36x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx6]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = CloosQRC350OL_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,139]), zeros(36,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 139), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx139] (double)');
assert(isreal(MDP) && all(size(MDP) == [36 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [36x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 6);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = CloosQRC350OL_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
