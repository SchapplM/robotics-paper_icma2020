% Inverse Dynamik für komplette Trajektorie für
% CloosQRC350OL
%
% Eingabe:
% Q [NTx6]
%   Trajektorie von Gelenkpositionen (NT Zeitschritte in den Zeilen)
% QD [NTx6]
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD [NTx6]
%   Trajektorie von Gelenkbeschleunigungen
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see CloosQRC350OL_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx6]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 22:05
% Revision: 9ee7546dde8543a81bf40e37a1400ef9d9e232c4 (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = CloosQRC350OL_invdynJ_fixb_mdp_slag_vp_traj(Q, QD, QDD, g, pkin, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),
%$cgargs  zeros(3,1), zeros(6,1), zeros(19,1)}
assert(isreal(Q) && all(size(Q,2) == 6), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp_traj: Q needs to be [NTx6] (double)');
assert(isreal(QD) && all(size(QD,2) == 6), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp_traj: QD needs to be [NTx6] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 6), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp_traj: QDD needs to be [NTx6] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp_traj: Kinematic parameters pkin have to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'CloosQRC350OL_invdynJ_fixb_mdp_slag_vp_traj: Dynamics parameter vector MDP has to be [19x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(Q));
for k = 1:size(Q,1)
  TAU(k,:) = CloosQRC350OL_invdynJ_fixb_mdp_slag_vp(Q(k,:)', QD(k,:)', QDD(k,:)', g, pkin, MDP);
end
