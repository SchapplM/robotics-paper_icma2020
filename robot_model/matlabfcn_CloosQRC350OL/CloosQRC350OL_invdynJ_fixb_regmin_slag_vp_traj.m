% Calculate time series of minimal parameter regressor of inv. dyn. joint torques for
% CloosQRC350OL
%
% Input:
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
%
% Output:
% RV_Traj [NTx139]
%   time series of regressor matrices as vectors
%   see CloosQRC350OL_invdynJ_fixb_regmin2vec.m

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV_Traj = CloosQRC350OL_invdynJ_fixb_regmin_slag_vp_traj(Q, QD, QDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),
%$cgargs  zeros(3,1), zeros(6,1)}
assert(isreal(Q) && all(size(Q,2) == 6), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp_traj: Q needs to be [NTx6] (double)');
assert(isreal(QD) && all(size(QD,2) == 6), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp_traj: QD needs to be [NTx6] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 6), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp_traj: QDD needs to be [NTx6] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_invdynJ_fixb_regmin_slag_vp_traj: Kinematic parameters pkin have to be [6x1] (double)');
  
%% Trajektorie der Regressor-Vektoren aufbauen
RV_Traj = NaN(size(Q,1), 139);
for ii = 1:size(Q,1)
  RV_Traj(ii,:) = CloosQRC350OL_invdynJ_fixb_regmin2vec( ...
    CloosQRC350OL_invdynJ_fixb_regmin_slag_vp(Q(ii,:)', QD(ii,:)', QDD(ii,:)', g, pkin) );
end
