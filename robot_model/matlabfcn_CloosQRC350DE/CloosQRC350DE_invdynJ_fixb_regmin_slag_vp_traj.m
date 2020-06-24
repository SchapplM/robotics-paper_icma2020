% Calculate time series of minimal parameter regressor of inv. dyn. joint torques for
% CloosQRC350DE
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,kDG]';
%
% Output:
% RV_Traj [NTx63]
%   time series of regressor matrices as vectors
%   see CloosQRC350DE_invdynJ_fixb_regmin2vec.m

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV_Traj = CloosQRC350DE_invdynJ_fixb_regmin_slag_vp_traj(Q, QD, QDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),
%$cgargs  zeros(3,1), zeros(7,1)}
assert(isreal(Q) && all(size(Q,2) == 6), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp_traj: Q needs to be [NTx6] (double)');
assert(isreal(QD) && all(size(QD,2) == 6), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp_traj: QD needs to be [NTx6] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 6), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp_traj: QDD needs to be [NTx6] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'CloosQRC350DE_invdynJ_fixb_regmin_slag_vp_traj: Kinematic parameters pkin have to be [7x1] (double)');
  
%% Trajektorie der Regressor-Vektoren aufbauen
RV_Traj = NaN(size(Q,1), 63);
for ii = 1:size(Q,1)
  RV_Traj(ii,:) = CloosQRC350DE_invdynJ_fixb_regmin2vec( ...
    CloosQRC350DE_invdynJ_fixb_regmin_slag_vp(Q(ii,:)', QD(ii,:)', QDD(ii,:)', g, pkin) );
end
