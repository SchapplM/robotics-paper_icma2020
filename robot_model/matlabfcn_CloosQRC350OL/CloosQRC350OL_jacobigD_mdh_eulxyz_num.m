% Jacobimatrix-Zeitableitung zu beliebigen Punkten eines Körpers für
% CloosQRC350OL
% Floating-Base-Notation: Jacobi-Matrix bezogen auf Basis und Gelenke
%
% Input:
% phi_base [3x1]
%   Base orientation in world frame. Expressed with XYZ-Euler angles
% phiD_base [3x1]
%   Time Derivative of Base Orientation in world frame.
%   Expressed with XYZ-Euler angles ("rpy")
% qJ [6x1]
%   Joint Angles [rad]
% qJD [6x1]
%   Joint Velocities [rad/s]
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: CloosQRC350OL_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
%
% Output:
% JgD_C [6x(6+6)]
%   time derivative of geometric body jacobian for the defined point with respect to
%   generalized coordinates:
%
% Quellen:
% [1] Aufzeichnungen Schappler (vom 22.06.2016)
% [2] BouyarmaneKhe2012: On the dynamics modeling of free-floating-base
%     articulated mechanisms and applications to humanoid whole-body dynamics
%     and control

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-20 08:27
% Revision: 6013df02bda2c1f6ebc95d3649839f696d960e41 (2020-06-19)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-06
% (C) Institut für Regelungstechnik, Leibniz Universität Hannover

function JgD_C = CloosQRC350OL_jacobigD_mdh_eulxyz_num(phi_base, phiD_base, qJ, qJD, link_index, r_i_i_C, ...
  pkin)
%% Init
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),
%$cgargs  uint8(zeros(1,1)),zeros(3,1),zeros(6,1)}
assert(isreal(phi_base) && all(size(phi_base) == [3 1]), ...
  'CloosQRC350OL_jacobigD_mdh_eulxyz_num: Base Euler angles have to be [3x1] (double)');
assert(isreal(phiD_base) && all(size(phiD_base) == [3 1]), ...
  'CloosQRC350OL_jacobigD_mdh_eulxyz_num: Base Euler angles time derivatives have to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobigD_mdh_eulxyz_num: Joint angles qJ have to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'CloosQRC350OL_jacobigD_mdh_eulxyz_num: Joint velocities qJD have to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
  'CloosQRC350OL_jacobigD_mdh_eulxyz_num: link_index has to be [1x1] uint8');
assert(isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
  'CloosQRC350OL_jacobigD_mdh_eulxyz_num: Position vector r_i_i_C has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobigD_mdh_eulxyz_num: Kinematic parameters pkin have to be [6x1] (double)');

if link_index > 7-1
  error('Index exceeds number of bodies');
end

%% Kinematik berechnen
T_c_mdh = CloosQRC350OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin);
v_mdh = CloosQRC350OL_structural_kinematic_parameters();
T_B_i = T_c_mdh(:,:,link_index+1);
% Vektor von der Basis (B) zum Punkt C (auf Körper i)
r_B_B_C = T_B_i(1:3,4) + T_B_i(1:3,1:3) * r_i_i_C; % r_B_B_i + R_B_i_C

% Rotationsmatrix vom Welt- ins Basis-KS
R_W_B = eulxyz2r(phi_base);

% Transformation zwischen Euler-Winkel-Zeitableitung und
% Winkelgeschwindigkeit der Basis im Welt-KS
T_angvel = eulxyzjac(phi_base);
TD_angvel = eulxyzjacD(phi_base, phiD_base);

omega_W_B = T_angvel*phiD_base;

%% Gelenk-Jacobi-Matrix
Jg_JC_B = CloosQRC350OL_jacobig_mdh_num(qJ, link_index, r_i_i_C, pkin);
JgD_JC_B = CloosQRC350OL_jacobigD_mdh_num(qJ, qJD, link_index, r_i_i_C, pkin);
% Rotation ins Welt-KS und Zeitableitung der Rotation
JgD_JTC = skew(omega_W_B)*R_W_B*Jg_JC_B(1:3,:) + R_W_B*JgD_JC_B(1:3,:);
JgD_JRC = skew(omega_W_B)*R_W_B*Jg_JC_B(4:6,:) + R_W_B*JgD_JC_B(4:6,:);

%% Basis-Jacobi-Matrix
% Siehe [1]
% [2], equ. (14). Zeitableitung durch Produktregel und Eulersche
% Differentiationsregel
Term5_Teil1 = skew(Jg_JC_B(1:3,:)*qJD) * (R_W_B'); % erster Term aus [1] Gl. (5)
Term5_Teil2 = skew(r_B_B_C)*((skew(omega_W_B)*R_W_B)'); % zweiter Term aus [1] Gl. (5)
Term5 = Term5_Teil1 + Term5_Teil2;
Term4_Teil1 = skew(omega_W_B) * R_W_B * skew(r_B_B_C) * (R_W_B'); % erster Term aus [1] Gl. (4)
Term4 = Term4_Teil1 + R_W_B*Term5;
Term3_Teil2 = R_W_B*skew(r_B_B_C)*R_W_B'*TD_angvel; % zweiter Term aus [1] Gl. (3)
Term3 = Term4*T_angvel + Term3_Teil2;

JgD_BTC = [zeros(3,3),     -Term3];
% [2], equ. (17). Zeitableitung
JgD_BRC = [zeros(3,3), TD_angvel];

%% Gesamt-Jacobi-Matrix zusammensetzen
JgD_C = [JgD_BTC, JgD_JTC; ... % [2], equ. (14)
         JgD_BRC, JgD_JRC]; % [2], equ. (17)
