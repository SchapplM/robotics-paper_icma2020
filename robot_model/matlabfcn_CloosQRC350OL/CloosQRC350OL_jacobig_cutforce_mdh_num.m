% Geometrische Jacobi-Matrix für vollständige Schnittkräfte für einen
% beliebigen Roboter in kinematischer Baumstruktur
% Suffix "_m": Modulare Funktion (keine Übergabe von Gelenkwinkeln, sondern
% von Transformationsmatrizen. Daher für beliebige Roboter verwendbar
% unabhängig von der Gelenk-Transformations-Notation
% 
% Die Matrix kann verwendet werden, um die aus einer externen Kraft
% resultierenden Schnittkräfte in allen Gelenken des Roboters zu bestimmen
% 
% Eingabe:
% qJ [6x1]
%   Gelenkkoordinaten
% v [NJx1] uint8
%   Vorgänger-Indizes der Gelenke
% link_index [1x1] uint8
%   Index des Körpers, auf dem der Punkt C liegt. 0=Basis
% r_i_i_C [3x1]
%   Punkt C, zu dem die Jacobi-Matrix berechnet wird (Angriffspunkt der
%   Kraft)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% 
% Ausgabe:
% Jg [6x(6*NL)]
%   Geometrische Jacobimatrix für vollständige Schnittkräfte in allen
%   Körper-KS basierend auf externer Kraft am gegebenen Punkt
%
% Quellen:
% [1] Ortmaier: Robotik I Skript

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-05
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_C = CloosQRC350OL_jacobig_cutforce_mdh_num(qJ, link_index, r_i_i_C, pkin)

%% Init
%#codegen
%$cgargs {zeros(6,1),uint8(zeros(1,1)),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'CloosQRC350OL_jacobig_cutforce_mdh_num: Joint angles qJ have to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
  'robot_tree_jacobig_m: link_index has to be [1x1] uint8');
assert(isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
  'CloosQRC350OL_jacobig_cutforce_mdh_num: Position vector r_i_i_C has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'CloosQRC350OL_jacobig_cutforce_mdh_num: Kinematic parameters pkin have to be [6x1] (double)');

T_c_mdh = CloosQRC350OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin);
[v,~,~,NL] = CloosQRC350OL_structural_kinematic_parameters();

% Initialisierung. Alle Spalten die nicht gesetzt werden haben keinen
% Einfluss.
Jg_C = zeros(6,6*NL);

% Variablen initialisieren (zur Kompilierbarkeit in Simulink)
% Zuordnung mit (:) weiter unten deshalb auch notwendig.
% Dadurch werden Fehler bei der Variablenzuweisung (falsche Dimension)
% leichter erkannt.
ax_0 = NaN(3,1);
r_0_j_C = NaN(3,1);
r_0_0_j = NaN(3,1);
R_0_j = NaN(3,3);
R_0_i = NaN(3,3);
r_0_0_i = NaN(3,1);
r_0_i_C = NaN(3,1);

%% Jacobi berechnen
R_0_i(1:3,1:3) = T_c_mdh(1:3, 1:3,link_index+1);
r_0_0_i(:) = T_c_mdh(1:3, 4, link_index+1);
r_0_i_C(:) = R_0_i * (r_i_i_C);


j = link_index; % Die Indizes j und k haben die Basis als 0.
for tmp = 1:NL
  % Ortsvektor des aktuellen Gelenkes "j"
  r_0_0_j(:) = T_c_mdh(1:3, 4, j+1);
  R_0_j(1:3,1:3) = T_c_mdh(1:3, 1:3, j+1);
  
  % Berechne Vektor vom aktuellen Gelenk zum betrachteten Punkt C
  r_0_j_i = -r_0_0_j + r_0_0_i;
  r_0_j_C(:) = r_0_j_i + r_0_i_C;

  % Geometrische Jacobi-Matrix für alle Achsen des Körper-KS bestimmen.
  % Normalerweise wird eine Schnittmomentkomponente auf ein Drehgelenk
  % projiziert.
  for i_xyz = 1:3
    ax_0(:) = T_c_mdh(1:3,i_xyz,j+1); % Koordinatenachse für Projektion
    
    % Spalte für Kräfte: Entspricht Schubgelenk, [1], Gl. (4.18)
    jt = ax_0;
    jr = zeros(3,1);
    Jg_C(:,(j)*6+i_xyz) = [jt; jr];
    
    % Spalte für Momente: Entspricht Drehgelenk, [1], Gl. (4.19)
    jt = cross(ax_0, r_0_j_C); % Hebelarm vom Gelenk zum Punkt
    jr = ax_0;
    Jg_C(:,(j)*6+i_xyz+3) = [jt; jr];
  end
  % Indizes tauschen: Kinematische Kette weiter Richtung Basis entlanggehen
  if j == 0%-1
    % Nächster wäre Basis. Ist bereits berechnet.
    return;
  end
  j = v(j); % Index mit Basis als 1
end
