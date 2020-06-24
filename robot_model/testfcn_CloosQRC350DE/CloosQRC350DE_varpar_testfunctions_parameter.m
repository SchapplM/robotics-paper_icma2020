% Einstellungen für Modultests von CloosQRC350DE generieren:
% Enthält Definitionen für Parameter und zufällige Konfigurationen
%
% Ausgabe:
% TSS
%   Struktur mit Einstellungen für die Modultests

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-03
% (C) Institut für Regelungstechnik, Universität Hannover

function TSS = CloosQRC350DE_varpar_testfunctions_parameter()

%% General Definitions
NQJ = 6; % Number of generalized coordinates
NJ = 6; % number of joints (number of rows in Denavit-Hartenberg table)
NL = 7; % number of links (rigid bodies, including the base)
robot_name = 'CloosQRC350DE'; % prefix in all function names and simulink models and blocks

%% Kinematik-Parameter definieren
% These parameters may be overwritten down in this script
a = rand(NJ,1); % Kinematic length of MDH notation
alpha = zeros(NJ,1); % Kinematic angle of MDH notation
d = rand(NJ,1); % Kinematic length of MDH notation
q_offset = zeros(NJ,1); % Kinematic angle of MDH notation
b = zeros(NJ,1); % Kinematic length of MDH notation
beta = zeros(NJ,1); % Kinematic angle of MDH notation

%% Zufällige Roboterkonfigurationen
% Werden in den Testfunktionen benutzt

n = 100;
% Gelenkwinkel und -zeitableitungen
q_min = -pi*ones(NQJ,1);q_max = pi*ones(NQJ,1);
Q = repmat(q_min',n,1) + rand(n,NQJ).*repmat(q_max'-q_min',n,1);
QD = (0.5-rand(n, NQJ))*pi;
QD(1:NQJ,:)=eye(NQJ);
QDD = (0.5-rand(n, NQJ))*pi;

% Gravitation
G = (0.5-rand(n, 3))*10;
G(:,~[0,0,1]) = 0;

% Basisposition
RB = (0.5-rand(n, 3))*10;

% Basisorientierung
% Zufällige Zeitableitung der Orientierungsdarstellung und ihrer
% Ableitungen
OB = (0.5-rand(n, 3))*pi;
OB(1,:) = 0;
OBD = (0.5-rand(n, 3))*pi;
OBDD = (0.5-rand(n, 3))*pi;

% Winkelgeschwindigkeit und Beschleunigung aus der Orientierungsdarstellung
% berechnen
VB = NaN(n,6);
AB = NaN(n,6);
VB(:,1:3) = (0.5-rand(n, 3))*10;
AB(:,1:3) = (0.5-rand(n, 3))*10;
for i = 1:n
  % Darstellung umwandeln.
  % TODO: Mit Transformationsmatrizen
  VB(i, 4:6) = eulxyzD2omega(OB(i,:)', OBD(i,:)');
  AB(i, 4:6) = eulxyzDD2omegaD(OB(i,:)', OBD(i,:)', OBDD(i,:)');
  
  % Probe:
  T_basevel = eulxyzjac(OB(i,:)');
  obd_test = T_basevel \ VB(i, 4:6)';
  
  if any(abs(obd_test-OBD(i,:)') > 1e-10)
    error('Orientierung und Basisgeschwindigkeit stimmen nicht überein');
  end
end




%% MDH-Parametereinträge auf Zufallswerte setzen
% Aus robot_matlabtmp_par.m
a1 = a(1);
a2 = a(2);
a3 = a(3);
a4 = a(4);
a5 = a(5);
a6 = a(6);
alpha1 = alpha(1);
alpha2 = alpha(2);
alpha3 = alpha(3);
alpha4 = alpha(4);
alpha5 = alpha(5);
alpha6 = alpha(6);
d1 = d(1);
d2 = d(2);
d3 = d(3);
d4 = d(4);
d5 = d(5);
d6 = d(6);
qoffset1 = q_offset(1);
qoffset2 = q_offset(2);
qoffset3 = q_offset(3);
qoffset4 = q_offset(4);
qoffset5 = q_offset(5);
qoffset6 = q_offset(6);
b1 = b(1);
b2 = b(2);
b3 = b(3);
b4 = b(4);
b5 = b(5);
b6 = b(6);
beta1 = beta(1);
beta2 = beta(2);
beta3 = beta(3);
beta4 = beta(4);
beta5 = beta(5);
beta6 = beta(6);

%% Werte für Kinematikparameter direkt eintragen

% Aus CloosQRC350DE/kinematic_parameter_values.m
L1 = 640*1e-3;
L2 = 250*1e-3;
L3 = 630*1e-3;
L4 = 196*1e-3;
L5 = 805*1e-3;
L6 = 100*1e-3;
kDG= 227/1200; % ration of the differential gear (empiric value)

% Aus CloosQRC350DE/parameter_kin_matlab.m
t1 = [L1; L2; L3; L4; L5; L6; kDG;];
pkin = t1;
if isempty(pkin)
  pkin = 0;%Platzhalter-Eingabe
end


%% MDH-Parametereinträge mit vorgegebenen Werten überschreiben
% Aus CloosQRC350DE/parameters_d_matlab.m
t1 = [L1; 0; 0; L5; 0; L6;];
d = t1;

% Aus CloosQRC350DE/parameters_a_matlab.m
t1 = [0; L2; L3; L4; 0; 0;];
a = t1;

% Aus CloosQRC350DE/parameters_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
theta = t1;

% Aus CloosQRC350DE/parameters_b_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
b = t1;

% Aus CloosQRC350DE/parameters_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
beta = t1;

% Aus CloosQRC350DE/parameters_alpha_matlab.m
t1 = -pi / 0.2e1;
t2 = [0; t1; 0; t1; pi / 0.2e1; t1;];
alpha = t2;

% Aus CloosQRC350DE/parameters_qoffset_matlab.m
t1 = [0; -pi / 0.2e1; 0; 0; 0; pi;];
q_offset = t1;

% Aus CloosQRC350DE/parameters_v_matlab.m
t1 = [0; 1; 2; 3; 4; 5;];
v = uint8(t1);

% Aus CloosQRC350DE/parameters_sigma_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
sigma = t1;

% Aus CloosQRC350DE/parameters_mu_matlab.m
t1 = [1; 1; 1; 1; 1; 1;];
mu = t1;

% Aus CloosQRC350DE/kinconstr_index_dependant_joints_matlab.m
t1 = [1; 0; 0; 0; 0; 1;];
Ind_depjoints = t1;
%% Dynamik-Parameter definieren
rSges = rand(NL,3); % All center of mass coordinates in body frames
m = rand(NL,1); % masses of all links (are positive due to rand() function)
Ic_pa = rand(NL,3); % inertia of all links around their center of mass in principal axes
Icges = NaN(NL,6); % inertial of all links around their center of mass in body frame
for i = 1:NL
  R_pa = eulxyz2r(rand(3,1)); % random principal axes
  % inertia tensor in body frame: make sure the eigenvalues are positive and the tensor is positive definite
  Icges(i,:) = inertiamatrix2vector(R_pa*diag(Ic_pa(i,:))*R_pa');
end

% Parameter reduzieren, falls durch Benutzereingabe gefordert.
% Notwendig, damit Dynamikmodell konsistent ist mit den Eingabeparametern
% Das betrifft nur die baryzentrischen Parameter (par1). Die Inertialparameter (par2) werden daraus berechnet.
[m,rSges,Icges]=CloosQRC350DE_dynamics_parameters_modification(pkin,m,rSges,Icges);
% Prüfe, ob die Reduktion konsistent ist. Es dürfen keine zyklischen
% Abhängigkeiten zwischen Masse, Schwerpunkt und Trägheit gesetzt sein.
[m2,rSges2,Icges2]=CloosQRC350DE_dynamics_parameters_modification(pkin,m,rSges,Icges);
if any(abs([m-m2;rSges(:)-rSges2(:);Icges(:)-Icges2(:)])>1e-10)
  error('Bei zweifacher Durchführung der Parameterreduktion ändern sich die Parameter. Nicht konsistent!');
end
% Inertialparameter aus baryzentrischen Parametern berechnen
[mrSges, ... % first moment of all links (mass times center of mass)
 Ifges] = ... % second moment of all links (inertia around body frame origins)
  inertial_parameters_convert_par1_par2(rSges, Icges, m);

%% Set Outputs
TSS = struct('type', 'Test Settings Structure');
% Allgemeine Definitionen
TSS.NQJ = NQJ;
TSS.NJ = NJ;
TSS.NL = NL;
TSS.Ind_depjoints = logical(Ind_depjoints); % Binärindizes der abhängigen Gelenke
% Kinematische Zwangsbedingungen
TSS.NQJ = NQJ;
% Kinematikparameter
TSS.a = a;
TSS.alpha = alpha;
TSS.d = d;
TSS.q_offset = q_offset;
TSS.b = b;
TSS.beta = beta;
TSS.v = v;
TSS.pkin = pkin;
TSS.theta = theta;
TSS.sigma = sigma;
TSS.mu = mu;
TSS.q_min = q_min;
TSS.q_max = q_max;
% Dynamikparameter
TSS.m = m;
TSS.rSges = rSges;
TSS.Icges = Icges;
TSS.mrSges = mrSges;
TSS.Ifges = Ifges;
% Zufällige Konfigurationen für Modultests
TSS.n = n;
TSS.Q = Q;
TSS.QD = QD;
TSS.QDD = QDD;
TSS.RB = RB;
TSS.OB = OB;
TSS.OBD = OBD;
TSS.OBDD = OBDD;
TSS.VB = VB;
TSS.AB = AB;
TSS.G = G;
