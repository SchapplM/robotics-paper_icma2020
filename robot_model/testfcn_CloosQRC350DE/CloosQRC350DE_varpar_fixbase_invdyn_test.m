% Test function for inverse dynamics from symbolic calculations
%
% Quellen:
% [KhalilDombre2002] Modeling, Identification and Control of Robots

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-23 21:15
% Revision: 3f22bf868ffa24e21e77a0fe3b46e78b2d6fdc1f (2020-06-23)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-03
% (C) Institut für Regelungstechnik, Universität Hannover

clc
clear

NQJ = 6;
NJ = 6;
NL = 7;
KINCONSTR = logical(1); %#ok<LOGL>
robot_name = 'CloosQRC350DE';

%% Parameter
TSS = CloosQRC350DE_varpar_testfunctions_parameter();
for f = fields(TSS)'
  eval(sprintf('%s=TSS.%s;',f{1},f{1}));
end
% Es werden fixed-base Funktionen getestet. Der Basis-Körper wird also mit
% NaN belegt und darf keinen Einfluss haben.

%% Prüfe, ob kinematische Zwangsbedingungen vorliegen
thrfac = 1; % Faktor zur Erhöhung der Schwellwerte für Fehlererkennung
NJvirt = length(find(sigma==2));
NJreal = NJ - NJvirt;
if KINCONSTR
  thrfac = 10; % bei ZB sind treten mehr Rechenfehler auf (da aufwändiger)
end

%% Test kinetic Energy
for i = 1:n
  q = Q(i,:)';
  qD = QD(i,:)';
  % calculate kinetic energy with optimized function
  T_func = CloosQRC350DE_energykin_fixb_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);
  % calculate kinetic energy with mass matrix
  Mq = CloosQRC350DE_inertiaJ_slag_vp1(q, ...
    pkin, m, rSges, Icges);
  T_M = 1/2 * qD' * Mq*qD;

  % compare both
  Delta = T_func-T_M;
  if any(abs(Delta) > thrfac*1e6*eps(T_func))
    error('Kinetic Energy from Mass Matrix does not match with directly calculated. %1.5e', ...
      max(abs(Delta)));
  end
end
fprintf('Tested kinetic energy/inertia for %d random joint angles for %s\n', ...
  n, robot_name);


%% Test coriolis vector and Matrix
for i = 1:n
  q = Q(i,:)';
  qD = QD(i,:)';

  % calculate coriolis vector directly from optimized function
  cq = CloosQRC350DE_coriolisvecJ_fixb_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);

  % calculate with coriolis matrix
  Cq = CloosQRC350DE_coriolismatJ_fixb_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);
  cq_mat = (Cq*qD);

  % compare both
  Delta = cq-cq_mat;
  relErr = (cq_mat./ cq)-1; % bei großen Absolutwerten ist der relative Fehler verlässlicher
  if any( abs(Delta) > thrfac*1e7*eps(1+max(abs(cq))) ) && max(abs(relErr)) > 1e-4
    error('coriolis torques from vector and matrix do not match. %1.5e', ...
      max(abs(Delta)));
  end

end
fprintf('Tested coriolis vector/matrix for %d random joint angles for %s\n', ...
  size(Q, 1), robot_name);

%% Teste Massenmatrix-Zeitableitung gegen Coriolis-Matrix
% [KhalilDombre2002] p. 198 9.3.3.3 e)
for i = 1:n
  q = Q(i,:)';
  qD = QD(i,:)';

  % Coriolis Matrix
  Cq = CloosQRC350DE_coriolismatJ_fixb_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);

  % Inertia Matrix time derivative
  MqD = CloosQRC350DE_inertiaDJ_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);

  % Calculate test expression
  Test = MqD-2*Cq;

  % Test skew symmetry of Test-Matrix
  for ii = 1:NQJ
    for jj = 1:NQJ
      % Muss für Diagonal- und Nebendiagonaleinträge Null sein
      test_iijj_abs = Test(ii,jj) + Test(jj,ii);
      if ii == jj
        % Diagonaleintrag muss Null sein
        test_iijj_rel = Test(ii,jj) ./ max(abs(MqD(:)));
      else
        % Nebendiagonaleinträge müssen antisymmetrisch sein
        test_iijj_rel = Test(ii,jj) ./ Test(jj,ii) + 1;
      end
      if abs(test_iijj_abs) > thrfac*1e9*eps(1+max(abs(Cq(:)))) && abs(test_iijj_rel) > 1e-4
        error('Test-Matrix is not skew-symmetric');
      end
    end
  end
end
fprintf('Tested Inertia time derivative/coriolis matrix for %d random joint angles for %s\n', ...
  size(Q, 1), robot_name);

%% Alle Funktionen mit Parametersatz 2 testen
for i = 1:n
  q = Q(i,:)';
  qD = QD(i,:)';
  qDD = QDD(i,:)';
  g_base = G(i,:)';
  % Gravitation
  taug_par1 = CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1(q, g_base, ...
    pkin, m, rSges);
  taug_par2 = CloosQRC350DE_gravloadJ_floatb_twist_slag_vp2(q, g_base, ...
    pkin, m, mrSges);
  Delta_g = taug_par1 - taug_par2;
  Delta_g_rel = Delta_g ./ taug_par2;
  if any( abs(Delta_g) > thrfac*1e6*eps(1+max(abs(taug_par1))) & abs(Delta_g_rel) > 1e-4)
    error('Gravity torques do not match between par1/par2.');
  end
  % Coriolis
  tauc_par1 = CloosQRC350DE_coriolisvecJ_fixb_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);
  tauc_par2 = CloosQRC350DE_coriolisvecJ_fixb_slag_vp2(q, qD, ...
    pkin, m, mrSges, Ifges);
  Delta_c = tauc_par1 - tauc_par2;
  Delta_c_rel = Delta_c ./ tauc_par2;
  if any( abs(Delta_c) > thrfac*1e8*eps(1+max(abs(tauc_par1))) & abs(Delta_c_rel) > 1e-4)
    error('Coriolis vectors do not match between par1/par2.');
  end
  Cq_par1 = CloosQRC350DE_coriolismatJ_fixb_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);
  Cq_par2 = CloosQRC350DE_coriolismatJ_fixb_slag_vp2(q, qD, ...
    pkin, m, mrSges, Ifges);
  Delta_Cq = Cq_par1(:)-Cq_par2(:);
  Delta_Cq_rel = Delta_Cq ./ Cq_par1(:);
  if any( abs(Delta_Cq) > thrfac*2e6*eps(1+max(abs(Cq_par1(:)))) & abs(Delta_Cq_rel) >1e-4 )
    error('Coriolis matrices do not match between par1/par2.');
  end
  % Inertia
  Mq_par1 = CloosQRC350DE_inertiaJ_slag_vp1(q, ...
    pkin, m, rSges, Icges);
  Mq_par2 = CloosQRC350DE_inertiaJ_slag_vp2(q, ...
    pkin, m, mrSges, Ifges);
  Delta_M = Mq_par1(:)-Mq_par2(:);
  Delta_M_rel = Delta_M ./ Mq_par2(:);
  if any( abs(Delta_M) > thrfac*1e6*eps(1+max(abs(Mq_par1(:)))) & abs(Delta_M_rel) > 1e-4)
    error('Inertia matrices do not match between par1/par2.');
  end
  % Inverse Dynamics
  tau_par1 = CloosQRC350DE_invdynJ_fixb_slag_vp1(q, qD, qDD, g_base, ...
    pkin, m, rSges, Icges);
  tau_par2 = CloosQRC350DE_invdynJ_fixb_slag_vp2(q, qD, qDD, g_base, ...
    pkin, m, mrSges, Ifges);
  Delta_tau = tau_par1-tau_par2;
  Delta_tau_rel = Delta_tau./tau_par2; % Lasse auch einen relativen Fehler zu, da bei einigen Systemen dieser Term anfällig für Rechenfehler ist
  if any( abs(Delta_tau) > thrfac*1e7*eps(1+max(abs(tau_par1))) & abs(Delta_tau_rel) > 1e-4)
    error('inverse dynamics joint torque vectors do not match between par1/par2.');
  end
  % Potentielle Energie
  U_par1 = CloosQRC350DE_energypot_fixb_slag_vp1(q, g_base, ...
    pkin, m, rSges);
  U_par2 = CloosQRC350DE_energypot_fixb_slag_vp2(q, g_base, ...
    pkin, m, mrSges);
  if any(abs(U_par1-U_par2) > thrfac*1e6*eps(U_par1))
    error('Potential energies do not match between par1/par2.');
  end
  % Kinetische Energie
  T_par1 = CloosQRC350DE_energykin_fixb_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);
  T_par2 = CloosQRC350DE_energykin_fixb_slag_vp2(q, qD, ...
    pkin, m, mrSges, Ifges);
  if any(abs(T_par1-T_par2) > thrfac*1e6*eps(1+T_par1))
    error('Kinetic energies do not match between par1/par2.');
  end
end
fprintf('Tested Functions with parameter sets 1 and 2 for %d random joint angles for %s\n', ...
  size(Q, 1), robot_name);

%% Komplettfunktion für inverse Dynamik prüfen
for i = 1:n
  q = Q(i,:)';
  qD = QD(i,:)';
  qDD = QDD(i,:)';
  g_base = G(i,:)';
  % Gravitation
  taug_id = CloosQRC350DE_invdynJ_fixb_slag_vp1(q, zeros(NQJ,1), zeros(NQJ,1), g_base, ...
    pkin, m, rSges, Icges);
  taug_sp = CloosQRC350DE_gravloadJ_floatb_twist_slag_vp1(q, g_base, ...
    pkin, m, rSges);
  Delta_g = taug_id - taug_sp;
  Delta_g_rel = Delta_g ./ taug_sp;
  if any(abs(Delta_g) > thrfac*1e6*eps(max(abs(taug_id))) & abs(Delta_g_rel) > 1e-4)
    error('gravity torques from inverse dynamics and dedicated function do not match. %1.5e', ...
      max(abs(Delta_g)));
  end
  % Coriolis
  tauc_id = CloosQRC350DE_invdynJ_fixb_slag_vp1(q, qD, zeros(NQJ,1), zeros(3,1), ...
    pkin, m, rSges, Icges);
  tauc_sp = CloosQRC350DE_coriolisvecJ_fixb_slag_vp1(q, qD, ...
    pkin, m, rSges, Icges);
  Delta_c = tauc_id - tauc_sp;
  Delta_c_rel = Delta_c ./ tauc_sp;
  if any(abs(Delta_c) > thrfac*5e6*eps(1+max(abs([tauc_id;tauc_sp]))) & abs(Delta_c_rel)>1e-4)
    error('Coriolis torques from inverse dynamics and dedicated function do not match. %1.5e', ...
      max(abs(Delta_c)));
  end
  % Inertia
  taua_id = CloosQRC350DE_invdynJ_fixb_slag_vp1(q, zeros(NQJ,1), qDD, zeros(3,1), ...
    pkin, m, rSges, Icges);
  Mq = CloosQRC350DE_inertiaJ_slag_vp1(q, ...
    pkin, m, rSges, Icges);
  taua_sp = Mq*qDD;
  Delta_a = taua_id - taua_sp;
  Delta_a_rel = Delta_c ./ Delta_a;
  if any(abs(Delta_a) > thrfac*1e6*eps(max(abs(taua_sp))) & abs(Delta_a_rel) > 1e-4)
    error('Inertial torques from inverse dynamics and dedicated function do not match. %1.5e', ...
      max(abs(Delta_a)));
  end
end
fprintf('Tested inverse dynamics for %d random joint angles for %s\n', ...
  size(Q, 1), robot_name);

%% Energetische Konsistenz prüfen
n = 5; % Unterschiedliche Anfangskonfigurationen
t_End = 1;
n_complete = 0;
for ii = 1:n
  q0 = Q(ii,:)';
  qD0 = QD(ii,:)';
  g_base = G(ii,:)';

  % Zustand des Systems (Position und Geschwindigkeit)
  x0 = [q0; qD0];
  % Wähle nur Gelenk-Transformationen, die zu realen Körpern gehören (nicht: virtuell)
  if KINCONSTR
    % System mit Zwangsbedingungen: Alle q-Größen sind bereits auf
    % Minimalkoordinaten bezogen. Eine Reduktion ist nicht notwendig.
    % Die sigma=2-Transformationen dienen nur zur Angleichung der Modelle.
    NQJvirt = 0;
    NQJreal = NQJ;
  else
    % System ohne Zwangsbedingungen. Eventuell Offene Baumstruktur eines
    % schließbaren Systems. Transformationen mit sigma=2 werden nicht
    % genommen. Bei normalen seriellen Ketten gibt es kein sigma=2
    NQJvirt = length(find(sigma==2));
    NQJreal = NQJ - NQJvirt;
  end
  selectEntryM = @(x,n) x(n,n);
  selectEntry = @(x,n) x(n);
  % Funktion zur Berechnung der Beschleunigung (Inverse Massenmatrix
  % multipliziert mit Beschleunigungsmoment)
  odeM = @(x) CloosQRC350DE_inertiaJ_slag_vp1(x(1:NQJ), ...
    pkin, m, rSges, Icges);
  odeInvDyn = @(x) CloosQRC350DE_invdynJ_fixb_slag_vp1(x(1:NQJ), ...
    x(NQJ+1:2*NQJ), zeros(NQJ,1), g_base, pkin, m, rSges, Icges);
  odefun = @(t, x) ([x(NQJ+1:2*NQJ); -selectEntryM(odeM(x),1:NQJreal) ...
    \ selectEntry(odeInvDyn(x),1:NQJreal); zeros(NQJvirt,1)]);
  % Testen der Funktionen
  odeM(x0); odeInvDyn(x0); odefun(0, x0); 
  % Numerische Integration konfigurieren und durchführen
  options = odeset('MaxStep',1e-3); % mit automatischer Schrittweite ist der Fehler zu groß
  try
    SolverOutput = ode45(odefun,[0 t_End],x0, options);
  catch err
    if KINCONSTR
      warning('Versuch %d/%d: Keine Lösung der direkten Dynamik möglich. Wahrscheinlich Verlassen des gültigen Gelenkraums', ii, n);
      continue
    else
      % Es gibt keine kinematischen Zwangsbedingungen, die direkte Dynamik
      % muss also funktionieren
      raise(err);
    end
  end
  t = SolverOutput.x;
  % Bestimme Zeitpunkt der ersten Grenzüberschreitung
  if KINCONSTR
    I_viol = false(length(t), 1);
    for j = 1:length(q0)
      I_viol = I_viol | SolverOutput.y(j,:)' < q_min(j) | SolverOutput.y(j,:)' > q_max(j);
    end
    I_end = find(I_viol, 1)-1;
  end
  if ~KINCONSTR || isempty(I_end)
    I_end = length(t);
  end
  % Berechne Energie aus den Ausgaben
  E = NaN(length(t), 3);
  tau_Acc = NaN(length(t), NQJ);
  for i = 1:length(t)
    q = SolverOutput.y(1:NQJ,i);
    qD = SolverOutput.y(NQJ+1:2*NQJ,i);
    E(i,1) = CloosQRC350DE_energykin_fixb_slag_vp1(q, qD, pkin, m, rSges, Icges);
    E(i,2) = CloosQRC350DE_energypot_fixb_slag_vp1(q, g_base, pkin, m, rSges);
    tau_Acc(i,:) = CloosQRC350DE_invdynJ_fixb_slag_vp1(q, qD, zeros(NQJ,1), g_base, pkin, m, rSges, Icges)';
  end
  E(:,3) = sum(E(:,1:2),2); % Gesamtenergie
  % Ergebnisse zeichnen
  change_current_figure(1);clf;
  sgtitle(sprintf('Energie-Test: %s', robot_name), 'interpreter', 'none'); grid on;
  subplot(2,2,1); hold on; grid on;
  hdl1=plot(t(1:I_end), E(1:I_end,:));ylabel('Energie [J]');
  set(gca, 'ColorOrderIndex', 1); hdl2=plot(t(I_end+1:end), E(I_end+1:end,:), '--');
  lt1={'T', 'U', 'Ges'}; lt2={'T (after limit viol.)', 'U (after limit viol.)', 'Ges (after limit viol.)'};
  if isempty(hdl2), legend(hdl1, lt1);
  else,             legend([hdl1;hdl2], [lt1(:)', lt2(:)']); end
  subplot(2,2,2); hold on; grid on;
  plot(t(1:I_end), SolverOutput.y(1:NQJ,1:I_end)'); ylabel('Gelenkkoord [rad], [m]');
  set(gca, 'ColorOrderIndex', 1); plot(t(I_end+1:end), SolverOutput.y(1:NQJ,I_end+1:end)', '--');
  subplot(2,2,3); hold on; grid on;
  plot(t(1:I_end), SolverOutput.y(NQJ+1:2*NQJ,1:I_end)'); ylabel('Geschwindigkeit [rad/s], [m/s]');
  set(gca, 'ColorOrderIndex', 1); plot(t(I_end+1:end), SolverOutput.y(NQJ+1:2*NQJ,I_end+1:end)', '--'); 
  subplot(2,2,4); hold on; grid on;
  plot(t(1:I_end), tau_Acc(1:I_end,:)); ylabel('Beschleunigungsmoment [Nm], [m]');
  set(gca, 'ColorOrderIndex', 1); plot(t(I_end+1:end), tau_Acc(I_end+1:end,:), '--');
  linkxaxes
  % Berechne den Energiefehler nur bezogen auf den Gelenkbereich, der
  % innerhalb der Definitionsgrenzen liegt (bei Systemen mit ZB)
  E_Delta = E(1,3) - E(I_end,3);
  E_Delta_Rel = E_Delta/max(abs(E(:))); % beziehen relativen Fehler auf maximale Teilkomponente
  if (~KINCONSTR && (abs(E_Delta) > 5e-4 && abs(E_Delta_Rel) > 5e-2) ) ... % Nutze absoluten und relativen Fehler für Systeme ohne Zwangsbedingungen
  || (KINCONSTR  && abs(E_Delta_Rel) > 5e-2) % mit Zwangsbed. sind die Systeme komplizierter. Nur relativer Fehler da mehr Rundungsfehler bei Rechnung möglich
    error('Das System ist nicht energetisch konsistent');
  end
  n_complete = n_complete + 1;
end
fprintf('Tested energetic consistency for free movement with %d starting configurations for %s. %d configurations not successful.\n', ...
  n_complete, robot_name, n-n_complete);


