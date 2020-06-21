% Validate the generated model from SA Jonas Diekmeyer (main.mw) with a
% second implementation from [1]
% Use different test scenarios to isolate potential errors in the model
% 
% Preliminary result:
% * test case 3 not working yet
% 
% Dependencies:
% * https://github.com/SchapplM/matlab_toolbox
% * https://github.com/SchapplM/robotics-dep-ext
% * https://github.com/SchapplM/robotics-toolbox
% 
% References:
% [1] https://github.com/SchapplM/robsynth-modelgen
% [2] Diekmeyer, Jonas: Identifikation der inversen Dynamik eines seriellen
% Roboters im geschlossenen Regelkreis, Studienarbeit, imes, LUH, 2018

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-06
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Initialize Path
this_path = fileparts( mfilename('fullpath') );
% path of the implementation from [2]
addpath(fullfile(this_path, '..'));
addpath(fullfile(this_path, '..', 'generation', 'kinematics'));
addpath(fullfile(this_path, '..', 'generation', 'dynamics'));
% Location of the second implementation [1]
addpath(fullfile(this_path, '..', 'matlabfcn_CloosQRC350DE'));
addpath(fullfile(this_path, '..', 'matlabfcn_CloosQRC350OL'));

% test scenarios:
% 1: statics only
% 2: dynamics only, only inertia parameters. CoM zero
% 3: dynamics only, no inertia (around CoM), but CoM unequal to zero
for ts = 1:3
  %% Define Kinematics Parameters and Test Settings
  n = 100;
  N = 6;
  q_min = -pi*ones(N,1);q_max = pi*ones(N,1);
  Q = repmat(q_min',n,1) + rand(n,N).*repmat(q_max'-q_min',n,1);
  QD = (0.5-rand(n, N))*pi;
  QDD = (0.5-rand(n, N))*pi;
  % Select single poses for statical evaluation
  Q(1,:) = [0,0,0,0,0,0]*pi/180;
  Q(2,:) = [0,0,-90,0,0,0]*pi/180;
  Q(3,:) = [0,90,-90,0,0,0]*pi/180;
  Q(4,:) = [0,90, 0,0,0,0]*pi/180;
  if ts == 1 % only statics
    QD(:,:)=0;
    QDD(:,:)=0;
  elseif ts == 2
    % alle velocities and accelerations set unequal to zero
  elseif ts == 3
    % Erste Konfigurationen nur mit Beschleunigung testen
    QD(1:N+1,:) = 0;
    QDD(1,:) = 0;
    QDD(2:N+1,1:N) = eye(N);
  end

  % kinematic parameters of the robot
  L1 = 640*1e-3;
  L2 = 250*1e-3;
  L3 = 630*1e-3;
  L4 = 196*1e-3;
  L5 = 805*1e-3;
  L6 = 100*1e-3;
  kDG= 227/1200; % ratio of the differential gear (empiric value, from [2])
  pkin = [L1,L2,L3,L4,L5,L6,kDG]';

  % Constant transformations to account for different frame definitions
  % (model from [2] uses standard DH, reference from [1] uses modified)
  T_mdl1_mdl2 = zeros(4,4,7);
  T_mdl1_mdl2(:,:,1) = eye(4);
  T_mdl1_mdl2(:,:,2) = trotx(pi/2)*transl([-L2;0;0]);
  T_mdl1_mdl2(:,:,3) = transl([-L3;0;0]);
  T_mdl1_mdl2(:,:,4) = trotx(pi/2)*transl([-L4;0;0]);
  T_mdl1_mdl2(:,:,5) = trotx(-pi/2);
  T_mdl1_mdl2(:,:,6) = trotx(pi/2);
  T_mdl1_mdl2(:,:,7) = eye(4);

  %% Define Dynamics Parameters and Test Settings
  g = 9.81; % only z component
  if ts == 2 || ts == 3
    g = 0;
  end
  % Define dynamics parameters
  rc_all_sdh = 0.3*(0.5-rand(N,3)); % All center of mass coordinates in body frames

  m = rand(N,1); % masses of all links (are positive due to rand() function)
  if ts == 3
    m = 1*ones(N,1);
    m([1 3:end]) = 0;
  end

  % Set parameters to zero with assumed symmetries (see Maple worksheet)
  % Components refer to SDH frames
  % rC2y=0,rC4x=0,rC4z=0,rC5x=0,rC5y=0,rC6x=0,rC6y=0
  rc_all_sdh(1,1) = 0;
  rc_all_sdh(1,3) = 0;
  rc_all_sdh(2,2) = 0;
  rc_all_sdh(4,1) = 0;
  rc_all_sdh(4,3) = 0;
  rc_all_sdh(5,1) = 0;
  rc_all_sdh(5,2) = 0;
  rc_all_sdh(6,1) = 0;
  rc_all_sdh(6,2) = 0;
  if ts == 2
    rc_all_sdh(1:6,:) = 0;
  end
  if ts == 3
    rc_all_sdh(1:6,:) = 0;
    rc_all_sdh(2,3) = 0.1;
  end
  Ic_pa = rand(N,3); % inertia of all links around their center of mass in principal axes
  Ic_all_sdh = NaN(N,6); % inertia of all links around their center of mass in body frame
  for i = 1:N
    R_pa = eulxyz2r(rand(3,1)); % random principal axes
    % inertia tensor in body frame: make sure the eigenvalues are positive and the tensor is positive definite
    Ic_all_sdh(i,:) = inertiamatrix2vector(R_pa*diag(Ic_pa(i,:))*R_pa');
  end
  % Set assumptions on products of inertia (again in sdh frames)
  % J2xy=0,J2yz=0,J4xy=0,J4xz=0,J4yz=0,J5xy=0,J5xz=0,J5yz=0,J6xy=0,J6xz=0,J6yz=0
  Ic_all_sdh(2,[4 6]) = 0;
  Ic_all_sdh(4,4:6) = 0;
  Ic_all_sdh(5,4:6) = 0;
  Ic_all_sdh(6,4:6) = 0;
  if ts == 3
    Ic_all_sdh(:,:) = 0;
  end
  % Get inertial parameters from barycentric parameters
  [mrc_all_sdh, ... % first moment of all links (mass times center of mass)
   If_all_sdh] = ... % second moment of all links (inertia around body frame origins)
    inertial_parameters_convert_par1_par2(rc_all_sdh, Ic_all_sdh, m);

  % Convert dynamics parameters from mdh to sdh frames
  rc_all_mdh = NaN*rc_all_sdh;
  Ic_all_mdh = NaN*Ic_all_sdh;
  for i = 1:N
    T_mdh_sdh_i = invtr(T_mdl1_mdl2(:,:,i+1));
    % Schwerpunkt mit Koordinatentransformation
    rc_all_mdh(i,:) = eye(3,4) * T_mdh_sdh_i*[rc_all_sdh(i,:)';1];
    % Trägheitstensor mit Basiswechsel (kein Steiner-Anteil notwendig, da auf
    % Schwerpunkt bezogen)
    R_m_s = T_mdh_sdh_i(1:3,1:3);
    Ic_sdh = inertiavector2matrix(Ic_all_sdh(i,:));
    Ic_mdh = R_m_s * Ic_sdh * R_m_s';
    Ic_all_mdh(i,:) = inertiamatrix2vector(Ic_mdh);
  end
  [mrc_all_mdh, If_all_mdh] = inertial_parameters_convert_par1_par2(rc_all_mdh, Ic_all_mdh, m);

  % Inertia of motor and gear
  JA_all = zeros(N,1);

  % friction parameters
  fc = zeros(N,1);
  fv = zeros(N,1);

  % Initialize robot class
  PS = struct('beta',  NaN(N,1), 'b', NaN(N,1), ...
              'alpha', NaN(N,1), 'a', NaN(N,1), ...
              'theta', NaN(N,1), 'd', NaN(N,1), ...
              'sigma', zeros(N,1), 'offset', NaN(N,1), ...
              'pkin', [], 'v', uint8(0:N-1)', ...
              'mu', ones(N,1), ...
              'NJ', N, 'NL', N+1, 'NQJ', N, ...
              'qmin', NaN(N,1), 'qmax', NaN(N,1), 'vmax', NaN(N,1), 'qref', zeros(N,1));

  [PS.beta, PS.b, PS.alpha, PS.a, PS.theta, PS.d, PS.offset] = CloosQRC350DE_pkin2mdhparam(pkin);
  PS.pkin = pkin;
  RS = SerRob(PS, 'CloosQRC350DE');
  RS.fill_fcn_handles(false);
  Ic_all_mdh_plot = Ic_all_mdh;
  Ic_all_mdh_plot(:,1:3) = 1; % for plotting the CoM as spheres
  RS.update_dynpar1([NaN;m], [NaN(1,3);rc_all_mdh], [NaN(1,6);Ic_all_mdh_plot]);

  %% Test Kinematics
  for i = 1:n
    q = Q(i,:)';
    % forward kinematics from reference implementation [1]
    Tc_mdl2 = CloosQRC350DE_fkine_fixb_rotmat_mdh_sym_varpar(q, pkin);
    jv = CloosQRC350DE_kinconstr_expl_mdh_num_varpar(q, pkin);
    % forward kinematics from implementation [2] for this paper
    Tc_mdl1 = zeros(4,4,7);
    Tc_mdl1(1:4,1:4,1) = eye(4);
    % Rotation matrices
    Tc_mdl1(1:3,1:3,2) = A01(q(1));
    Tc_mdl1(1:3,1:3,3) = A02(q(1), q(2));
    Tc_mdl1(1:3,1:3,4) = A03(q(1), q(2), q(3));
    Tc_mdl1(1:3,1:3,5) = A04(q(1), q(2), q(3), q(4));
    Tc_mdl1(1:3,1:3,6) = A05(q(1), q(2), q(3), q(4), q(5));
    Tc_mdl1(1:3,1:3,7) = A06(q(1), q(2), q(3), q(4), q(5), q(6));
    % body frame origins
    Tc_mdl1(1:3,4,2) = r01(L1, L2, q(1));
    Tc_mdl1(1:3,4,3) = r02(L1, L2, L3, q(1), q(2));
    Tc_mdl1(1:3,4,4) = r03(L1, L2, L3, L4, q(1), q(2), q(3));
    Tc_mdl1(1:3,4,5) = r04(L1, L2, L3, L4, L5, q(1), q(2), q(3));
    Tc_mdl1(1:3,4,6) = r05(L1, L2, L3, L4, L5, q(1), q(2), q(3));
    Tc_mdl1(1:3,4,7) = r06(L1, L2, L3, L4, L5, L6, q(1), q(2), q(3), q(4), q(5));
    Tc_mdl1(4,4,:) = 1;
    % compare
    for j = 1:7
      Tdiff = (Tc_mdl1(:,:,j)*T_mdl1_mdl2(:,:,j)) \ Tc_mdl2(:,:,j);
      T_test = Tdiff-eye(4);
      phi_diff = r2eulxyz(Tdiff(1:3,1:3));
      if any(abs(T_test(:))>1e-10)
        error('Kinematics do not match between the two implementations at frame %d', j-1);
      end
    end
  end
  fprintf('Scenario %d: Tested the inverse kinematics for %d random configurations\n', ts, n);

  %% Test Dynamics with different implementations
  mpv = Minimalparameter_fcn(pkin, m, rc_all_sdh, Ic_all_sdh, JA_all, fc, fv)';
  for i = 1:n
    q = Q(i,:)';
    qD = QD(i,:)';
    qDD = QDD(i,:)';

    % Test against second implementation
    tau_mdl1 = InvDyn_fcn(q, qD, qDD, g, pkin, m, rc_all_sdh, Ic_all_sdh, JA_all)';
    tau_mdl2 = CloosQRC350DE_invdynJ_fixb_slag_vp2(q, qD, qDD, [0;0;-g], pkin, [NaN;m], ...
      [NaN(1,3);mrc_all_mdh], [NaN(1,6);If_all_mdh]);
    tau_test_mdl2 = tau_mdl2 - tau_mdl1;
    if any(abs(tau_test_mdl2) > 1e-10)
      disp(tau_mdl2');
      disp(tau_mdl1');
      error('Joint torques do not match with second implementation for config. %d', i);
    end

    % Test against regressor form
    Phi = RegressorMatrix(qDD', qD', g, 1000, pkin, q');
    tau_from_reg = Phi * mpv;
    tau_test_reg = tau_from_reg - tau_mdl1;
    if any(abs(tau_test_reg) > 1e-10)
      disp('from regressor:');
      disp(tau_from_reg');
      disp('from function:');
      disp(tau_mdl1');
      % Plot Robot
      Tc_mdl2 = CloosQRC350DE_fkine_fixb_rotmat_mdh_sym_varpar(q, pkin);
      % Plot robot in test configuration
      s_plot = struct( 'ks', [1:RS.NJ, RS.NJ+2], 'straight', 0, 'mode', 3);
      figure(1);clf;
      hold on; grid on;
      xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
      view(3);
      RS.plot( q, s_plot );
      title(sprintf('Robot in test configuration %d', i));
      % Plot SDH frames (relevant for some parameters)
      for j = 1:N
        T_sdh_j = Tc_mdl2(:,:,j+1)*invtr(T_mdl1_mdl2(:,:,j+1));
        trplot(T_sdh_j, 'frame', sprintf('sdh,%d',j), 'rgb', 'length', 0.30)
      end
      error('Joint torques do not match with regressor form for config. %d. abs error %1.3e', ...
        i, max(abs(abs(tau_test_reg))));
    end
  end
  fprintf('Scenario %d: Tested the inverse dynamics for %d random configurations\n', ts, n);

  %% Test Case for statics only
  if ts == 1
    continue
  end
  %% Get Inertia Matrix from Regressor Form
  for i = 1:n
    q = Q(i,:)';
    MM_fromreg = zeros(N,N);
    MM_fromfcn = zeros(N,N);
    % calculate inertia matrix column-wise
    for j = 1:N
      qDD_test = zeros(N,1);
      qDD_test(j) = 1; % unit acceleration in joint j to build column
      Phi_j = RegressorMatrix(qDD_test', zeros(1,N), 0, 1000, pkin, q');
      tau_j = Phi_j * mpv;
      MM_fromreg(:,j) = tau_j;
      tau_j = InvDyn_fcn(q, zeros(N,1), qDD_test, 0, pkin, m, rc_all_sdh, Ic_all_sdh, JA_all)';
      MM_fromfcn(:,j) = tau_j;
    end
    % inertia matrix from [1]
    MM_mdl2 = CloosQRC350DE_inertiaJ_slag_vp2(q, pkin, [NaN;m], ...
      [NaN(1,3);mrc_all_mdh], [NaN(1,6);If_all_mdh]);
    % compare
    test_MM_mdl2 = MM_mdl2 - MM_fromfcn;
    if any(abs(test_MM_mdl2(:)) > 1e-10)
      error('Inertia matrix does not match with reference implementation for config. %d', i);
    end
    test_MM_reg = MM_fromreg - MM_fromfcn;
    if any(abs(test_MM_reg(:)) > 1e-10)
      error('Inertia matrix does not match with regressor form for config. %d', i);
    end
  end
  fprintf('Scenario %d: Tested the inertia matrix for %d random configurations\n', ts, n);
end
