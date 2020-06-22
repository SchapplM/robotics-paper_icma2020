% % Author M.Sc. Björn Volkmann

% % The model reduction needs to be run for process 1 and 2 prior to this
% script to compare the reduced models to the models derived fom optimal
% excitation

addpath(genpath('../'));
clear
close all
clc
tic

%% Variablen und Parameter
% Robot Parameters
g = 9.81; % Gravitation
l = [640;250;630;196;805;100]/1000; % Kinematic
ks = 1000; % Steifigkeit der Coulomb-Reibung (M_Reib = fc * tanh(ks*dq))

%% Load Data
fprintf('Load Data \n');

%load Trajectory A
load('Trajectory_A.mat');

qA = q_meas.average_rect;
dqA = dq_meas.average_rect;
ddqA = ddq_meas.average_rect;
tauA = reshape(tau_meas.average, [], 1);
tauAA = abs(tau_meas.average);
VarianceA = tau_meas.variance;
timeA = time_vec_save;

C_a = RegressorMatrix(ddqA,dqA,g,ks,l,qA);%Designmatrix


%load Trajectory B
load('Trajectory_B.mat');

qB = q_meas.average_rect;
dqB = dq_meas.average_rect;
ddqB = ddq_meas.average_rect;
tauB = reshape(tau_meas.average, [], 1);
tauBB = abs(tau_meas.average);
VarianceB = tau_meas.variance;
timeB = time_vec_save;

C_b = RegressorMatrix(ddqB,dqB,g,ks,l,qB);%Designmatrix


%load Trajectory C
load('Trajectory_C.mat');

qC = q_meas.average_rect;
dqC = dq_meas.average_rect;
ddqC = ddq_meas.average_rect;
tauC = reshape(tau_meas.average, [], 1);
tauCC = abs(tau_meas.average);
VarianceC = tau_meas.variance;
timeC = time_vec_save;

C_c = RegressorMatrix(ddqC,dqC,g,ks,l,qC);%Designmatrix


% load Process 1
load('Validation_Data_Process1.mat');

qP1 = deg2rad(deg_q);
dqP1 = deg2rad(deg_dq);
ddqP1 = deg2rad(deg_ddq);
tauP1 = reshape(tau, [], 1);
tauP1P1 = abs(tau);
timeP1 = time_vec_save;

C_P1 = RegressorMatrix(ddqP1,dqP1,g,ks,l,qP1);%Designmatrix


% load Process 2
load('Validation_Data_Process2.mat');

qP2 = deg2rad(deg_q);
dqP2 = deg2rad(deg_dq);
ddqP2 = deg2rad(deg_ddq);
tauP2 = reshape(tau, [], 1);
tauP2P2 = abs(tau);
timeP2 = time_vec_save;

C_P2 = RegressorMatrix(ddqP2,dqP2,g,ks,l,qP2);%Designmatrix


fprintf('Data Loaded \n');

%% Parameter estimation

theta_A = WLS_method(C_a, tauA, VarianceA);

theta_B = WLS_method(C_b, tauB, VarianceB);


%% Calculate Model Error

[rA, ~] = size(tauA);%number of measurement points in trajectoryA
[rB, ~] = size(tauB);%number of measurement points in trajectoryB
[rC, ~] = size(tauC);%number of measurement points in trajectoryC

%Model error for both models in their respective trajectory
%
%   naming convention: e[TRAJECTORY]_[MODEL]
%
eA_A.axes = transpose( mean( reshape( abs(tauA - C_a*theta_A), [], 6)));%absolute
eB_B.axes = transpose( mean( reshape( abs(tauB - C_b*theta_B), [], 6)));%absolute

eA_A.axes_rel = (rA/6) * eA_A.axes(:, 1) ./ transpose( sum(tauAA) );%relative
eB_B.axes_rel = (rB/6) * eB_B.axes(:, 1) ./ transpose( sum(tauBB) );%relative

%Model error for both models in trajectory C
eC_A.axes = transpose( mean( reshape( abs(tauC - C_c*theta_A), [], 6)));%absolute
eC_B.axes = transpose( mean( reshape( abs(tauC - C_c*theta_B), [], 6)));%absolute

eC_A.axes_rel = (rC/6) * eC_A.axes(:, 1) ./ transpose( sum(tauCC) );%relative
eC_B.axes_rel = (rC/6) * eC_B.axes(:, 1) ./ transpose( sum(tauCC) );%relative

% Variance
%model A
sigma_A = zeros(size(theta_A));
W = diag(reshape(1./VarianceA, [], 1));
cov_A = inv( transpose(C_a) * W * C_a );
for k = 1:length(theta_A)
    
    sigma_A(k,1) = 100 * sqrt(abs(cov_A(k,k)))/abs(theta_A(k));
    
end

%model B
sigma_B = zeros(size(theta_B));
W = diag(reshape(1./VarianceB, [], 1));
cov_B = inv( transpose(C_b) * W * C_b );
for k = 1:length(theta_A)
    
    sigma_B(k,1) = 100 * sqrt(abs(cov_B(k,k)))/abs(theta_B(k));
    
end


%% Applying the models A and B to the processes 1 and 2

[rP1, ~] = size(tauP1);%number of measurement points in process 1
[rP2, ~] = size(tauP2);%number of measurement points in process 2

% % % Model A

% % Process 1
eP1_A.axes =  transpose( mean( reshape( abs(tauP1 - C_P1*theta_A), [], 6)));%absolute
eP1_A.axes_rel = (rP1/6) * eP1_A.axes(:, 1) ./ transpose( sum(tauP1P1) );%relative

% % Process 2
eP2_A.axes =  transpose( mean( reshape( abs(tauP2 - C_P2*theta_A), [], 6)));%absolute
eP2_A.axes_rel = (rP2/6) * eP2_A.axes(:, 1) ./ transpose( sum(tauP2P2) );%relative


% % % Model B

% % Process 1
eP1_B.axes =  transpose( mean( reshape( abs(tauP1 - C_P1*theta_B), [], 6)));%absolute
eP1_B.axes_rel = (rP1/6) * eP1_B.axes(:, 1) ./ transpose( sum(tauP1P1) );%relative

% % Process 2
eP2_B.axes =  transpose( mean( reshape( abs(tauP2 - C_P2*theta_B), [], 6)));%absolute
eP2_B.axes_rel = (rP2/6) * eP2_B.axes(:, 1) ./ transpose( sum(tauP2P2) );%relative

%save results
save('optimal_exitation', 'eA_A', 'eC_A' ,'eB_B', 'eC_B' ,'eP1_A', 'eP2_A' ,'eP1_B', 'eP2_B', 'theta_A', 'theta_B', 'sigma_A', 'sigma_B');

%% Ploting

% Model A
tauC = transpose(reshape(tauC, [], 6));
tauC_vgl = transpose( reshape( C_c*theta_A, [], 6) );

%Plotten der Drehmomente
figure(1)
sgtitle('Model A for Trajectory C')
for plott = 1:6 
    subplot(3,2,plott)
    plot(timeA, [tauC(plott,:); tauC_vgl(plott,:)])
    title(strcat('Joint ',num2str(plott)))
    xlabel('time (s)')
    ylabel('torque (Nm)')
end
legend('training', 'validation', 'Location','northeast')

% subplot(3,2,2)
% plot(timeA, [tauC(2,:); tauC_vgl(2,:)])
% title('Axis 2')
% xlabel('Time (s)')
% ylabel('Torque (Nm)')
% legend('Measurement', 'Model')
% 
% subplot(3,2,3)
% plot(timeA, [tauC(3,:); tauC_vgl(3,:)])
% title('Axis 3')
% xlabel('Time (s)')
% ylabel('Torque (Nm)')
% legend('Measurement', 'Model')
% 
% subplot(3,2,4)
% plot(timeA, [tauC(4,:); tauC_vgl(4,:)])
% title('Axis 4')
% xlabel('Time (s)')
% ylabel('Torque (Nm)')
% legend('Measurement', 'Model')
% 
% subplot(3,2,5)
% plot(timeA, [tauC(5,:); tauC_vgl(5,:)])
% title('Axis 5')
% xlabel('Time (s)')
% ylabel('Torque (Nm)')
% legend('Measurement', 'Model')
% 
% subplot(3,2,6)
% plot(timeA, [tauC(6,:); tauC_vgl(6,:)])
% title('Axis 6')
% xlabel('Time (s)')
% ylabel('Torque (Nm)')
% legend('Measurement', 'Model')

tauC_vgl = transpose( reshape( C_c*theta_B, [], 6) );

%Plotten der Drehmomente
figure(2)
sgtitle('Model B for Trajectory C')
for plott = 1:6 
    subplot(3,2,plott)
    plot(timeB, [tauC(plott,:); tauC_vgl(plott,:)])
    title(strcat('Joint ',num2str(plott)))
    xlabel('time (s)')
    ylabel('torque (Nm)')
end
legend('training', 'validation', 'Location','northeast')


rmpath(genpath('../'));