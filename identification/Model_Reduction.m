
% Bjoern Volkmann, bjoern.volkmann@imes.uni-hannover.de, 2020-06
% (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover

addpath(genpath('../'));
clear
close all
clc
tic


%% Parameter

errorTol = 1.05;% Error Tolerance
nAxes = 6; %number of joints
nReibParam = 2; %number of friction parameters

% % % which process is to be used % % %
% process1 = Process 1 //PnP movement
% process2 = Process 2 //square
evaluated_process = 'process1';

keep_friction = true;% wether the friction parameters can be removed from the model or not // true=parameters will be kept -- false=parameters can be removed

% Robot Parameters
g = 9.81; % Gravitation
l = [640;250;630;196;805;100]/1000; % Kinematic
ks = 1000; % Steifigkeit der Coulomb-Reibung (M_Reib = fc * tanh(ks*dq))

           
%% Load Data
fprintf('Load Data \n');
%Trainings Set
if strcmpi(evaluated_process, 'process1')
    load('Identification_Data_Process1.mat');
elseif strcmpi(evaluated_process, 'process2')
    load('Identification_Data_Process2.mat');
else
    fprintf('Process not defined \n');
    return
end

qT = q_meas.average;
dqT = dq_meas.average;
ddqT = ddq_meas.average;
tauT = reshape(tau_meas.average, [], 1);
tauTT = abs(tau_meas.average);
VarianceT = tau_meas.variance;
timeT = time_vec_save;


C_train = RegressorMatrix(ddqT,dqT,g,ks,l,qT);%Designmatrix

%Validierungs Set
if strcmpi(evaluated_process, 'process1')
    load('Validation_Data_Process1.mat');
elseif strcmpi(evaluated_process, 'process2')
    load('Validation_Data_Process2.mat');
else
    fprintf('Process not defined \n');
    return
end

qV = deg2rad(deg_q);
dqV = deg2rad(deg_dq);
ddqV = deg2rad(deg_ddq);
tauV = reshape(tau, [], 1);
tauVV = abs(tau);
timeV = time_vec_save;


C_val = RegressorMatrix(ddqV,dqV,g,ks,l,qV);%Designmatrix

clear deg_q deg_dq deg_ddq tau time_vec_save q_meas dq_meas ddq_meas tau_meas
fprintf('Data Loaded \n');

%% Initialisation of Variables

%Number of parameters
[~, numParameter] = size(C_train);

%Initialization of sets E and I
ignoredIndex = false(1, numParameter); %set E
%the parameter is nonnegligible = true
%the parameter is negligible = false
if keep_friction == true
    ignoredIndex(end-11:end) = true;%friction is nonnegligible
    numParameter_temp = numParameter-12;
else
    numParameter_temp = numParameter;
end

index = 1:numParameter;%set I


%arrays for the model error
eT.ges = zeros(1, numParameter_temp+1);% Index T for the training set
eV.ges = zeros(1, numParameter_temp+1);% Index V for the validation set
eT.axes = zeros(6, numParameter_temp+1);
eV.axes = zeros(6, numParameter_temp+1);
eT.ges_rel = zeros(1, numParameter_temp+1);
eV.ges_rel = zeros(1, numParameter_temp+1);
eT.axes_rel = zeros(6, numParameter_temp+1);
eV.axes_rel = zeros(6, numParameter_temp+1);
theta_red = zeros(numParameter, numParameter_temp);

%matrix to document the success of the model reduction for all axes in each iteration
success = false(nAxes,numParameter_temp);

%array to document the index of the tested parameter in each iteration
testedIndex = zeros(1,numParameter_temp);


%% Parameter Estimation

fprintf('Initial parameter estimation \n');
[rT, ~] = size(tauT);%number of measurement points in the training set
[rV, ~] = size(tauV);%number of measurement points in the validation set

% %inition parameter estimation using the WLS method with the full model
% w = reshape(VarianceT, [], 1);
% w = 1./w;
% WC = zeros(size(C_train));
% for row = 1:rT
%     WC(row,:) = w(row)*C_train(row,:);
% end
% CT_C = transpose(C_train) * WC;% C^T*W*C
% CT_y = transpose(C_train) * (w .* tauT);% C^T*W*y
% 
% rang = rank(CT_C, 1e-8);
% 
% %Test for full rank of C^T*W*C
% if rang < numParameter
%    fprintf('ERROR: The Watrix C^T*W*C is not of full rank \n');
%     return
% end
% 
% %estimate parameters
% theta = CT_C \ CT_y;

theta = WLS_method(C_train, tauT, VarianceT);


%calculate absolute and relative model error e_0 and e^*_0 for the full model
eT.axes(:, 1) = transpose( mean( reshape( abs(tauT - C_train*theta), [], 6)));%absolute
eV.axes(:, 1) = transpose( mean( reshape( abs(tauV - C_val*theta), [], 6)));%absolute

eT.axes_rel(:, 1) = (rT/6) * eT.axes(:, 1) ./ transpose( sum(tauTT) );%relative
eV.axes_rel(:, 1) = (rV/6) * eV.axes(:, 1) ./ transpose( sum(tauVV) );%relative

%averages
eT.ges(1) = sum(eT.axes(:, 1))/6;
eV.ges(1) = sum(eV.axes(:, 1))/6;
eT.ges_rel(1) = sum(eT.axes_rel(:, 1))/6;
eV.ges_rel(1) = sum(eV.axes_rel(:, 1))/6;


%save theta_0
theta_red(:,1) = theta;


%calculate condition number of the full model
Kondition.voll = cond( C_train * diag(1./sqrt(sum(C_train.^2,1))) );


%Retry Zähler
retry = 0;

%% Model Reduction
fprintf('Start model reduction \n');
%Modell iterativ reduzieren
for i = 1:numParameter_temp
    
    %Construct design matrix
    C_train_temp = C_train(:,index);
    
    %%find parameter with lowest sensitivity
    sens = transpose( mean( abs(C_train_temp) ) );
    %ignore nonnegligible parameters
    sens(ignoredIndex) = 2*max(sens);

    %search for minimum
    [~, ind] = min(sens);

    %delete Parameter
    col2keep = index;%bis jetzt beibehaltene Spalten
    col2keep(ind) = [];%identifizierte Spalte entfernen
    C_train_temp = C_train(:, col2keep);
    C_val_temp = C_val(:, col2keep);

    %estimate parameters
%     WC = zeros(size(C_train_temp));
%     for row = 1:rT
%         WC(row,:) = w(row)*C_train_temp(row,:);
%     end
%     CT_C = transpose(C_train_temp) * WC; %C^T*C ist Gram-Matrix von C
%     CT_y = transpose(C_train_temp) * (w .* tauT);
% 
%     theta = CT_C \ CT_y;
    
    theta = WLS_method(C_train_temp, tauT, VarianceT);

    %calculate model errors
    eT.axes(:, i+1) = transpose( mean( reshape( abs(tauT - C_train_temp*theta), [], 6)));
    eV.axes(:, i+1) = transpose( mean( reshape( abs(tauV - C_val_temp*theta), [], 6)));

    eT.axes_rel(:, i+1) = (rT/6) * eT.axes(:, i+1) ./ transpose( sum(tauTT) );
    eV.axes_rel(:, i+1) = (rV/6) * eV.axes(:, i+1) ./ transpose( sum(tauVV) );

    eT.ges(i+1) = sum(eT.axes(:, i+1))/6;
    eV.ges(i+1) = sum(eV.axes(:, i+1))/6;
    eT.ges_rel(i+1) = sum(eT.axes_rel(:, i+1))/6;
    eV.ges_rel(i+1) = sum(eV.axes_rel(:, i+1))/6;

    fprintf('Iteration: %d == Remove colomn %d with sensitivity %1.4f :', i, index(ind), sens(ind));
    testedIndex(i) = index(ind);%Abspeichern welcher index in der Iteration überprüft wurde

    %test for successful model reduction
    success(1,i) = eV.axes_rel(1,i+1) / eV.axes_rel(1,1) <= errorTol;
    success(2,i) = eV.axes_rel(2,i+1) / eV.axes_rel(2,1) <= errorTol;
    success(3,i) = eV.axes_rel(3,i+1) / eV.axes_rel(3,1) <= errorTol;
    success(4,i) = eV.axes_rel(4,i+1) / eV.axes_rel(4,1) <= errorTol;
    success(5,i) = eV.axes_rel(5,i+1) / eV.axes_rel(5,1) <= errorTol;
    success(6,i) = eV.axes_rel(6,i+1) / eV.axes_rel(6,1) <= errorTol;
    if any(~success(:,i)) %not successful
        
        retry = retry +1;
        fprintf('failure \\ %d\n', retry)
        ignoredIndex(ind) = true;
        
    else %successful
        
        %Wenn rel. Fehler die Grenze noch nicht überschritten hat,
        %entfernte Spalte merken
        fprintf('success \\ %d\n', retry)
        index(ind) = [];
        ignoredIndex(ind) = [];
        theta_red(index, i+1-retry) = theta;%Save estimated parameters on succsessful model reduction
        
    end

    
end

fprintf('Model reduction complete \n');
%remove unused entries
theta_red(:, i+2-retry:end) = [];

clear ind


%% Calulate Modell Prediction with Reduced Model


fprintf('%d parameters removed\n', numParameter-length(index))

fprintf('\tNOT removed columns:\t')
fprintf('%d \t', index)
fprintf('\n')

%calculate torque
tauT_vgl = C_train(:,index)*theta_red(index,end);%Trainingsset
tauV_vgl = C_val(:,index)*theta_red(index,end);%Validierungsset

%Condition number
Kondition.red = cond( C_train(:,index) * diag(1./sqrt(sum(C_train(:,index).^2,1))) );

fprintf('Condition number of the full model: %1.2f \n', Kondition.voll)
fprintf('Condition number of the reduced model: %1.2f \n', Kondition.red)

%berechnen von relativem und absolutem Fehler
eT.result = transpose( mean( reshape( abs(tauT - tauT_vgl), [], 6)));
eV.result = transpose( mean( reshape( abs(tauV - tauV_vgl), [], 6)));
eT.result_rel = (rT/6) * eT.result ./ transpose( sum(tauTT) );
eV.result_rel = (rV/6) * eV.result ./ transpose( sum(tauVV) );


%reshape zum späteren plotten
tauT_vgl = transpose(reshape(tauT_vgl, [], 6));
tauV_vgl = transpose(reshape(tauV_vgl, [], 6));
tauT = transpose(reshape(tauT, [], 6));
tauV = transpose(reshape(tauV, [], 6));


fprintf('Results for theta:\n')
%===Printe Sensitivität
fprintf('\tSensitivity:\t\t')
fprintf('%f \t', mean( abs(C_train) ))
fprintf('\n')
fprintf('\n')


%Volles Modell
fprintf('\ttheta (full):\t\t')
fprintf('%f \t', theta_red(:,1))
fprintf('\n')


%reduziert
fprintf('\ttheta (reduced):\t')
fprintf('%f \t', theta_red(:,end))
fprintf('\n')


%save results
save(strcat('ModelReduction_', evaluated_process),'eT', 'eV', 'theta_red');

%% Plotten der Ergebnisse
fprintf('Plotting results ... \n');


%Plot der Reduktion mit relativem Fehler gesamt
figure(1)
plot(0:i ,eT.ges_rel ,'DisplayName','Training-Error')
title('Average relativ model error')
hold on
plot(0:i ,eV.ges_rel ,'DisplayName','Validation-Error')

legend('Location','northeast')
ymax = ylim;
xlabel('Iterationen')
ylabel('Durchschnittlicher relativer Approximationsfehler')
hold off


%Plot des absoluten Fehlerverlaufs pro Achse
figure(2)
sgtitle('Average Absolute Model Error')
for plott = 1:6 
    subplot(3,2,plott)
    plot(0:i, eT.axes(plott,:))
    title(strcat('Joint ',num2str(plott)))
    hold on
    plot(0:i, eV.axes(plott,:))
    xlabel('iteration $k$','interpreter','latex')
    ylabel('Absolute Model Error (Nm)')
    hold off
end
legend('training', 'validation', 'Location','northeast')

% subplot(3,2,2)
% plot(0:i, eT.axes(2,:))%, 'DisplayName','training')
% title('Joint 2')
% hold on
% plot(0:i, eV.axes(2,:))%, 'DisplayName','Validierung')
% % legend('Location','northeast')
% xlabel('iteration $k$','interpreter','latex')
% ylabel('Absolute Model Error (Nm)')
% hold off
% 
% subplot(3,2,3)
% plot(0:i, eT.axes(3,:))%, 'DisplayName','Training')
% title('Joint 3')
% hold on
% plot(0:i, eV.axes(3,:))%, 'DisplayName','Validierung')
% % legend('Location','northeast')
% xlabel('iteration $k$','interpreter','latex')
% ylabel('Absolute Model Error (Nm)')
% hold off
% 
% subplot(3,2,4)
% plot(0:i, eT.axes(4,:))%, 'DisplayName','Training')
% title('Joint 4')
% hold on
% plot(0:i, eV.axes(4,:))%, 'DisplayName','Validierung')
% % legend('Location','northeast')
% xlabel('iteration $k$','interpreter','latex')
% ylabel('Absolute Model Error (Nm)')
% hold off
% 
% subplot(3,2,5)
% plot(0:i, eT.axes(5,:))%, 'DisplayName','Training')
% title('Joint 5')
% hold on
% plot(0:i, eV.axes(5,:))%, 'DisplayName','Validierung')
% % legend('Location','northeast')
% xlabel('iteration $k$','interpreter','latex')
% ylabel('Absolute Model Error (Nm)')
% hold off
% 
% subplot(3,2,6)
% plot(0:i, eT.axes(6,:))%, 'DisplayName','Training')
% title('Joint 6')
% hold on
% plot(0:i, eV.axes(6,:))%, 'DisplayName','Validierung')
% % legend('Location','northeast')
% xlabel('iteration $k$','interpreter','latex')
% ylabel('Absolute Model Error (Nm)')
% hold off


%Plot des relativen Fehlerverlaufs pro Achse
figure(3)
sgtitle('Relative Model Error')
for plott = 1:6 
    subplot(3,2,plott)
    plot(0:i, eT.axes_rel(plott,:))
    title(strcat('Joint ',num2str(plott)))
    hold on
    plot(0:i, eV.axes_rel(plott,:))
    xlabel('iteration $k$','interpreter','latex')
    ylabel('Relative Model Error (Nm)')
    hold off
end
legend('training', 'validation', 'Location','northeast')

% subplot(3,2,2)
% plot(0:i, eT.axes_rel(2,:), 'DisplayName','Training')
% title('Durchschnittlicher relativer Modellfehler Achse 2')
% hold on
% plot(0:i, eV.axes_rel(2,:), 'DisplayName','Validierung')
% legend('Location','northeast')
% xlabel('Iterationen')
% ylabel('Durchschnittlicher Approximationsfehler (Nm)')
% hold off
% 
% subplot(3,2,3)
% plot(0:i, eT.axes_rel(3,:), 'DisplayName','Training')
% title('Durchschnittlicher relativer Modellfehler Achse 3')
% hold on
% plot(0:i, eV.axes_rel(3,:), 'DisplayName','Validierung')
% legend('Location','northeast')
% xlabel('Iterationen')
% ylabel('Durchschnittlicher Approximationsfehler (Nm)')
% hold off
% 
% subplot(3,2,4)
% plot(0:i, eT.axes_rel(4,:), 'DisplayName','Training')
% title('Durchschnittlicher relativer Modellfehler Achse 4')
% hold on
% plot(0:i, eV.axes_rel(4,:), 'DisplayName','Validierung')
% legend('Location','northeast')
% xlabel('Iterationen')
% ylabel('Durchschnittlicher Approximationsfehler (Nm)')
% hold off
% 
% subplot(3,2,5)
% plot(0:i, eT.axes_rel(5,:), 'DisplayName','Training')
% title('Durchschnittlicher relativer Modellfehler Achse 5')
% hold on
% plot(0:i, eV.axes_rel(5,:), 'DisplayName','Validierung')
% legend('Location','northeast')
% xlabel('Iterationen')
% ylabel('Durchschnittlicher Approximationsfehler (Nm)')
% hold off
% 
% subplot(3,2,6)
% plot(0:i, eT.axes_rel(6,:), 'DisplayName','Training')
% title('Durchschnittlicher relativer Modellfehler Achse 6')
% hold on
% plot(0:i, eV.axes_rel(6,:), 'DisplayName','Validierung')
% legend('Location','northeast')
% xlabel('Iterationen')
% ylabel('Durchschnittlicher Approximationsfehler (Nm)')
% hold off


%Plotten der Drehmomente
figure(4)
sgtitle('Model Prediction')
for plott = 1:6 
    subplot(3,2,plott)
    plot(timeV, [tauV(plott,:); tauV_vgl(plott,:)])
    title(strcat('Joint ',num2str(plott)))
    xlabel('time (s)')
    ylabel('torque (Nm)')
end
legend('measurement', 'model')

% subplot(3,2,2)
% plot(timeV, [tauV(2,:); tauV_vgl(2,:)])
% title('Vergleich Messung/LS Achse 2')
% xlabel('Zeit (s)')
% ylabel('Drehmoment (Nm)')
% legend('Messung', 'LS')
% 
% subplot(3,2,3)
% plot(timeV, [tauV(3,:); tauV_vgl(3,:)])
% title('Vergleich Messung/LS Achse 3')
% xlabel('Zeit (s)')
% ylabel('Drehmoment (Nm)')
% legend('Messung', 'LS')
% 
% subplot(3,2,4)
% plot(timeV, [tauV(4,:); tauV_vgl(4,:)])
% title('Vergleich Messung/LS Achse 4')
% xlabel('Zeit (s)')
% ylabel('Drehmoment (Nm)')
% legend('Messung', 'LS')
% 
% subplot(3,2,5)
% plot(timeV, [tauV(5,:); tauV_vgl(5,:)])
% title('Vergleich Messung/LS Achse 5')
% xlabel('Zeit (s)')
% ylabel('Drehmoment (Nm)')
% legend('Messung', 'LS')
% 
% subplot(3,2,6)
% plot(timeV, [tauV(6,:); tauV_vgl(6,:)])
% title('Vergleich Messung/LS Achse 6')
% xlabel('Zeit (s)')
% ylabel('Drehmoment (Nm)')
% legend('Messung', 'LS')

toc
rmpath(genpath('../'));