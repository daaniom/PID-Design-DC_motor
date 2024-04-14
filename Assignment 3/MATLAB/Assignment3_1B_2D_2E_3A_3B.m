%%
%% 
% This file consist of the calculations for questions 1B, 2D, 2E, 3A and 3B
% Other questions that needed calculations have a seperate file
clc;
clear all;
close all;
figures=true; % Figures for 1B,2D,2E and 3A
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Question 1B %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts=0.01;
% Values such that system becomes unstable z = TsK+1 (pole must lie in unit circle)
k_up = (-1-1)/-Ts;
k_low = (1-1)/-Ts;
if figures == true
    figure(1)
    hold on
    x = [0.8 0.7];   % adjust length and location of arrow 
    y = [0.5 0.5];
    annotation('textarrow', x, y, 'FontSize', 13, 'Linewidth', 2)
    annotation('textbox', [0.68 0.3 0.7 0.27], 'EdgeColor', 'none', 'String', 'Increasing K', 'FontSize', 12, 'Linewidth', 2)
    K=k_low+1
    a=tf([-Ts],[1,K*Ts-1],Ts);
    hold on;
    pzmap(a)
    hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
    hm(2).MarkerSize = 25;                      % 3Zero5 Marker
    hm(3).MarkerSize = 25;                      % 3Pole5 Marker
    K= (k_up-k_low)/2 % between upper and lower
    b=tf([-Ts],[1,K*Ts-1],Ts);
    pzmap(b)
    hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
    hm(2).MarkerSize = 25;                      % 0Zero5 Marker
    hm(3).MarkerSize = 25;                      % 3Pole5 Marker
    K= (k_up)-1
    c=tf([-Ts],[1,K*Ts-1],Ts);
    pzmap(c)
    hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
    hm(2).MarkerSize = 25;                      % 3Zero2 Marker
    hm(3).MarkerSize = 25;                      % 3Pole5 Marker
    %quiver(p1(1),p1(2),dp(1),dp(2),0)
    hold off
end

k_values = [(k_low)+1:k_up-1];
poles_DT = -Ts*k_values+1;
poles_CT = -log(poles_DT)/Ts;

if figures==true 
    figure(2)
    plot(k_values,poles_CT,'LineWidth', 2, 'Color', 'k');
    grid on;
    xlabel('K values')
    ylabel('\omega [rad/s]')
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Question 2D %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 1;
B = Ts;
C = -1;
D = 0;
Q_R = logspace(-6,6,1000);
R = 1;
Q = Q_R*R;

L_expression = -1*(Q_R+sqrt(Q_R.^2+4*Q_R))./(Q_R+sqrt(Q_R.^2+4*Q_R)+2);
L_LQE = zeros(1,length(Q_R));
for i = 1:length(Q_R)
    L_LQE(i) = dlqr(A',(A'*C'),Q(i),R)';
end
error_L = L_expression-L_LQE;

if figures==true 
    figure(3)
    subplot(1,2,1)
    semilogx(Q_R,L_expression,'LineWidth', 2, 'Color', 'k');
    hold on;
    semilogx(Q_R,L_LQE,'LineWidth', 2, 'Color', '[0.7, 0.7, 0.7]');
    grid on;
    xlabel('Q/R')
    ylabel('L')
    legend('L_{\infty}','L_{LQE}')
    hold off;
    subplot(1,2,2)
    semilogx(Q_R,error_L,'LineWidth', 2, 'Color', 'k');
    grid on;
    xlabel('Q/R')
    ylabel('Error = L_{\infty}- L_{LQE}')
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Question 2E %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_LQE = 1 -1*(Q_R+sqrt(Q_R.^2+4*Q_R))./(Q_R+sqrt(Q_R.^2+4*Q_R)+2); % discrete pole
s_LQE = -log(z_LQE)/Ts; % continuous pole
if figures == true
    figure(4)
    subplot(1,2,1)
    semilogx(Q_R,z_LQE,'LineWidth', 2, 'Color', 'k');
    xlabel('Q/R')
    ylabel('Closed-loop LQE pole in Discrete time')
    grid on;
    subplot(1,2,2)
    semilogx(Q_R,s_LQE,'LineWidth', 2, 'Color', 'k');
    xlabel('Q/R')
    ylabel('Closed-loop LQE pole in Continuous time')
    grid on;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Question 3A %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the relative path to the CSV file
relativePath_distance_10cm = 'distance_10cm.csv';
relativePath_distance_15cm = 'distance_15cm.csv';
relativePath_distance_20cm = 'distance_20cm.csv';
relativePath_distance_25cm = 'distance_25cm.csv';

% Get the current working directory
currentDir = pwd;
currentDir_Assignment3_3A_Data = fullfile(currentDir,'Assignment3_3A_data');
% Construct the full file path by combining the current directory and the relative path
csvfile_10cm = fullfile(currentDir_Assignment3_3A_Data, relativePath_distance_10cm);
csvfile_15cm = fullfile(currentDir_Assignment3_3A_Data, relativePath_distance_15cm);
csvfile_20cm = fullfile(currentDir_Assignment3_3A_Data, relativePath_distance_20cm);
csvfile_25cm = fullfile(currentDir_Assignment3_3A_Data, relativePath_distance_25cm);
%csvfile = 'C:\Users\toonm\OneDrive\Documenten\KULeuven\2023-2024\Sem1\Systems and Control Theory\Control theory Assigment\Assignment 1\excitationSignal_pulseTrain.csv';
labels = strsplit(fileread(csvfile_10cm), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
measurements_10cm = dlmread(csvfile_10cm, ',', 2, 0); % Data follows the labels
measurements_15cm = dlmread(csvfile_15cm, ',', 2, 0); % Data follows the labels
measurements_20cm = dlmread(csvfile_20cm, ',', 2, 0); % Data follows the labels
measurements_25cm = dlmread(csvfile_25cm, ',', 2, 0); % Data follows the labels

% 10cm
sensor_measurements_10cm = measurements_10cm(:,2).';
actual_distances_10cm = 0.1*ones(1,length(sensor_measurements_10cm));
covariance_matrix_10cm = cov(sensor_measurements_10cm, actual_distances_10cm);
covariance_10cm = covariance_matrix_10cm(1,1);
% 15cm
sensor_measurements_15cm = measurements_15cm(:,2).';
actual_distances_15cm = 0.1*ones(1,length(sensor_measurements_15cm));
covariance_matrix_15cm = cov(sensor_measurements_15cm, actual_distances_15cm);
covariance_15cm = covariance_matrix_15cm(1,1);
% 20cm
sensor_measurements_20cm = measurements_20cm(:,2).';
actual_distances_20cm = 0.1*ones(1,length(sensor_measurements_20cm));
covariance_matrix_20cm = cov(sensor_measurements_20cm, actual_distances_20cm);
covariance_20cm = covariance_matrix_20cm(1,1);
% 25cm
sensor_measurements_25cm = measurements_25cm(:,2).';
actual_distances_25cm = 0.1*ones(1,length(sensor_measurements_25cm));
covariance_matrix_25cm = cov(sensor_measurements_25cm, actual_distances_25cm);
covariance_25cm = covariance_matrix_25cm(1,1);

mean_covariance = (covariance_10cm+covariance_15cm+covariance_20cm+covariance_25cm)/4

K_estimation = dlqr(1,0.01,8*10^-6,8*10^-6)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Question 3B %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentDir_Assignment3_3B_Data = fullfile(currentDir,'Assignment3_3B_data');
endtime = 3.5;    % this time becomes smaller because for higher values of K the time measured becomes smaller due to the instabillity of these values
time = linspace(0,endtime,endtime/0.01+1);
N = length(time);
%%%% 3B_k1 %%%%
csvfile_3B_k1 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k1.csv');
measurements_3B_k1 = dlmread(csvfile_3B_k1, ',', 2, 0); % Data follows the labels

referenceStep_k1 = measurements_3B_k1(:,2);
indexFirst_k1 = find(referenceStep_k1~=-0.1,1,"first")-50;
indexLast_k1 = indexFirst_k1+N-1;
referenceStep_k1 = referenceStep_k1(indexFirst_k1:indexLast_k1);
responseStep_k1 = -1*measurements_3B_k1(indexFirst_k1:indexLast_k1,3);  % measured = -x state ==> state x = - measured

voltageA_k1 = measurements_3B_k1(indexFirst_k1:indexLast_k1,4);
voltageB_k1 = measurements_3B_k1(indexFirst_k1:indexLast_k1,5);
%%%% 3B_k2 %%%%
csvfile_3B_k2 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k2.csv');
measurements_3B_k2 = dlmread(csvfile_3B_k2, ',', 2, 0); % Data follows the labels

referenceStep_k2 = measurements_3B_k2(:,2);
indexFirst_k2 = find(referenceStep_k2~=-0.1,1,"first")-50;
indexLast_k2 = indexFirst_k2+N-1;
referenceStep_k2 = referenceStep_k2(indexFirst_k2:indexLast_k2);
responseStep_k2 = -1*measurements_3B_k2(indexFirst_k2:indexLast_k2,3);  % measured = -x state ==> state x = - measured

voltageA_k2 = measurements_3B_k2(indexFirst_k2:indexLast_k2,4);
voltageB_k2 = measurements_3B_k2(indexFirst_k2:indexLast_k2,5);
%%%% 3B_k5 %%%%
csvfile_3B_k5 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k5.csv');
measurements_3B_k5 = dlmread(csvfile_3B_k5, ',', 2, 0); % Data follows the labels

referenceStep_k5 = measurements_3B_k5(:,2);
indexFirst_k5 = find(referenceStep_k5~=-0.1,1,"first")-50;
indexLast_k5 = indexFirst_k5+N-1;
referenceStep_k5 = referenceStep_k5(indexFirst_k5:indexLast_k5);
responseStep_k5 = -1*measurements_3B_k5(indexFirst_k5:indexLast_k5,3);  % measured = -x state ==> state x = - measured

voltageA_k5 = measurements_3B_k5(indexFirst_k5:indexLast_k5,4);
voltageB_k5 = measurements_3B_k5(indexFirst_k5:indexLast_k5,5);
%%%% 3B_k10 %%%%
csvfile_3B_k10 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k10.csv');
measurements_3B_k10 = dlmread(csvfile_3B_k10, ',', 2, 0); % Data follows the labels

referenceStep_k10 = measurements_3B_k10(:,2);
indexFirst_k10 = find(referenceStep_k10~=-0.1,1,"first")-50;
indexLast_k10 = indexFirst_k10+N-1;
referenceStep_k10 = referenceStep_k10(indexFirst_k10:indexLast_k10);
responseStep_k10 = -1*measurements_3B_k10(indexFirst_k10:indexLast_k10,3);  % measured = -x state ==> state x = - measured

voltageA_k10 = measurements_3B_k10(indexFirst_k10:indexLast_k10,4);
voltageB_k10 = measurements_3B_k10(indexFirst_k10:indexLast_k10,5);
%%%% 3B_k25 %%%%
csvfile_3B_k25 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k25.csv');
measurements_3B_k25 = dlmread(csvfile_3B_k25, ',', 2, 0); % Data follows the labels

referenceStep_k25 = measurements_3B_k25(:,2);
indexFirst_k25 = find(referenceStep_k25~=-0.1,1,"first")-50;
indexLast_k25 = indexFirst_k25+N-1;
referenceStep_k25 = referenceStep_k25(indexFirst_k25:indexLast_k25);
responseStep_k25 = -1*measurements_3B_k25(indexFirst_k25:indexLast_k25,3);  % measured = -x state ==> state x = - measured

voltageA_k25 = measurements_3B_k25(indexFirst_k25:indexLast_k25,4);
voltageB_k25 = measurements_3B_k25(indexFirst_k25:indexLast_k25,5);
%%%% 3B_k50 %%%%
csvfile_3B_k50 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k50.csv');
measurements_3B_k50 = dlmread(csvfile_3B_k50, ',', 2, 0); % Data follows the labels

referenceStep_k50 = measurements_3B_k50(:,2);
indexFirst_k50 = find(referenceStep_k50~=-0.1,1,"first")-50;
indexLast_k50 = indexFirst_k50+N-1;
referenceStep_k50 = referenceStep_k50(indexFirst_k50:indexLast_k50);
responseStep_k50 = -1*measurements_3B_k50(indexFirst_k50:indexLast_k50,3);  % measured = -x state ==> state x = - measured

voltageA_k50 = measurements_3B_k50(indexFirst_k50:indexLast_k50,4);
voltageB_k50 = measurements_3B_k50(indexFirst_k50:indexLast_k50,5);
%%%% 3B_k75 %%%%
csvfile_3B_k75 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k75.csv');
measurements_3B_k75 = dlmread(csvfile_3B_k75, ',', 2, 0); % Data follows the labels

referenceStep_k75 = measurements_3B_k75(:,2);
indexFirst_k75 = find(referenceStep_k75~=-0.1,1,"first")-50;
indexLast_k75 = indexFirst_k75+N-1;
referenceStep_k75 = referenceStep_k75(indexFirst_k75:indexLast_k75);
responseStep_k75 = -1*measurements_3B_k75(indexFirst_k75:indexLast_k75,3);  % measured = -x state ==> state x = - measured

voltageA_k75 = measurements_3B_k75(indexFirst_k75:indexLast_k75,4);
voltageB_k75 = measurements_3B_k75(indexFirst_k75:indexLast_k75,5);
%%%% 3B_k100 %%%%
csvfile_3B_k100 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k100.csv');
measurements_3B_k100 = dlmread(csvfile_3B_k100, ',', 2, 0); % Data follows the labels

referenceStep_k100 = measurements_3B_k100(:,2);
indexFirst_k100 = find(referenceStep_k100~=-0.1,1,"first")-50;
indexLast_k100 = indexFirst_k100+N-1;
referenceStep_k100 = referenceStep_k100(indexFirst_k100:indexLast_k100);
responseStep_k100 = -1*measurements_3B_k100(indexFirst_k100:indexLast_k100,3);  % measured = -x state ==> state x = - measured

voltageA_k100 = measurements_3B_k100(indexFirst_k100:indexLast_k100,4);
voltageB_k100 = measurements_3B_k100(indexFirst_k100:indexLast_k100,5);
%%%% 3B_k125 %%%%
csvfile_3B_k125 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k125.csv');
measurements_3B_k125 = dlmread(csvfile_3B_k125, ',', 2, 0); % Data follows the labels

referenceStep_k125 = measurements_3B_k125(:,2);
indexFirst_k125 = find(referenceStep_k125~=-0.1,1,"first")-50;
indexLast_k125 = indexFirst_k125+N-1;
referenceStep_k125 = referenceStep_k125(indexFirst_k125:indexLast_k125);
responseStep_k125 = -1*measurements_3B_k125(indexFirst_k125:indexLast_k125,3);  % measured = -x state ==> state x = - measured

voltageA_k125 = measurements_3B_k125(indexFirst_k125:indexLast_k125,4);
voltageB_k125 = measurements_3B_k125(indexFirst_k125:indexLast_k125,5);
%%%% 3B_k150 %%%%
csvfile_3B_k150 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k150.csv');
measurements_3B_k150 = dlmread(csvfile_3B_k150, ',', 2, 0); % Data follows the labels

referenceStep_k150 = measurements_3B_k150(:,2);
indexFirst_k150 = find(referenceStep_k150~=-0.1,1,"first")-50;
indexLast_k150 = indexFirst_k150+N-1;
referenceStep_k150 = referenceStep_k150(indexFirst_k150:indexLast_k150);
responseStep_k150 = -1*measurements_3B_k150(indexFirst_k150:indexLast_k150,3);  % measured = -x state ==> state x = - measured

voltageA_k150 = measurements_3B_k150(indexFirst_k150:indexLast_k150,4);
voltageB_k150 = measurements_3B_k150(indexFirst_k150:indexLast_k150,5);
%%%% 3B_k175 %%%%
csvfile_3B_k175 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k175.csv');
measurements_3B_k175 = dlmread(csvfile_3B_k175, ',', 2, 0); % Data follows the labels

referenceStep_k175 = measurements_3B_k175(:,2);
indexFirst_k175 = find(referenceStep_k175~=-0.1,1,"first")-50;
indexLast_k175 = indexFirst_k175+N-1;
referenceStep_k175 = referenceStep_k175(indexFirst_k175:indexLast_k175);
responseStep_k175 = -1*measurements_3B_k175(indexFirst_k175:indexLast_k175,3);  % measured = -x state ==> state x = - measured

voltageA_k175 = measurements_3B_k175(indexFirst_k175:indexLast_k175,4);
voltageB_k175 = measurements_3B_k175(indexFirst_k175:indexLast_k175,5);
%%%% 3B_k199 %%%%
csvfile_3B_k199 = fullfile(currentDir_Assignment3_3B_Data,'Assignment3_3b_k199.csv');
measurements_3B_k199 = dlmread(csvfile_3B_k199, ',', 2, 0); % Data follows the labels

referenceStep_k199 = measurements_3B_k199(:,2);
indexFirst_k199 = find(referenceStep_k199~=-0.1,1,"first")-50;
indexLast_k199 = indexFirst_k199+N-1;
referenceStep_k199 = referenceStep_k199(indexFirst_k199:indexLast_k199);
responseStep_k199 = -1*measurements_3B_k199(indexFirst_k199:indexLast_k199,3);  % measured = -x state ==> state x = - measured

voltageA_k199 = measurements_3B_k199(indexFirst_k199:indexLast_k199,4);
voltageB_k199 = measurements_3B_k199(indexFirst_k199:indexLast_k199,5);
%%%%%%%%%% plot results %%%%%%%%%%
figure(4) % all values for k
hold on;
plot(time,referenceStep_k1,'LineWidth',2,'Color','k')
plot(time,responseStep_k1,'LineWidth',2)
plot(time,responseStep_k2,'LineWidth',2)
plot(time,responseStep_k5,'LineWidth',2)
plot(time,responseStep_k10,'LineWidth',2)
plot(time,responseStep_k25,'LineWidth',2)
plot(time,responseStep_k50,'LineWidth',2)
plot(time,responseStep_k75,'LineWidth',2)
plot(time,responseStep_k100,'LineWidth',2,'Color', [0, 0.5, 0.5])
plot(time,responseStep_k125,'LineWidth',2,'Color','m')
plot(time,responseStep_k150,'LineWidth',2,'Color',[0.5, 0.5, 0.5])
plot(time,responseStep_k175,'LineWidth',2,'Color',[1, 0, 0])
plot(time,responseStep_k199,'LineWidth',2,'Color',[1,1,0])
legend('Step: 0.1m \rightarrow 0.25m','K=1','K=2','K=5','K=10','K=25','K=50','K=75','K=100','K=125','K=150','K=175','K=199','FontSize',10);
xlabel('Time [s]','FontSize',15)
ylabel('Position [m]','FontSize',15)
grid on;
hold off;

figure(5)  %selected values for K
hold on;
plot(time,referenceStep_k1,'LineWidth',2,'Color','k')
plot(time,responseStep_k1,'LineWidth',2)
plot(time,responseStep_k2,'LineWidth',2)
plot(time,responseStep_k25,'LineWidth',2)
plot(time,responseStep_k100,'LineWidth',2)
plot(time,responseStep_k125,'LineWidth',2)
plot(time,responseStep_k199,'LineWidth',2)
legend('Step: 0.1m \rightarrow 0.25m','K=1','K=2','K=25','K=100','K=125','K=199','FontSize',15);
xlabel('Time [s]','FontSize',15)
ylabel('Position [m]','FontSize',15)
grid on;
hold off;

figure(6) %selected values for K
subplot(1,2,1)
% motor A
hold on;
plot(time,voltageA_k1,'LineWidth',2)
plot(time,voltageA_k2,'LineWidth',2)
plot(time,voltageA_k25,'LineWidth',2)
plot(time,voltageA_k100,'LineWidth',2)
plot(time,voltageA_k125,'LineWidth',2)
plot(time,voltageA_k199,'LineWidth',2)
legend('K=1','K=2','K=25','K=100','K=125','K=199','FontSize',15);
xlabel('Time [s]','FontSize',15)
ylabel('Voltage motor A [V]','FontSize',15)
grid on;
hold off;
% motor B
subplot(1,2,2)
hold on;
plot(time,voltageB_k1,'LineWidth',2)
plot(time,voltageB_k2,'LineWidth',2)
plot(time,voltageB_k25,'LineWidth',2)
plot(time,voltageB_k100,'LineWidth',2)
plot(time,voltageB_k125,'LineWidth',2)
plot(time,voltageB_k199,'LineWidth',2)
legend('K=1','K=2','K=25','K=100','K=125','K=199','FontSize',15);
xlabel('Time [s]','FontSize',15)
ylabel('Voltage motor B [V]','FontSize',15)
grid on;
hold off;

figure(7)
hold on;
plot(time,voltageA_k1,'LineWidth',2)
plot(time,voltageA_k2,'LineWidth',2)
plot(time,voltageA_k25,'LineWidth',2)
plot(time,voltageA_k100,'LineWidth',2)
plot(time,voltageA_k125,'LineWidth',2)
plot(time,voltageA_k199,'LineWidth',2)
legend('K=1','K=2','K=25','K=100','K=125','K=199','FontSize',15);
xlabel('Time [s]','FontSize',15)
ylabel('Voltage motor A [V]','FontSize',15)
grid on;
hold off;

figure(8)
hold on;
plot(time,voltageA_k1,'LineWidth',2)
plot(time,voltageA_k2,'LineWidth',2)
plot(time,voltageA_k5,'LineWidth',2)
plot(time,voltageA_k10,'LineWidth',2)
plot(time,voltageA_k25,'LineWidth',2)
legend('K=1','K=2','K=5','K=10','K=25','FontSize',15);
xlabel('Time [s]','FontSize',15)
ylabel('Voltage motor A [V]','FontSize',15)
grid on;
hold off;
