close all;
clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Identification of DC motor, while cart in the air %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the excitation signal
% Define the relative path to the CSV file
relativePath_6V = 'excitationSignal_pulseTrain_6V.csv';
% Get the current working directory
currentDir = pwd;
figures =0;
% Construct the full file path by combining the current directory and the relative path
csvfile = fullfile(currentDir, relativePath_6V);
%csvfile = 'C:\Users\toonm\OneDrive\Documenten\KULeuven\2023-2024\Sem1\Systems and Control Theory\Control theory Assigment\Assignment 1\excitationSignal_pulseTrain.csv';
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels

indexFirst = find(data(:,3)~=0,1,"first")-50;% we started with zero state in Arduino
indexLast = find(data(:,3)~=0,1,"last")+50; % so it's clear we want from -6V back to 0V
uA = data(indexFirst:indexLast,3);
omegaA = data(indexFirst:indexLast,4);
thetaA = data(indexFirst:indexLast,5);

uB = data(indexFirst:indexLast,6);
omegaB = data(indexFirst:indexLast,7);
thetaB = data(indexFirst:indexLast,8);
% a) applied excitation signal
N = length(uA);
fs = 100;
Ts = 1/fs;
t = (0:N-1)'*Ts;
num_periods = 5;
pointsPerPeriod = (N)/num_periods;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if figures==1
figure(1)
% motor A
hold on;
subplot(2,2,1);
plot(t,uA,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('Voltage_{A} [V]', FontSize= 15);
subplot(2,2,3);
plot(t,omegaA,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\omega_{A} [rad/s]', FontSize= 15);
% motor B
subplot(2,2,2);
plot(t,uB,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('Voltage_{B} [V]', FontSize= 15);
subplot(2,2,4);
plot(t,omegaB,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\omega_{B} [rad/s]', FontSize= 15);
end

% to get the overlap the plot for different periods to appreciate the noise
% separate the different periods (one per column)
uA_matrix = reshape(uA,pointsPerPeriod,num_periods);
omegaA_matrix = reshape(omegaA,pointsPerPeriod,num_periods);
thetaA_matrix = reshape(thetaA,pointsPerPeriod,num_periods);

uB_matrix = reshape(uB,pointsPerPeriod,num_periods);
omegaB_matrix = reshape(omegaB,pointsPerPeriod,num_periods);
thetaB_matrix = reshape(thetaB,pointsPerPeriod,num_periods);

% lets compute the mean of the signals across the periods to have a point of comparison to assess the noise
uA_mean = mean(uA_matrix,2);
omegaA_mean = mean(omegaA_matrix,2);
thetaA_mean = mean(thetaA_matrix,2);

uB_mean = mean(uB_matrix,2);
omegaB_mean = mean(omegaB_matrix,2);
thetaB_mean = mean(thetaB_matrix,2);

duA_matrix = uA_matrix - repmat(uA_mean,1,num_periods);
domegaA_matrix = omegaA_matrix - repmat(omegaA_mean,1,num_periods);
dthetaA_matrix = thetaA_matrix - repmat(thetaA_mean,1,num_periods);

duB_matrix = uB_matrix - repmat(uB_mean,1,num_periods);
domegaB_matrix = omegaB_matrix - repmat(omegaB_mean,1,num_periods);
dthetaB_matrix = thetaB_matrix - repmat(thetaB_mean,1,num_periods);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for motor A!
if figures==1
figure(2)
hold on;
subplot(3,2,1);
plot(t(1:pointsPerPeriod),uA_matrix,'LineWidth',2);
grid on;
axis([0 t(pointsPerPeriod) min(uA) max(uA)]);
xlabel('Time [s]',FontSize= 15);
ylabel('Voltage_{a} [V]', FontSize= 15);

subplot(3,2,2);
plot(t(1:pointsPerPeriod),duA_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\DeltaV_{a} [V]', FontSize= 15);

subplot(3,2,3);
plot(t(1:pointsPerPeriod),omegaA_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\omega_{a} [rad/s]', FontSize= 15);

subplot(3,2,4);
plot(t(1:pointsPerPeriod),domegaA_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\Delta\omega_{a} [rad/s]', FontSize= 15);

subplot(3,2,5);
plot(t(1:pointsPerPeriod),thetaA_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\theta_{a} [rad]', FontSize= 15);

subplot(3,2,6);
plot(t(1:pointsPerPeriod),dthetaA_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\Delta\theta_{a} [rad]', FontSize= 15);

% this is for motor B!
figure(3)
hold on;
subplot(3,2,1);
plot(t(1:pointsPerPeriod),uB_matrix,'LineWidth',2);
grid on;
axis([0 t(pointsPerPeriod) min(uB) max(uB)]);
xlabel('Time [s]',FontSize= 15);
ylabel('Voltage_{B} [V]', FontSize= 15);

subplot(3,2,2);
plot(t(1:pointsPerPeriod),duB_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\DeltaV_{B} [V]', FontSize= 15);

subplot(3,2,3);
plot(t(1:pointsPerPeriod),omegaB_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\omega_{B} [rad/s]', FontSize= 15);

subplot(3,2,4);
plot(t(1:pointsPerPeriod),domegaB_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\Delta\omega_{B} [rad/s]', FontSize= 15);

subplot(3,2,5);
plot(t(1:pointsPerPeriod),thetaB_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\theta_{B} [rad]', FontSize= 15);

subplot(3,2,6);
plot(t(1:pointsPerPeriod),dthetaB_matrix,'LineWidth',2);
grid on;
axis tight;
xlabel('Time [s]',FontSize= 15);
ylabel('\Delta\theta_{B} [rad]', FontSize= 15);
end
%% b) Obtaining system parameters
% the empirical transfer-function estimate
% for robustness we will use the mean vectors x_mean and F_mean
N = numel(uA_mean); % number of samples
f = [0:N-1]'*(fs/N); % arrays of frequencies

% compute empiric fourier transform of the signals
% Motor A
UA_f =  fft(uA_mean(1:end-1)); %fft(uA_mean);
omegaA_f = fft(omegaA_mean(1:end-1)); %fft(omegaA_mean); 
% Motor B
UB_f =  fft(uB_mean(1:end-1)); %fft(uB_mean);
omegaB_f = fft(omegaB_mean(1:end-1)); %fft(omegaB_mean); 

%cut the data at the nyquist frequency
indices = 2:floor(numel(f)/2); % the first element corresponds to 0Hz that mean nothing in a log plot
f = f(indices);
% Motor A
UA_f = UA_f(indices);
omegaA_f = omegaA_f(indices);
% Motor B
UB_f = UB_f(indices);
omegaB_f = omegaB_f(indices);

% compute the empiric frequency response
% Motor A
A_FRF_emp = omegaA_f./UA_f;
A_mag_emp = 20*log10(abs(A_FRF_emp));
A_phs_emp = 180/pi*unwrap(angle(A_FRF_emp)); 
A_phs_emp = 360*ceil(-A_phs_emp(1)/360) + A_phs_emp; % just a trick to avoid not-meaningful 360deg differences at the start
% Motor B
B_FRF_emp = omegaB_f./UB_f;
B_mag_emp = 20*log10(abs(B_FRF_emp));
B_phs_emp = 180/pi*unwrap(angle(B_FRF_emp)); 
B_phs_emp = 360*ceil(-B_phs_emp(1)/360) + B_phs_emp; % just a trick to avoid not-meaningful 360deg differences at the start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if figures==1
figure(4)
% Motor A
hold on;
subplot(2,2,1);
semilogx(f, A_mag_emp, 'LineWidth', 1);
grid on;
xlabel('f [Hz]');
xlim([f(1) f(end)]);
ylabel('|FRF|_{A} [dB]');
subplot(2,2,3);
semilogx(f,A_phs_emp, 'LineWidth', 1);
grid on;
axis tight;
xlabel('f  [Hz]');
ylabel('\phi(FRF)_{A} [^\circ]');
xlim([f(1) f(end)]);
axis tight;
% Motor B
hold on;
subplot(2,2,2);
semilogx(f, B_mag_emp, 'LineWidth', 1);
grid on;
xlabel('f [Hz]');
xlim([f(1) f(end)]);
ylabel('|FRF|_{B} [dB]');
subplot(2,2,4);
semilogx(f,B_phs_emp, 'LineWidth', 1);
grid on;
axis tight;
xlabel('f  [Hz]');
ylabel('\phi(FRF)_{B} [^\circ]');
xlim([f(1) f(end)]);
axis tight;
end
%% Least-squares without data filtering
% assume the model to have the same shape of the discrete model we obtained
% on paper: H(s) = A/(z^2 + Bz)
%   k = 3:end
%   k-1 = 2:end-1
%   k-2 = 1:end-2

% collect the signals appearing in the difference equation
A_phi = [-omegaA(2:end-1), uA(1:end-2)];
B_phi = [-omegaB(2:end-1), uB(1:end-2)];

% perform the fit to get the desired parameters
A_theta_LS = A_phi\omegaA(3:end);
B_theta_LS = B_phi\omegaB(3:end);

% build the identified model
A_num_LS = [A_theta_LS(2)];
A_den_LS = [1, A_theta_LS(1), 0];
A_sys_LS = tf(A_num_LS, A_den_LS, Ts);

B_num_LS = [B_theta_LS(2)];
B_den_LS = [1, B_theta_LS(1), 0];
B_sys_LS = tf(B_num_LS, B_den_LS, Ts);

% natural frequecy, damping ratio and poles of the system
[A_wn_LS,A_zeta_LS, A_p_LS] = damp(A_sys_LS);  
[B_wn_LS,B_zeta_LS, B_p_LS] = damp(B_sys_LS); 

% compute the frequency response of the identified model
% Motor A
A_FRF_LS = squeeze(freqresp(A_sys_LS,2*pi*f));
A_mag_LS = 20*log10(abs(A_FRF_LS));
A_phs_LS = 180/pi*unwrap(angle(A_FRF_LS)); 
A_phs_LS = 360*ceil(-A_phs_LS(1)/360) + A_phs_LS;
% Motor B
B_FRF_LS = squeeze(freqresp(B_sys_LS,2*pi*f));
B_mag_LS = 20*log10(abs(B_FRF_LS));
B_phs_LS = 180/pi*unwrap(angle(B_FRF_LS)); 
B_phs_LS = 360*ceil(-B_phs_LS(1)/360) + B_phs_LS;

% time domain evaluation of the model 
A_omega_LS = lsim(A_sys_LS,uA,t);
A_error_LS = abs(omegaA - A_omega_LS);

B_omega_LS = lsim(B_sys_LS,uB,t);
B_error_LS = abs(omegaB - B_omega_LS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if figures==1
figure(5)
% Motor A
hold on;
subplot(2,2,1);
semilogx(f,A_mag_emp, f, A_mag_LS);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('|FRF|_{A}  [dB]');
legend('empirical', 'estimated','Location','SouthWest');
subplot(2,2,3);
semilogx(f, A_phs_emp, f, A_phs_LS);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('\phi(FRF)_{A}  [^\circ]');
legend('empirical', 'estimated','Location','SouthWest');
% Motor B
hold on;
subplot(2,2,2);
semilogx(f,B_mag_emp, f, B_mag_LS);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('|FRF|_{B}  [dB]');
legend('empirical', 'estimated','Location','SouthWest');
subplot(2,2,4);
semilogx(f, B_phs_emp, f, B_phs_LS);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('\phi(FRF)_{B}  [^\circ]');
legend('empirical', 'estimated','Location','SouthWest');

figure(6)
% Motor A
subplot(2,2,1);
plot(t,omegaA,t,A_omega_LS);
grid on;
legend('empirical','estimated','Location','SouthWest');
xlabel('time [s]');
ylabel('\omega_{A} [rad/s]');
axis tight;
subplot(2,2,3)
grid on;
plot(t,A_error_LS);
legend('error_{LS}');
xlabel('time [s]');
ylabel('error \omega_{A} [rad/s]');
axis tight;
% Motor B
subplot(2,2,2);
plot(t,omegaB,t,B_omega_LS);
grid on;
legend('empirical','estimated','Location','SouthWest');
xlabel('time [s]');
ylabel('\omega_{B} [rad/s]');
axis tight;
subplot(2,2,4)
grid on;
plot(t,B_error_LS);
legend('error_{LS}');
xlabel('time [s]');
ylabel('error \omega_{B} [rad/s]');
axis tight;

figure(7)
subplot(1,2,1);
pzmap(A_sys_LS);
title('Pole-Zero Map (Motor A)');
subplot(1,2,2);
pzmap(B_sys_LS);
title('Pole-Zero Map (Motor B)')
end
%% Least-squares with a low-pass filter applied to the input and output data
% cutOff_frequency
order = 6;    
magnitude = -3; % frequency where magnitude of system = -3dB
% Motor A
[A_d,A_ix] = min(abs(A_mag_LS-magnitude));
A_f_cutOff = f(A_ix);
[A_B_filt,A_A_filt] = butter(order, A_f_cutOff/(fs/2));
% Motor B
[B_d,B_ix] = min(abs(B_mag_LS-magnitude));
B_f_cutOff = f(B_ix);
[B_B_filt,B_A_filt] = butter(order, B_f_cutOff/(fs/2));

% apply the filter to both input and output
uA_filt = filter(A_B_filt, A_A_filt, uA);   
omegaA_filt = filter(A_B_filt, A_A_filt, omegaA);

uB_filt = filter(B_B_filt, B_A_filt, uB);   
omegaB_filt = filter(B_B_filt, B_A_filt, omegaB);

% repeat the identification
A_phi_filt = [-omegaA_filt(2:end-1), uA_filt(1:end-2)];
A_theta_filt = A_phi_filt\omegaA_filt(3:end);

B_phi_filt = [-omegaB_filt(2:end-1), uB_filt(1:end-2)];
B_theta_filt = B_phi_filt\omegaB_filt(3:end);

% build the identified model with filter
A_num_filt = [A_theta_filt(2)];
A_den_filt = [1, A_theta_filt(1), 0];
A_sys_filt = tf(A_num_filt, A_den_filt, Ts);

B_num_filt = [B_theta_filt(2)];
B_den_filt = [1, B_theta_filt(1), 0];
B_sys_filt = tf(B_num_filt, B_den_filt, Ts);

% natural frequecy, damping ratio and poles of the system
[A_wn_filt,A_zeta_filt, A_p_filt] = damp(A_sys_filt); 
[B_wn_filt,B_zeta_filt, B_p_filt] = damp(B_sys_filt); 

% compute the frequency response of the identified model
A_FRF_filt = squeeze(freqresp(A_sys_filt,2*pi*f));
A_mag_filt = 20*log10(abs(A_FRF_filt));
A_phs_filt = 180/pi*unwrap(angle(A_FRF_filt)); 
A_phs_filt = 360*ceil(-A_phs_filt(1)/360) + A_phs_filt;

B_FRF_filt = squeeze(freqresp(B_sys_filt,2*pi*f));
B_mag_filt = 20*log10(abs(B_FRF_filt));
B_phs_filt = 180/pi*unwrap(angle(B_FRF_filt)); 
B_phs_filt = 360*ceil(-B_phs_filt(1)/360) + B_phs_filt;

% time domain evaluation of the filterd model 
A_omega_filt = lsim(A_sys_filt,uA,t);
A_error_filt = abs(omegaA - A_omega_filt);

B_omega_filt = lsim(B_sys_filt,uB,t);
B_error_filt = abs(omegaB - B_omega_filt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if figures==1
figure(8)
% Motor A
hold on;
subplot(2,2,1);
semilogx(f,A_mag_emp, f, A_mag_LS, f, A_mag_filt);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('|FRF|_{A}  [dB]');
legend('empirical', 'estimated','low-pass','Location','SouthWest');
subplot(2,2,3);
semilogx(f, A_phs_emp, f, A_phs_LS, f, A_phs_filt);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('\phi(FRF)_{A}  [^\circ]');
legend('empirical', 'estimated','low-pass','Location','SouthWest');
% Motor B
hold on;
subplot(2,2,2);
semilogx(f,B_mag_emp, f, B_mag_LS, f, B_mag_filt);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('|FRF|_{B}  [dB]');
legend('empirical', 'estimated','low-pass','Location','SouthWest');
subplot(2,2,4);
semilogx(f, B_phs_emp, f, B_phs_LS, f, B_phs_filt);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('\phi(FRF)_{B}  [^\circ]');
legend('empirical', 'estimated','low-pass','Location','SouthWest');

figure(9)
% Motor A
subplot(2,2,1);
plot(t,omegaA,t,A_omega_LS,t,A_omega_filt);
grid on;
legend('empirical', 'estimated','low-pass','Location','SouthWest');
xlabel('time [s]');
ylabel('\omega_{A} [rad/s]');
axis tight;
subplot(2,2,3);
grid on;
plot(t,A_error_LS,t,A_error_filt);
legend('error_{LS}','error_{Filt}');
xlabel('time [s]');
ylabel('error \omega_{A} [rad/s]');
axis tight;
% Motor B
subplot(2,2,2);
plot(t,omegaB,t,B_omega_LS,t,B_omega_filt);
grid on;
legend('empirical', 'estimated','low-pass','Location','SouthWest');
xlabel('time [s]');
ylabel('\omega_{B} [rad/s]');
axis tight;
subplot(2,2,4);
grid on;
plot(t,B_error_LS,t,B_error_filt);
legend('error_{LS}','error_{Filt}');
xlabel('time [s]');
ylabel('error \omega_{B} [rad/s]');
axis tight;

figure(10)
subplot(1,2,1);
hold on;
pzmap(A_sys_LS,A_sys_filt);
title('Pole-Zero Map (Motor A)')
legend('estimated','low-pass','Location','SouthWest');
subplot(1,2,2);
hold on;
pzmap(B_sys_LS,B_sys_filt);
title('Pole-Zero Map (Motor B)')
legend('estimated','low-pass','Location','SouthWest');

end
%% Iterative identification and input and output data filtering 
difference = 10^-12;
% Motor A
A_theta_iter = A_theta_filt;
A_currentDen = A_theta_iter(1);
A_previousDen = 0;
A_denValues = [A_currentDen,A_previousDen];
A_iteration = 0;

while abs(A_denValues(1)-A_denValues(2))>difference
    A_iteration = A_iteration+1;
    A_denValues(2) = A_theta_iter(1);  %setting previous
    
    % filter signals using current denominator
    uA_iter = filter(1,[1,A_theta_iter(1),0],uA);
    omegaA_iter = filter(1,[1,A_theta_iter(1),0],omegaA);
    % repeat the identification
    A_phi_iter = [-omegaA_iter(2:end-1), uA_iter(1:end-2)];
    A_theta_iter = A_phi_iter\omegaA_iter(3:end);

    A_denValues(1) = A_theta_iter(1);  % setting current
    if A_iteration>100;
        break;
    end
end

% Motor B
B_theta_iter = B_theta_filt;
B_currentDen = B_theta_iter(1);
B_previousDen = 0;
B_denValues = [B_currentDen,B_previousDen];
B_iteration = 0;

while abs(B_denValues(1)-B_denValues(2))>difference
    B_iteration = B_iteration+1;
    B_denValues(2) = B_theta_iter(1);  %setting previous
    
    % filter signals using current denominator
    uB_iter = filter(1,[1,B_theta_iter(1),0],uB);
    omegaB_iter = filter(1,[1,B_theta_iter(1),0],omegaB);
    % repeat the identification
    B_phi_iter = [-omegaB_iter(2:end-1), uB_iter(1:end-2)];
    B_theta_iter = B_phi_iter\omegaB_iter(3:end);

    B_denValues(1) = B_theta_iter(1);  % setting current
    if B_iteration>100;
        break;
    end
end

A_num_iter = [A_theta_iter(2)];
A_den_iter = [1, A_theta_iter(1), 0];
A_sys_iter = tf(A_num_iter, A_den_iter, Ts);

B_num_iter = [B_theta_iter(2)];
B_den_iter = [1, B_theta_iter(1), 0];
B_sys_iter = tf(B_num_iter, B_den_iter, Ts);

% natural frequecy, damping ratio and poles of the system
[A_wn_iter,A_zeta_iter, A_p_iter] = damp(A_sys_iter) 
[B_wn_iter,B_zeta_iter, B_p_iter] = damp(B_sys_iter) 

% compute the frequency response of the identified model
A_FRF_iter = squeeze(freqresp(A_sys_iter,2*pi*f));
A_mag_iter = 20*log10(abs(A_FRF_iter));
A_phs_iter = 180/pi*unwrap(angle(A_FRF_filt)); 
A_phs_iter = 360*ceil(-A_phs_iter(1)/360) + A_phs_iter;

B_FRF_iter = squeeze(freqresp(B_sys_iter,2*pi*f));
B_mag_iter = 20*log10(abs(B_FRF_iter));
B_phs_iter = 180/pi*unwrap(angle(B_FRF_filt)); 
B_phs_iter = 360*ceil(-B_phs_iter(1)/360) + B_phs_iter;

% time domain evaluation of the filterd model 
A_omega_iter = lsim(A_sys_iter,uA,t);
A_error_iter = abs(omegaA - A_omega_iter);

B_omega_iter = lsim(B_sys_iter,uB,t);
B_error_iter = abs(omegaB - B_omega_iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if figures==1
figure(11)
% Motor A
hold on;
subplot(2,2,1);
semilogx(f,A_mag_emp, f, A_mag_LS, f, A_mag_filt, f, A_mag_iter);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('|FRF|_{A}  [dB]');
legend('empirical', 'estimated','low-pass','iterative','Location','SouthWest');
subplot(2,2,3);
semilogx(f, A_phs_emp, f, A_phs_LS, f, A_phs_filt, f, A_phs_iter);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('\phi(FRF)_{A}  [^\circ]');
legend('empirical', 'estimated','low-pass','iterative','Location','SouthWest');
% Motor B
hold on;
subplot(2,2,2);
semilogx(f,B_mag_emp, f, B_mag_LS, f, B_mag_filt, f, B_mag_iter);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('|FRF|_{B}  [dB]');
legend('empirical', 'estimated','low-pass','iterative','Location','SouthWest');
subplot(2,2,4);
semilogx(f, B_phs_emp, f, B_phs_LS, f, B_phs_filt, f, B_phs_iter);
grid on;
xlim([f(1) f(end)]);
xlabel('f  [Hz]');
ylabel('\phi(FRF)_{B}  [^\circ]');
legend('empirical', 'estimated','low-pass','iterative','Location','SouthWest');

figure(12)
% Motor A
subplot(2,2,1);
plot(t,omegaA,t,A_omega_LS,t,A_omega_filt,t,A_omega_iter);
grid on;
legend('empirical', 'estimated','low-pass','iterative','Location','SouthWest');
xlabel('time [s]');
ylabel('\omega_{A} [rad/s]');
axis tight;
subplot(2,2,3);
grid on;
plot(t,A_error_LS,t,A_error_filt,t,A_error_iter);
legend('error_{LS}','error_{Filt}','error_{Iter}');
xlabel('time [s]');
ylabel('error \omega_{A} [rad/s]');
axis tight;
% Motor B
subplot(2,2,2);
plot(t,omegaB,t,B_omega_LS,t,B_omega_filt,t,B_omega_iter);
grid on;
legend('empirical', 'estimated','low-pass','iterative','Location','SouthWest');
xlabel('time [s]');
ylabel('\omega_{B} [rad/s]');
axis tight;
subplot(2,2,4);
grid on;
plot(t,B_error_LS,t,B_error_filt,t,B_error_iter);
legend('error_{LS}','error_{Filt}','error_{Iter}');
xlabel('time [s]');
ylabel('error \omega_{B} [rad/s]');
axis tight;

figure(13)
hold on;
subplot(1,2,1);
pzmap(A_sys_LS,A_sys_filt,A_sys_iter);
title('Pole-Zero Map (Motor A)')
legend('estimated','low-pass','iterative','Location','SouthWest');
hold on;
subplot(1,2,2);
pzmap(B_sys_LS,B_sys_filt,B_sys_iter);
title('Pole-Zero Map (Motor B)')
legend('estimated','low-pass','iterative','Location','SouthWest');
end
%% Validate your identified model experimentally
% Loading the step input
% Define the relative path to the CSV file
relativePath_6V = 'step_input_6V.csv';

% Get the current working directory
currentDir = pwd;

% Construct the full file path by combining the current directory and the relative path
csvfile = fullfile(currentDir, relativePath_6V);
%csvfile = 'C:\Users\toonm\OneDrive\Documenten\KULeuven\2023-2024\Sem1\Systems and Control Theory\Control theory Assigment\Assignment 1\excitationSignal_pulseTrain.csv';
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels

indexFirst = find(data(:,3)~=0,1,"first")-50;% we started with zero state in Arduino
indexLast = find(data(:,3)~=0,1,"last")-247; % so x-axis gets bigger

uA_step = data(indexFirst:indexLast,3);
omegaA_step = data(indexFirst:indexLast,4);

uB_step = data(indexFirst:indexLast,6);
omegaB_step = data(indexFirst:indexLast,7);

N_step = length(uA_step);
t_step = [0:N_step-1]'*Ts;

% without filter
A_omega_step_LS = lsim(A_sys_LS,uA_step,t_step);
A_error_step_LS = (omegaA_step - A_omega_step_LS);

B_omega_step_LS = lsim(B_sys_LS,uB_step,t_step);
B_error_step_LS = (omegaB_step - B_omega_step_LS);

% with low-pass filter
A_omega_step_filt = lsim(A_sys_filt,uA_step,t_step);
A_error_step_filt = (omegaA_step - A_omega_step_filt);

B_omega_step_filt = lsim(B_sys_filt,uB_step,t_step);
B_error_step_filt = (omegaB_step - B_omega_step_filt);

% with iterative procedure
A_omega_step_iter = lsim(A_sys_iter,uA_step,t_step);
A_error_step_iter = (omegaA_step - A_omega_step_iter);

B_omega_step_iter = lsim(B_sys_iter,uB_step,t_step);
B_error_step_iter = (omegaB_step - B_omega_step_iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if figures==1
figure(14)
% total plot
% motor A
subplot(1,2,1);
plot(t_step,omegaA_step,t_step,A_omega_step_LS,t_step,A_omega_step_filt,t_step,A_omega_step_iter,'LineWidth',1);
grid on;
ylim([0 14]);
xlabel('Time [s]');
ylabel('\omega_{A} [rad/s]');
legend('empirical', 'estimated','low-pass','iterative','Location','Southeast');
% motor B
subplot(1,2,2);
plot(t_step,omegaB_step,t_step,B_omega_step_LS,t_step,B_omega_step_filt,t_step,B_omega_step_iter,'LineWidth',1);
grid on;
ylim([0 14]);
xlabel('Time [s]');
ylabel('\omega_{B} [rad/s]');
legend('empirical', 'estimated','low-pass','iterative','Location','Southeast');
% zoomed-in
figure(15)
subplot(1,2,1);
plot(t_step,omegaA_step,t_step,A_omega_step_LS,t_step,A_omega_step_filt,t_step,A_omega_step_iter,'LineWidth',1);
grid on;
ylim([12.2 13]);
xlabel('Time [s]');
ylabel('\omega_{A} [rad/s]');
legend('empirical', 'estimated','low-pass','iterative','Location','Southeast');
% motor B
subplot(1,2,2);
plot(t_step,omegaB_step,t_step,B_omega_step_LS,t_step,B_omega_step_filt,t_step,B_omega_step_iter,'LineWidth',1);
grid on;
ylim([12.2 12.75]);
xlabel('Time [s]');
ylabel('\omega_{B} [rad/s]');
legend('empirical', 'estimated','low-pass','iterative','Location','Southeast');

% error

figure(16)
% motor A
subplot(1,2,1);
plot(t_step,A_error_step_LS,t_step,A_error_step_filt,t_step,A_error_step_iter);
grid on;
ylim([-2 0.5]);
xlabel('Time [s]');
ylabel('error \omega_{A} [rad/s]');
legend('error_{LS}','error_{Filt}','error_{Iter}','Location','southeast');
% motor B
subplot(1,2,2);
plot(t_step,B_error_step_LS,t_step,B_error_step_filt,t_step,B_error_step_iter);
grid on;
ylim([-2 0.5]);
xlabel('Time [s]');
ylabel('error \omega_{B} [rad/s]');
legend('error_{LS}','error_{Filt}','error_{Iter}','Location','Southeast');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the cart while it is on the ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the excitation signal
% Define the relative path to the CSV file
relativePath_6V = 'excitationSignal_pulseTrain_6V_onGround.csv';

csvfile = fullfile(currentDir, relativePath_6V);
%csvfile = 'C:\Users\toonm\OneDrive\Documenten\KULeuven\2023-2024\Sem1\Systems and Control Theory\Control theory Assigment\Assignment 1\excitationSignal_pulseTrain.csv';
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels

indexFirst = find(data(:,3)~=0,1,"first")-50;% we started with zero state in Arduino
indexLast = find(data(:,3)~=0,1,"last")+50; % so it's clear we want from -6V back to 0V
uA_onGround = data(indexFirst:indexLast,3);
omegaA_onGround = data(indexFirst:indexLast,4);
thetaA_onGround = data(indexFirst:indexLast,5);

uB_onGround = data(indexFirst:indexLast,6);
omegaB_onGround = data(indexFirst:indexLast,7);
thetaB_onGround = data(indexFirst:indexLast,8);
if figures==1
figure();
grid on
hold on;
n1 = 50;
n2= 150;
axis([min(t(1:n2-n1+1)), max(t(1:n2-n1+1)), 0, 16]); % Set both x and y-axis limits
plot(t(1:n2-n1+1), omegaA_onGround(n1:n2), 'k', 'LineWidth', 2); % 'k' represents black color

plot(t(1:n2-n1+1), lsim(A_sys_iter, uA_onGround(n1:n2)), '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2); % [0.5, 0.5, 0.5] represents gray color
legend('Measured response with cart on the ground', 'Simulated response using identified model', 'FontSize', 12);

xlabel('Time (s)');
ylabel('\omega_A (rad/s)');
hold off;


figure();
hold on;
grid on
axis([min(t(1:n2-n1+1)), max(t(1:n2-n1+1)), 0, 4]); % Set both x and y-axis limits
diff=lsim(A_sys_iter, uA_onGround(n1:n2))-omegaA_onGround(n1:n2);
plot(t(1:n2-n1+1), diff, 'k', 'LineWidth', 2); % 'k' represents black color
xlabel('Time (s)');
ylabel('\omega_A (rad/s)');
hold off;
uB_onGround = data(indexFirst:indexLast,6);
omegaB_onGround = data(indexFirst:indexLast,7);
thetaB_onGround = data(indexFirst:indexLast,8);
end
if figures==1
figure(17);

subplot(2, 1, 1);
hold on;
n1 = 1;
n2= 150;
plot(t(n1:n2), omegaA_onGround(n1:n2), 'LineWidth', 2);
plot(t(n1:n2), omegaA(n1:n2), 'LineWidth', 2);
plot(t(n1:n2), uA_onGround(n1:n2), 'LineWidth', 2);
grid on;
axis tight;
title('First Subplot');

subplot(2, 1, 2);
hold on;
plot(t(n1:n2), omegaA_onGround(n1:n2)-omegaA(n1:n2), 'LineWidth', 2);

grid on;
axis tight;
title('Second Subplot');
end
% Adjust overall figure title, labels, etc. if needed
A_phi_onGround = [-omegaA_onGround(2:end-1), uA_onGround(1:end-2)];
B_phi_onGround = [-omegaB_onGround(2:end-1), uB_onGround(1:end-2)];
% perform the fit to get the desired parameters
A_theta_LS_onGround = A_phi_onGround\omegaA_onGround(3:end);
B_theta_LS_onGround = B_phi_onGround\omegaB_onGround(3:end);
% build the identified model
A_num_LS_onGround = [A_theta_LS_onGround(2)];
A_den_LS_onGround = [1, A_theta_LS_onGround(1), 0];
A_sys_LS_onGround = tf(A_num_LS_onGround, A_den_LS_onGround, Ts);

B_num_LS_onGround = [B_theta_LS_onGround(2)];
B_den_LS_onGround = [1, B_theta_LS_onGround(1), 0];
B_sys_LS_onGround = tf(B_num_LS_onGround, B_den_LS_onGround, Ts);

% natural frequecy, damping ratio and poles of the system
[A_wn_LS_onGround,A_zeta_LS_onGround, A_p_LS_onGround] = damp(A_sys_LS_onGround);  
[B_wn_LS_onGround,B_zeta_LS_onGround, B_p_LS_onGround] = damp(B_sys_LS_onGround); 


% Print the text
disp(['Omega A {In air, without filtering}: ']);
A_sys_LS
disp(['Omega A {On ground, without filtering}: ']);
A_sys_LS_onGround
A_FRF_LS_onGround = squeeze(freqresp(A_sys_LS_onGround,2*pi*f));
A_mag_LS_onGround = 20*log10(abs(A_FRF_LS_onGround));

% Motor B
B_FRF_LS_onGround = squeeze(freqresp(B_sys_LS_onGround,2*pi*f));
B_mag_LS_onGround = 20*log10(abs(B_FRF_LS_onGround));
order = 6;    
magnitude = -3; % frequency where magnitude of system = -3dB
% Motor A
[A_d_onGround,A_ix_onGround] = min(abs(A_mag_LS_onGround-magnitude));
A_f_cutOff_onGround = f(A_ix_onGround);
[A_B_filt_onGround,A_A_filt_onGround] = butter(order, A_f_cutOff_onGround/(fs/2));
% Motor B
[B_d_onGround,B_ix_onGround] = min(abs(B_mag_LS_onGround-magnitude));
B_f_cutOff_onGround = f(B_ix_onGround);
[B_B_filt_onGround,B_A_filt_onGround] = butter(order, B_f_cutOff_onGround/(fs/2));

% apply the filter to both input and output
uA_filt_onGround = filter(A_B_filt_onGround, A_A_filt_onGround, uA_onGround);   
omegaA_filt_onGround = filter(A_B_filt_onGround, A_A_filt_onGround, omegaA_onGround);

uB_filt_onGround= filter(B_B_filt_onGround, B_A_filt_onGround, uB_onGround);   
omegaB_filt_onGround = filter(B_B_filt_onGround, B_A_filt_onGround, omegaB_onGround);

% repeat the identification
A_phi_filt_onGround = [-omegaA_filt_onGround(2:end-1), uA_filt_onGround(1:end-2)];
A_theta_filt_onGround = A_phi_filt_onGround\omegaA_filt_onGround(3:end);

B_phi_filt_onGround = [-omegaB_filt_onGround(2:end-1), uB_filt_onGround(1:end-2)];
B_theta_filt_onGround = B_phi_filt_onGround\omegaB_filt_onGround(3:end);

% build the identified model with filter
A_num_filt_onGround = [A_theta_filt_onGround(2)];
A_den_filt_onGround = [1, A_theta_filt_onGround(1), 0];
A_sys_filt_onGround = tf(A_num_filt_onGround, A_den_filt_onGround, Ts);

B_num_filt_onGround = [B_theta_filt_onGround(2)];
B_den_filt_onGround = [1, B_theta_filt_onGround(1), 0];
B_sys_filt_onGround = tf(B_num_filt_onGround, B_den_filt_onGround, Ts);

% Print the text
disp(['Omega A {In air, with filtering}: ']);
A_sys_filt
disp(['Omega A {On ground, withfiltering}: ']);
A_sys_filt_onGround

A_theta_iter_onGround = A_theta_filt_onGround;
A_currentDen_onGround = A_theta_iter_onGround(1);
A_previousDen_onGround = 0;
A_denValues_onGround = [A_currentDen_onGround,A_previousDen_onGround];
A_iteration = 0;

while abs(A_denValues_onGround(1)-A_denValues_onGround(2))>difference
    A_iteration = A_iteration+1;
    A_denValues_onGround(2) = A_theta_iter_onGround(1);  %setting previous
    
    % filter signals using current denominator
    uA_iter_onGround = filter(1,[1,A_theta_iter_onGround(1),0],uA_onGround);
    omegaA_iter_onGround = filter(1,[1,A_theta_iter_onGround(1),0],omegaA_onGround);
    % repeat the identification
    A_phi_iter_onGround = [-omegaA_iter_onGround(2:end-1), uA_iter_onGround(1:end-2)];
    A_theta_iter_onGround = A_phi_iter_onGround\omegaA_iter_onGround(3:end);

    A_denValues_onGround(1) = A_theta_iter_onGround(1);  % setting current
    if A_iteration>100
        break;
    end
end

% Motor B
B_theta_iter_onGround = B_theta_filt_onGround;
B_currentDen_onGround = B_theta_iter_onGround(1);
B_previousDen_onGround = 0;
B_denValues_onGround = [B_currentDen_onGround,B_previousDen_onGround];
B_iteration = 0;

while abs(B_denValues_onGround(1)-B_denValues_onGround(2))>difference
    B_iteration = B_iteration+1;
    B_denValues_onGround(2) = B_theta_iter_onGround(1);  %setting previous
    
    % filter signals using current denominator
    uB_iter_onGround = filter(1,[1,B_theta_iter_onGround(1),0],uB_onGround);
    omegaB_iter_onGround = filter(1,[1,B_theta_iter_onGround(1),0],omegaB_onGround);
    % repeat the identification
    B_phi_iter_onGround = [-omegaB_iter_onGround(2:end-1), uB_iter_onGround(1:end-2)];
    B_theta_iter_onGround = B_phi_iter_onGround\omegaB_iter_onGround(3:end);

    B_denValues_onGround(1) = B_theta_iter_onGround(1);  % setting current
    if B_iteration>100;
        break;
    end
end

A_num_iter_onGround = [A_theta_iter_onGround(2)];
A_den_iter_onGround = [1, A_theta_iter_onGround(1), 0];
A_sys_iter_onGround = tf(A_num_iter_onGround, A_den_iter_onGround, Ts);

B_num_iter_onGround = [B_theta_iter_onGround(2)];
B_den_iter_onGround = [1, B_theta_iter_onGround(1), 0];
B_sys_iter_onGround = tf(B_num_iter_onGround, B_den_iter_onGround, Ts);
% Print the text
disp(['Omega A {In air, with iterative filtering}: ']);
A_sys_iter
disp(['Omega A {On ground, with ilterative filtering}: ']);
A_sys_iter_onGround
N_step = length(uA_onGround);
t_step = [0:N_step-1]'*Ts;
A_omega_step_iter_onGround = lsim(A_sys_iter_onGround,uA_onGround,t_step);
if figures==1
figure(24)
hold on
plot(t,A_omega_step_iter_onGround)
plot(t,uA_onGround)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the superposition principle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the excitation signal
% Define the relative path to the CSV file
relativePath_4V = 'step_input_4V.csv';
relativePath_6V = 'step_input_6V.csv';
relativePath_10V = 'step_input_10V.csv';
% Get the current working directory
currentDir = pwd;

% Construct the full file path by combining the current directory and the relative path
csvfile = fullfile(currentDir, relativePath_4V);
%csvfile = 'C:\Users\toonm\OneDrive\Documenten\KULeuven\2023-2024\Sem1\Systems and Control Theory\Control theory Assigment\Assignment 1\excitationSignal_pulseTrain.csv';
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels

indexFirst = find(data(:,3)~=0,1,"first")-50;% we started with zero state in Arduino
indexLast = find(data(:,3)~=0,1,"last"); % so it's clear we want from -6V back to 0V
uA_4V = data(indexFirst:indexLast,3);
omegaA_4V = data(indexFirst:indexLast,4);
thetaA_4V = data(indexFirst:indexLast,5);

csvfile = fullfile(currentDir, relativePath_6V);

labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels

indexFirst = find(data(:,3)~=0,1,"first")-50;
indexLast = find(data(:,3)~=0,1,"last");
uA_6V = data(indexFirst:indexLast,3);
omegaA_6V = data(indexFirst:indexLast,4);
thetaA_6V = data(indexFirst:indexLast,5);

csvfile = fullfile(currentDir, relativePath_10V);
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels

indexFirst = find(data(:,3)~=0,1,"first")-50;% we started with zero state in Arduino
indexLast = find(data(:,3)~=0,1,"last"); % so it's clear we want from -6V back to 0V
uA_10V = data(indexFirst:indexLast,3);
omegaA_10V = data(indexFirst:indexLast,4);
thetaA_10V = data(indexFirst:indexLast,5);
if figures==1
figure(19);
hold on;
N_step = length(uA_step);
t_step = [0:100-40]*Ts;
plot(t_step ,omegaA_10V(40:100), 'LineWidth', 2, 'Color', 'k', 'DisplayName', '\omega_A(10V)'); % Black line
plot(t_step ,omegaA_6V(40:100) + omegaA_4V(40:100), 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5], 'DisplayName', '\omega_A(4V) + \omega_A(6V)'); % Gray line
hold off; 
xlabel('Time (s)');
ylabel('\omega_A (rad/s)');
legend('show');
grid on;
axis tight;
end

