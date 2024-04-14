
A_sys = tf([0.5251], [1, - 0.7338,0], 0.01);
B_sys = tf([0.473], [1, - 0.7594,0], 0.01);
Ts = 0.01;
currentDir= pwd();
% Rename the system , shorter name makes it easier
%A_sys = A_sys_iter_onGround
%B_sys = B_sys_iter_onGround


%% Choose a design value for the phase margin
% The PI-compensator is designed according to the frequency response 
% of the system following the method described in the cours slides
Pm_design=90;
w_eval = logspace(1,log10((1/Ts)*3.14),500);

[mag_A,phase_A,wout_A] = bode(A_sys,w_eval);
[Gm_A,Pm_A,Wcg_A,Wcp_A] = margin(A_sys);

[mag_B,phase_B,wout_B] = bode(B_sys,w_eval);
[Gm_B,Pm_B,Wcg_B,Wcp_B] = margin(B_sys);

% Calculate the Phase margin of the current open loop system for A and B
Pm_A;
Pm_B;
%Choose the expected value for the lag angle of the PID
lag_Phase_PI=15;
%% Calculation of the new crossover frequency
targetPhase = -180+Pm_design+lag_Phase_PI;
phase_A = reshape(phase_A,[],1);
phase_B = reshape(phase_B,[],1);
w_crossNew_A = interp1(phase_A,w_eval,targetPhase)
w_crossNew_B = interp1(phase_B,w_eval,targetPhase)
Ti_A = tand(90-lag_Phase_PI)/w_crossNew_A
Ti_B = tand(90-lag_Phase_PI)/w_crossNew_B

% Transforming PI compensator from the form (K/s)(s+1/Ti) to
% (kp*s+ki)/s, multiplication with K is done afterwards since it is still
% to be determined
kp_A = 1;
ki_A = kp_A/Ti_A;
PI_A = c2d(tf([kp_A ki_A],[1 0]),Ts);

kp_B = 1;
ki_B = kp_B/Ti_B;
PI_B = c2d(tf([kp_B ki_B],[1 0]),Ts);
% Calculate the forward loop with a K=1
FL_A = PI_A*A_sys;
FL_B = PI_B*B_sys;


% Calculate K that is needed to have the unity gain at the new crossover
% frequency
[mag_FL_A,phase_FL_A,wout_FL_A] = bode(FL_A,w_eval);
mag_FL_A = reshape(mag_FL_A,[],1);
kp_A = 1/interp1(w_eval,mag_FL_A,w_crossNew_A);
ki_A = kp_A/Ti_A;
PI_A = c2d(tf([kp_A ki_A],[1 0]),Ts);


FL_B = PI_B*B_sys;
[mag_FL_B,phase_FL_B,wout_FL_B] = bode(FL_B,w_eval);
mag_FL_B = reshape(mag_FL_B,[],1);
kp_B = 1/interp1(w_eval,mag_FL_B,w_crossNew_B);
ki_B = kp_B/Ti_B;
PI_B = c2d(tf([kp_B ki_B],[1 0]),Ts);
%Forward loop with the right gain

FL_A = A_sys*PI_A
figure(500);
margin(FL_A);
FL_B = B_sys*PI_B;
%figure(6546512)
%margin(FL_A)
CL_A = feedback(FL_A,1)
%figure(5000);
%bode(CL_A);

[mag, phase, w] = bode(CL_A);
w_3dB = w(find(abs(mag - 0.707) == min(abs(mag - 0.707))));

fprintf('Frequency at which -3dB is crossed: %.2f rad/s\n', w_3dB);
CL_B = feedback(FL_B,1)

%% Plotting the figures for case without any added weight
relativePath ='6rad_unweighted.csv'; 
csvfile = fullfile(currentDir, relativePath);
%csvfile = 'C:\Users\toonm\OneDrive\Documenten\KULeuven\2023-2024\Sem1\Systems and Control Theory\Control theory Assigment\Assignment 1\60dgFM.csv';
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels

indexFirst = find(data(:,3)~=0,1,"first")-10;% we started with zero state in Arduino
indexLast = indexFirst+90; 
omegaref_noWeight = data(indexFirst:indexLast,8);
omegaA_noWeight = data(indexFirst:indexLast,6);
errorA_noWeight= data(indexFirst:indexLast,4);
uA_noWeight= data(indexFirst:indexLast,2);

relativePath ='6rad_weighted.csv'; 
csvfile = fullfile(currentDir, relativePath);
%csvfile = 'C:\Users\toonm\OneDrive\Documenten\KULeuven\2023-2024\Sem1\Systems and Control Theory\Control theory Assigment\Assignment 1\60dgFM.csv';
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels

indexFirst = find(data(:,3)~=0,1,"first")-10;% we started with zero state in Arduino
indexLast = indexFirst+90; 
omegaref_addedWeight = data(indexFirst:indexLast,8);
omegaA_addedWeight = data(indexFirst:indexLast,6);
errorA_addedWeight= data(indexFirst:indexLast,4);
uA_addedWeight= data(indexFirst:indexLast,2);


N_step = length(uA_noWeight);
length(uA_noWeight)
t_step = [0:N_step-1]'*Ts;
color = [0.5, 0, 0];

% Plot of response to a step input
figure(5417);

plot(t_step, omegaref_noWeight, '--', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Reference step signal');
hold on;
plot(t_step, omegaA_noWeight, '-','Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', 'Measured \omega_{A}');
plot(t_step, lsim(CL_A, omegaref_noWeight, t_step), '-', 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Simulated \omega_{A}');
hold off;
xlabel('Time (s)');
ylabel('\omega (rad/s)');
legend('show','Location', 'SouthEast');
grid on;
saveas(gcf, 'plots/noDisturbance_stepresponse.png');

% Plot of the tracking error to a step reference
figure(5411);
hold on;
plot(t_step, omegaref_noWeight-omegaA_noWeight, '-','Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', 'Measured \omega_{A}');
plot(t_step, omegaref_noWeight-lsim(CL_A, omegaref_noWeight, t_step), '-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Simulated \omega_{A}');
hold off;
xlabel('Time (s)');
ylabel('\omega (rad/s)');
legend('show','Location', 'NorthEast');
grid on;
saveas(gcf, 'plots/noDisturbance_trackingerror.png');

% Plot of the control signal simulated and measured
figure(5421);
sim_trackerror = omegaref_noWeight - lsim(CL_A, omegaref_noWeight, t_step);

plot(t_step, uA_noWeight, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', 'Measured control singal');
hold on;
plot(t_step, lsim(PI_A, sim_trackerror, t_step), '-','Color', 'k', 'LineWidth', 2, 'DisplayName', 'Simulated control signal');

hold off;
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('show','Location', 'NorthEast');
grid on;
saveas(gcf, 'plots/noDisturbance_controlsignal.png');
figure(55);

plot(t_step, omegaref_noWeight, '--', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Reference step signal');
hold on;

plot(t_step, omegaA_noWeight, '-','Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', 'Measured \omega_{A} without disturbance');
plot(t_step, lsim(CL_A, omegaref_noWeight, t_step), '-', 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Simulated \omega_{A}');
plot(t_step, omegaA_addedWeight, '-','Color', color, 'LineWidth', 2, 'DisplayName', 'Measured \omega_{A} with disturbance');
hold off;
xlabel('Time (s)');
ylabel('\omega (rad/s)');
legend('show','Location', 'SouthEast');
grid on;
saveas(gcf, 'plots/Disturbance_stepresponse.png');

% Plot of the tracking error to a step reference
figure(5455555);
plot(t_step, omegaref_noWeight-omegaA_noWeight, '-','Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', 'Measured \omega_{A} without disturbance');
hold on;
plot(t_step, omegaref_noWeight-omegaA_addedWeight, '-','Color', color, 'LineWidth', 2, 'DisplayName', 'Measured \omega_{A} with disturbance');
plot(t_step, omegaref_noWeight-lsim(CL_A, omegaref_noWeight, t_step), '-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Simulated \omega_{A}');
hold off;
xlabel('Time (s)');
ylabel('\omega (rad/s)');
legend('show','Location', 'NorthEast');
grid on;
saveas(gcf, 'plots/Disturbance_trackingerror.png');

% Plot of the control signal simulated and measured
figure(54215555);
sim_trackerror = omegaref_noWeight - lsim(CL_A, omegaref_noWeight, t_step);

plot(t_step, uA_noWeight, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', 'Measured control singal');
hold on;
plot(t_step, uA_addedWeight, '-', 'Color', color, 'LineWidth', 2, 'DisplayName', 'Measured control singal with disturbance');
plot(t_step, lsim(PI_A, sim_trackerror, t_step), '-','Color', 'k', 'LineWidth', 2, 'DisplayName', 'Simulated control signal');

hold off;
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('show','Location', 'NorthEast');
grid on;
saveas(gcf, 'plots/Disturbance_controlsignal.png');


