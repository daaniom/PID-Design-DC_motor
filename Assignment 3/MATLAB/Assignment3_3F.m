clc;
close all
clear all
currentDir = pwd;
currentDir_Assignment3_3F_Data = fullfile(currentDir,'Assignment3_3F_data');

%% Calculation closed-loop pole estimator for different speeds
A = 1;
B = 0.01;
C = -1;
D = 0;
speeds = [0.1 1 2 3 4 5 6];
K = 2;  % this value is chosen
Ts = 0.01;
P_DT= 1-K*Ts;
P_CT = - log(P_DT)/Ts;
P_CT_differentSpeed = P_CT*speeds;
P_DT_differentSpeed = exp(-P_CT_differentSpeed*Ts);

L_values = P_DT_differentSpeed - 1

L_control = zeros(size(P_DT_differentSpeed)); 

for i = 1:length(P_DT_differentSpeed)
    p = P_DT_differentSpeed(i);
    result = place(A', C', p);
    L_control(i) = result;
end
L_control

%% results
fontsize_Value = 12;

%%%% 10 times slower %%%%%
csvfile_3F_01 = fullfile(currentDir_Assignment3_3F_Data,'Assignment3_3F_01.csv');
measurements_3F_01 = dlmread(csvfile_3F_01, ',', 2, 0); % Data follows the labels

indexFirst_01 = find(measurements_3F_01(:,2)~=0,1,"first")-1;
%indexLast_01 = indexFirst_01+N-1;

measuredDistance_01 = -1*(measurements_3F_01(indexFirst_01:end,3));
estimatedDistance_01 = measurements_3F_01(indexFirst_01:end,12);

time_01 = linspace(0,(length(measuredDistance_01)-1)*Ts,length(measuredDistance_01));

%%%% 1 times faster %%%%
csvfile_3F_1 = fullfile(currentDir_Assignment3_3F_Data,'Assignment3_3F_1.csv');
measurements_3F_1 = dlmread(csvfile_3F_1, ',', 2, 0); % Data follows the labels

indexFirst_1 = find(measurements_3F_1(:,2)~=0,1,"first")-1;
%indexLast_1 = indexFirst_1+N-1;

measuredDistance_1 = -1*(measurements_3F_1(indexFirst_1:end,3));
estimatedDistance_1 = measurements_3F_1(indexFirst_1:end,12);

time_1 = linspace(0,(length(measuredDistance_1)-1)*Ts,length(measuredDistance_1));

%%%% 2 times faster %%%%
csvfile_3F_2 = fullfile(currentDir_Assignment3_3F_Data,'Assignment3_3F_2.csv');
measurements_3F_2 = dlmread(csvfile_3F_2, ',', 2, 0); % Data follows the labels

indexFirst_2 = find(measurements_3F_2(:,2)~=0,1,"first")-1;
%indexLast_1 = indexFirst_1+N-1;

measuredDistance_2 = -1*(measurements_3F_2(indexFirst_2:end,3));
estimatedDistance_2 = measurements_3F_2(indexFirst_2:end,12);

time_2 = linspace(0,(length(measuredDistance_2)-1)*Ts,length(measuredDistance_2));

%%%% 3 times faster %%%%
csvfile_3F_3 = fullfile(currentDir_Assignment3_3F_Data,'Assignment3_3F_3.csv');
measurements_3F_3 = dlmread(csvfile_3F_3, ',', 2, 0); % Data follows the labels

indexFirst_3 = find(measurements_3F_3(:,2)~=0,1,"first")-1;
%indexLast_1 = indexFirst_1+N-1;

measuredDistance_3 = -1*(measurements_3F_3(indexFirst_3:end,3));
estimatedDistance_3 = measurements_3F_3(indexFirst_3:end,12);

time_3 = linspace(0,(length(measuredDistance_3)-1)*Ts,length(measuredDistance_3));

%%%% 4 times faster %%%%
csvfile_3F_4 = fullfile(currentDir_Assignment3_3F_Data,'Assignment3_3F_4.csv');
measurements_3F_4 = dlmread(csvfile_3F_4, ',', 2, 0); % Data follows the labels

indexFirst_4 = find(measurements_3F_4(:,2)~=0,1,"first")-1;
%indexLast_1 = indexFirst_1+N-1;

measuredDistance_4 = -1*(measurements_3F_4(indexFirst_4:end,3));
estimatedDistance_4 = measurements_3F_4(indexFirst_4:end,12);

time_4 = linspace(0,(length(measuredDistance_4)-1)*Ts,length(measuredDistance_4));

%%%% 5 times faster %%%%
csvfile_3F_5 = fullfile(currentDir_Assignment3_3F_Data,'Assignment3_3F_5.csv');
measurements_3F_5 = dlmread(csvfile_3F_5, ',', 2, 0); % Data follows the labels

indexFirst_5 = find(measurements_3F_5(:,2)~=0,1,"first")-1;
%indexLast_1 = indexFirst_1+N-1;

measuredDistance_5 = -1*(measurements_3F_5(indexFirst_5:end,3));
estimatedDistance_5 = measurements_3F_5(indexFirst_5:end,12);

time_5 = linspace(0,(length(measuredDistance_5)-1)*Ts,length(measuredDistance_5));

%%%% 6 times faster %%%%
csvfile_3F_6 = fullfile(currentDir_Assignment3_3F_Data,'Assignment3_3F_6.csv');
measurements_3F_6 = dlmread(csvfile_3F_6, ',', 2, 0); % Data follows the labels

indexFirst_6 = find(measurements_3F_6(:,2)~=0,1,"first")-1;
%indexLast_1 = indexFirst_1+N-1;

measuredDistance_6 = -1*(measurements_3F_6(indexFirst_6:end,3));
estimatedDistance_6 = measurements_3F_6(indexFirst_6:end,12);

time_6 = linspace(0,(length(measuredDistance_6)-1)*Ts,length(measuredDistance_6));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
hold on;
plot(time_01,measuredDistance_01,'LineWidth',2,'Color','k')
plot(time_01,estimatedDistance_01,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('10 times slower','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;

figure(2)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
hold on;
plot(time_1,measuredDistance_1,'LineWidth',2,'Color','k')
plot(time_1,estimatedDistance_1,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('Same gain','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;

figure(3)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
hold on;
plot(time_2,measuredDistance_2,'LineWidth',2,'Color','k')
plot(time_2,estimatedDistance_2,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('2 times faster','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;

figure(4)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
hold on;
plot(time_3,measuredDistance_3,'LineWidth',2,'Color','k')
plot(time_3,estimatedDistance_3,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('3 times faster','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;

figure(5)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
hold on;
plot(time_4,measuredDistance_4,'LineWidth',2,'Color','k')
plot(time_4,estimatedDistance_4,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('4 times faster','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;

figure(6)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
hold on;
plot(time_5,measuredDistance_5,'LineWidth',2,'Color','k')
plot(time_5,estimatedDistance_5,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('5 times faster','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;

figure(7)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
hold on;
plot(time_6,measuredDistance_6,'LineWidth',2,'Color','k')
plot(time_6,estimatedDistance_6,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('6 times faster','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;


%%
endtime = 2.5;    
time = linspace(0,endtime,endtime/0.01+1);
N = length(time);

time_2 = time_2(1:N);
measuredDistance_2 =measuredDistance_2(1:N);
estimatedDistance_2 = estimatedDistance_2 (1:N);

time_6 = time_6(1:N);
measuredDistance_6 =measuredDistance_6(1:N);
estimatedDistance_6 = estimatedDistance_6 (1:N);

figure(8)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
subplot(2,1,1)
hold on;
plot(time_2,measuredDistance_2,'LineWidth',2,'Color','k')
plot(time_2,estimatedDistance_2,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('2 times faster','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;
subplot(2,1,2)
hold on;
plot(time_6,measuredDistance_6,'LineWidth',2,'Color','k')
plot(time_6,estimatedDistance_6,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('6 times faster','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;