clc;
close all
clear all
currentDir = pwd;
currentDir_Assignment3_3E_Data = fullfile(currentDir,'Assignment3_3E_data');

endtime = 1.4;    
time = linspace(0,endtime,endtime/0.01+1);
N = length(time);
fontsize_Value = 12;

%%%% 3E_Q/R = 0.01 %%%%
csvfile_3E_QR_001 = fullfile(currentDir_Assignment3_3E_Data,'Assignment3_3E_Q_8_R_6.csv');
measurements_3E_QR_001 = dlmread(csvfile_3E_QR_001, ',', 2, 0); % Data follows the labels

indexFirst_QR_001 = find(measurements_3E_QR_001(:,2)~=0,1,"first")-1;
indexLast_QR_001 = indexFirst_QR_001+N-1;

measuredDistance_QR_001 = -1*(measurements_3E_QR_001(indexFirst_QR_001:indexLast_QR_001,3));
estimatedDistance_QR_001 = measurements_3E_QR_001(indexFirst_QR_001:indexLast_QR_001,12);

%%%% 3E_Q/R = 0.1 %%%%
csvfile_3E_QR_01 = fullfile(currentDir_Assignment3_3E_Data,'Assignment3_3E_Q_7_R_6.csv');
measurements_3E_QR_01 = dlmread(csvfile_3E_QR_01, ',', 2, 0); % Data follows the labels

indexFirst_QR_01 = find(measurements_3E_QR_01(:,2)~=0,1,"first")-1;
indexLast_QR_01 = indexFirst_QR_01+N-1;

measuredDistance_QR_01 = -1*(measurements_3E_QR_01(indexFirst_QR_01:indexLast_QR_01,3));
estimatedDistance_QR_01 = measurements_3E_QR_01(indexFirst_QR_01:indexLast_QR_01,12);

%%%% 3E_Q/R = 1 %%%%
csvfile_3E_QR_1 = fullfile(currentDir_Assignment3_3E_Data,'Assignment3_3E_Q_6_R_6.csv');
measurements_3E_QR_1 = dlmread(csvfile_3E_QR_1, ',', 2, 0); % Data follows the labels

indexFirst_QR_1 = find(measurements_3E_QR_1(:,2)~=0,1,"first")-1;
indexLast_QR_1 = indexFirst_QR_1+N-1;

measuredDistance_QR_1 = -1*(measurements_3E_QR_1(indexFirst_QR_1:indexLast_QR_1,3));
estimatedDistance_QR_1 = measurements_3E_QR_1(indexFirst_QR_1:indexLast_QR_1,12);

%%%% 3E_Q/R = 10 %%%%
csvfile_3E_QR_10 = fullfile(currentDir_Assignment3_3E_Data,'Assignment3_3E_Q_5_R_6.csv');
measurements_3E_QR_10 = dlmread(csvfile_3E_QR_10, ',', 2, 0); % Data follows the labels

indexFirst_QR_10 = find(measurements_3E_QR_10(:,2)~=0,1,"first")-1;
indexLast_QR_10 = indexFirst_QR_10+N-1;

measuredDistance_QR_10 = -1*(measurements_3E_QR_10(indexFirst_QR_10:indexLast_QR_10,3));
estimatedDistance_QR_10 = measurements_3E_QR_10(indexFirst_QR_10:indexLast_QR_10,12);

%%%% 3E_Q/R = 100 %%%%
csvfile_3E_QR_100 = fullfile(currentDir_Assignment3_3E_Data,'Assignment3_3E_Q_4_R_6.csv');
measurements_3E_QR_100 = dlmread(csvfile_3E_QR_100, ',', 2, 0); % Data follows the labels

indexFirst_QR_100 = find(measurements_3E_QR_100(:,2)~=0,1,"first")-1;
indexLast_QR_100 = indexFirst_QR_100+N-1;

measuredDistance_QR_100 = -1*(measurements_3E_QR_100(indexFirst_QR_100:indexLast_QR_100,3));
estimatedDistance_QR_100 = measurements_3E_QR_100(indexFirst_QR_100:indexLast_QR_100,12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
subplot(3,1,1)
hold on;
plot(time,measuredDistance_QR_001,'LineWidth',2,'Color','k')
plot(time,estimatedDistance_QR_001,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('$\frac{Q}{R}=0.01$','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;
subplot(3,1,2)
hold on;
plot(time,measuredDistance_QR_01,'LineWidth',2,'Color','k')
plot(time,estimatedDistance_QR_01,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('$\frac{Q}{R}=0.1$','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;
subplot(3,1,3)
hold on;
plot(time,measuredDistance_QR_1,'LineWidth',2,'Color','k')
plot(time,estimatedDistance_QR_1,'LineWidth',2, 'Color','[0.5,0.5,0.5]','LineStyle','--')
title('$\frac{Q}{R}=1$','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',fontsize_Value)
grid on;
hold off;

figure(2)
figureSize = [100, 100, 1200, 600]; % [left bottom width height]
set(gcf, 'Position', figureSize);
subplot(2,1,1)
hold on;
plot(time,measuredDistance_QR_10,'LineWidth',2,'Color','k')
plot(time,estimatedDistance_QR_10,'LineWidth',2, 'Color','[0.3,0.3,0.3]','LineStyle','--')
title('$\frac{Q}{R}=10$','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]','FontSize',fontsize_Value)
grid on;
hold off;
subplot(2,1,2)
hold on;
plot(time,measuredDistance_QR_100,'LineWidth',2,'Color','k')
plot(time,estimatedDistance_QR_100,'LineWidth',2, 'Color','[0.3,0.3,0.3]','LineStyle','--')
title('$\frac{Q}{R}=100$','Interpreter', 'latex','FontSize',fontsize_Value)
legend('Distance measurement','Position estimate','FontSize',fontsize_Value,'Location','Southeast');
xlabel('Time [s]','FontSize',fontsize_Value)
ylabel('Position [m]','FontSize',fontsize_Value)
grid on;
hold off;