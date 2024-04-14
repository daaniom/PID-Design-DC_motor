clc;
close all
clear all
currentDir = pwd;
currentDir_Assignment3_3C_Data = fullfile(currentDir,'Assignment3_3C_data');

assignment3_3C = true;
assignment3_3D = false;

endtime = 3.5;    % this time becomes smaller because for higher values of K the time measured becomes smaller due to the instabillity of these values
time = linspace(0,endtime,endtime/0.01+1);
N = length(time);
differentStart = 0;
if assignment3_3C
    differentStart = 50;
end
if assignment3_3D
    differentStart = 0
end

%% Assignment 3 3C:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   R = constant (=8*10^-6) & Q = changes %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 3C_Q_4_R_6 %%%%
% Q = 8*10-4 & R = 8*10-6
csvfile_3C_Q_4_R_6 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_4_R_6.csv');
measurements_3C_Q_4_R_6 = dlmread(csvfile_3C_Q_4_R_6, ',', 2, 0); % Data follows the labels

indexFirst_Q4R6 = find(measurements_3C_Q_4_R_6(:,2)~=0,1,"first")-1-differentStart;
indexLast_Q4R6 = indexFirst_Q4R6+N-1;

stateCovariance_Q4R6 = measurements_3C_Q_4_R_6(indexFirst_Q4R6:indexLast_Q4R6,5);
estimatorGain_Q4R6 = measurements_3C_Q_4_R_6(indexFirst_Q4R6:indexLast_Q4R6,6);

%%%% 3C_Q_5_R_6 %%%%
% Q = 8*10-5 & R = 8*10-6
csvfile_3C_Q_5_R_6 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_5_R_6.csv');
measurements_3C_Q_5_R_6 = dlmread(csvfile_3C_Q_5_R_6, ',', 2, 0); % Data follows the labels

indexFirst_Q5R6 = find(measurements_3C_Q_5_R_6(:,2)~=0,1,"first")-1-differentStart;
indexLast_Q5R6 = indexFirst_Q5R6+N-1;

stateCovariance_Q5R6 = measurements_3C_Q_5_R_6(indexFirst_Q5R6:indexLast_Q5R6,5);
estimatorGain_Q5R6 = measurements_3C_Q_5_R_6(indexFirst_Q5R6:indexLast_Q5R6,6);

%%%% 3C_Q_6_R_6 %%%%
% Q = 8*10-6 & R = 8*10-6
csvfile_3C_Q_6_R_6 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_6_R_6.csv');
measurements_3C_Q_6_R_6 = dlmread(csvfile_3C_Q_6_R_6, ',', 2, 0); % Data follows the labels

indexFirst_Q6R6 = find(measurements_3C_Q_6_R_6(:,2)~=0,1,"first")-1-differentStart;
indexLast_Q6R6 = indexFirst_Q6R6+N-1;

stateCovariance_Q6R6 = measurements_3C_Q_6_R_6(indexFirst_Q6R6:indexLast_Q6R6,5);
estimatorGain_Q6R6 = measurements_3C_Q_6_R_6(indexFirst_Q6R6:indexLast_Q6R6,6);

%%%% 3C_Q_7_R_6 %%%%
% Q = 8*10-7 & R = 8*10-6
csvfile_3C_Q_7_R_6 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_7_R_6.csv');
measurements_3C_Q_7_R_6 = dlmread(csvfile_3C_Q_7_R_6, ',', 2, 0); % Data follows the labels

indexFirst_Q7R6 = find(measurements_3C_Q_7_R_6(:,2)~=-0,1,"first")-1-differentStart;
indexLast_Q7R6 = indexFirst_Q7R6+N-1;

stateCovariance_Q7R6 = measurements_3C_Q_7_R_6(indexFirst_Q7R6:indexLast_Q7R6,5);
estimatorGain_Q7R6 = measurements_3C_Q_7_R_6(indexFirst_Q7R6:indexLast_Q7R6,6);

%%%% 3C_Q_8_R_6 %%%%
% Q = 8*10-8 & R = 8*10-6
csvfile_3C_Q_8_R_6 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_8_R_6.csv');
measurements_3C_Q_8_R_6 = dlmread(csvfile_3C_Q_8_R_6, ',', 2, 0); % Data follows the labels

indexFirst_Q8R6 = find(measurements_3C_Q_8_R_6(:,2)~=0,1,"first")-1-differentStart;
indexLast_Q8R6 = indexFirst_Q8R6+N-1;

stateCovariance_Q8R6 = measurements_3C_Q_8_R_6(indexFirst_Q8R6:indexLast_Q8R6,5);
estimatorGain_Q8R6 = measurements_3C_Q_8_R_6(indexFirst_Q8R6:indexLast_Q8R6,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   R = changes & Q = constant (=8*10^-6) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 3C_Q_6_R_4 %%%%
% Q = 8*10-6 & R = 8*10-4
csvfile_3C_Q_6_R_4 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_6_R_4.csv');
measurements_3C_Q_6_R_4 = dlmread(csvfile_3C_Q_6_R_4, ',', 2, 0); % Data follows the labels

indexFirst_Q6R4 = find(measurements_3C_Q_6_R_4(:,2)~=0,1,"first")-1-differentStart;
indexLast_Q6R4 = indexFirst_Q6R4+N-1;

stateCovariance_Q6R4 = measurements_3C_Q_6_R_4(indexFirst_Q6R4:indexLast_Q6R4,5);
estimatorGain_Q6R4 = measurements_3C_Q_6_R_4(indexFirst_Q6R4:indexLast_Q6R4,6);

%%%% 3C_Q_6_R_5 %%%%
% Q = 8*10-6 & R = 8*10-5
csvfile_3C_Q_6_R_5 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_6_R_5.csv');
measurements_3C_Q_6_R_5 = dlmread(csvfile_3C_Q_6_R_5, ',', 2, 0); % Data follows the labels

indexFirst_Q6R5 = find(measurements_3C_Q_6_R_5(:,2)~=0,1,"first")-1-differentStart;
indexLast_Q6R5 = indexFirst_Q6R5+N-1;

stateCovariance_Q6R5 = measurements_3C_Q_6_R_5(indexFirst_Q6R5:indexLast_Q6R5,5);
estimatorGain_Q6R5 = measurements_3C_Q_6_R_5(indexFirst_Q6R5:indexLast_Q6R5,6);

%%%% 3C_Q_6_R_6 %%%%
% Q = 8*10-6 & R = 8*10-6
% Values already

%%%% 3C_Q_6_R_7 %%%%
% Q = 8*10-6 & R = 8*10-7
csvfile_3C_Q_6_R_7 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_6_R_7.csv');
measurements_3C_Q_6_R_7 = dlmread(csvfile_3C_Q_6_R_7, ',', 2, 0); % Data follows the labels

indexFirst_Q6R7 = find(measurements_3C_Q_6_R_7(:,2)~=0,1,"first")-1-differentStart;
indexLast_Q6R7 = indexFirst_Q6R7+N-1;

stateCovariance_Q6R7 = measurements_3C_Q_6_R_7(indexFirst_Q6R7:indexLast_Q6R7,5);
estimatorGain_Q6R7 = measurements_3C_Q_6_R_7(indexFirst_Q6R7:indexLast_Q6R7,6);

%%%% 3C_Q_6_R_8 %%%%
% Q = 8*10-6 & R = 8*10-8
csvfile_3C_Q_6_R_8 = fullfile(currentDir_Assignment3_3C_Data,'Assignment3_3C_Q_6_R_8.csv');
measurements_3C_Q_6_R_8 = dlmread(csvfile_3C_Q_6_R_8, ',', 2, 0); % Data follows the labels

indexFirst_Q6R8 = find(measurements_3C_Q_6_R_8(:,2)~=0,1,"first")-1-differentStart;
indexLast_Q6R8 = indexFirst_Q6R8+N-1;

stateCovariance_Q6R8 = measurements_3C_Q_6_R_8(indexFirst_Q6R8:indexLast_Q6R8,5);
estimatorGain_Q6R8 = measurements_3C_Q_6_R_8(indexFirst_Q6R8:indexLast_Q6R8,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if assignment3_3C
    figure(1)
    subplot(1,2,1)
    hold on;
    plot(time,stateCovariance_Q4R6,'LineWidth',2)
    plot(time,stateCovariance_Q5R6,'LineWidth',2)
    plot(time,stateCovariance_Q6R6,'LineWidth',2)
    plot(time,stateCovariance_Q7R6,'LineWidth',2)
    plot(time,stateCovariance_Q8R6,'LineWidth',2)
    legend('Q=8*10^{-4}','Q=8*10^{-5}','Q=8*10^{-6}','Q=8*10^{-7}','Q=8*10^{-8}','FontSize',12);
    xlabel('Time [s]','FontSize',15)
    ylabel('$\hat{P}_{k|k}$ $[m^{2}]$', 'Interpreter', 'latex','FontSize',15)
    grid on;
    hold off;
    subplot(1,2,2)
    hold on;
    plot(time,stateCovariance_Q6R4,'LineWidth',2)
    plot(time,stateCovariance_Q6R5,'LineWidth',2)
    plot(time,stateCovariance_Q6R6,'LineWidth',2)
    plot(time,stateCovariance_Q6R7,'LineWidth',2)
    plot(time,stateCovariance_Q6R8,'LineWidth',2)
    legend('R=8*10^{-4}','R=8*10^{-5}','R=8*10^{-6}','R=8*10^{-7}','R=8*10^{-8}','FontSize',12);
    xlabel('Time [s]','FontSize',15)
    ylabel('$\hat{P}_{k|k}$ $[m^{2}]$', 'Interpreter', 'latex','FontSize',15)
    grid on;
    hold off
    
    
    figure(2)
    subplot(1,2,1)
    hold on;
    plot(time,estimatorGain_Q4R6,'LineWidth',2)
    plot(time,estimatorGain_Q5R6,'LineWidth',2)
    plot(time,estimatorGain_Q6R6,'LineWidth',2)
    plot(time,estimatorGain_Q7R6,'LineWidth',2)
    plot(time,estimatorGain_Q8R6,'LineWidth',2)
    legend('Q=8*10^{-4}','Q=8*10^{-5}','Q=8*10^{-6}','Q=8*10^{-7}','Q=8*10^{-8}','FontSize',12);
    xlabel('Time [s]','FontSize',15)
    ylabel('$L_k$ $[m^{2}]$', 'Interpreter', 'latex','FontSize',15)
    grid on;
    hold off;
    
    subplot(1,2,2)
    hold on;
    plot(time,estimatorGain_Q6R4,'LineWidth',2)
    plot(time,estimatorGain_Q6R5,'LineWidth',2)
    plot(time,estimatorGain_Q6R6,'LineWidth',2)
    plot(time,estimatorGain_Q6R7,'LineWidth',2)
    plot(time,estimatorGain_Q6R8,'LineWidth',2)
    legend('R=8*10^{-4}','R=8*10^{-5}','R=8*10^{-6}','R=8*10^{-7}','R=8*10^{-8}','FontSize',12);
    xlabel('Time [s]','FontSize',15)
    ylabel('$L_k$ $[m^{2}]$', 'Interpreter', 'latex','FontSize',15)
    grid on;
    hold off
end

%% Assignment 3 3D:
% To create a KalmanExperiment object the data files from assignment 3_3C
% need to be adjusted so that the columns have the correct name to use the
% createfromQRC3 method. These adjusted datafiles are stored in the folder
% Assignment3_3D.data

if assignment3_3D
    R = 8*10^-6;
    Q4R6 = KalmanExperiment.createfromQRC3();
    Q5R6 = KalmanExperiment.createfromQRC3();
    Q6R6 = KalmanExperiment.createfromQRC3();
    Q7R6 = KalmanExperiment.createfromQRC3();
    Q8R6 = KalmanExperiment.createfromQRC3();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q = 8*10^-4 and R = 8*10^-6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q4R6.t = time';
    Q4R6.x = Q4R6.x(indexFirst_Q4R6:indexLast_Q4R6);
    Q4R6.P = Q4R6.P(indexFirst_Q4R6:indexLast_Q4R6);
    Q4R6.R = Q4R6.R(indexFirst_Q4R6:indexLast_Q4R6);
    Q4R6.y = Q4R6.y(indexFirst_Q4R6:indexLast_Q4R6);
    Q4R6.nu = Q4R6.nu(indexFirst_Q4R6:indexLast_Q4R6);
    Q4R6.S = Q4R6.S(indexFirst_Q4R6:indexLast_Q4R6);
    
    Q4R6.analyzeconsistency();
    gcf_Q4R6 = gcf;
    % Find all axes objects in the figure
    axesHandles = findobj(gcf_Q4R6, 'Type', 'axes');
    titles = cell(1, length(axesHandles));
    % Get the title of each subplot
    for i = 1:length(axesHandles)
       titleString = get(axesHandles(i), 'Title');
       titles{i}=titleString.String;
    end
    titles{1} = strcat(titles{1}, '  [Q = 8*10^{-4} m^2 & R = 8*10^{-6} m^2]');
    titles{2} = strcat(titles{2}, '  [Q = 8*10^{-4} m^2 & R = 8*10^{-6} m^2]');
    
    axesHandles(1).Title.String = titles{1};
    axesHandles(2).Title.String = titles{2};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q = 8*10^-5 and R = 8*10^-6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q5R6.t = time';
    Q5R6.x = Q5R6.x(indexFirst_Q5R6:indexLast_Q5R6);
    Q5R6.P = Q5R6.P(indexFirst_Q5R6:indexLast_Q5R6);
    Q5R6.R = Q5R6.R(indexFirst_Q5R6:indexLast_Q5R6);
    Q5R6.y = Q5R6.y(indexFirst_Q5R6:indexLast_Q5R6);
    Q5R6.nu = Q5R6.nu(indexFirst_Q5R6:indexLast_Q5R6);
    Q5R6.S = Q5R6.S(indexFirst_Q5R6:indexLast_Q5R6);
    
    Q5R6.analyzeconsistency();
    gcf_Q5R6 = gcf;
    % Find all axes objects in the figure
    axesHandles = findobj(gcf_Q5R6, 'Type', 'axes');
    titles = cell(1, length(axesHandles));
    % Get the title of each subplot
    for i = 1:length(axesHandles)
       titleString = get(axesHandles(i), 'Title');
       titles{i}=titleString.String;
    end
    titles{1} = strcat(titles{1}, '  [Q = 8*10^{-5} m^2 & R = 8*10^{-6} m^2]');
    titles{2} = strcat(titles{2}, '  [Q = 8*10^{-5} m^2 & R = 8*10^{-6} m^2]');
    
    axesHandles(1).Title.String = titles{1};
    axesHandles(2).Title.String = titles{2};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q = 8*10^-6 and R = 8*10^-6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q6R6.t = time';
    Q6R6.x = Q6R6.x(indexFirst_Q6R6:indexLast_Q6R6);
    Q6R6.P = Q6R6.P(indexFirst_Q6R6:indexLast_Q6R6);
    Q6R6.R = Q6R6.R(indexFirst_Q6R6:indexLast_Q6R6);
    Q6R6.y = Q6R6.y(indexFirst_Q6R6:indexLast_Q6R6);
    Q6R6.nu = Q6R6.nu(indexFirst_Q6R6:indexLast_Q6R6);
    Q6R6.S = Q6R6.S(indexFirst_Q6R6:indexLast_Q6R6);
    
    Q6R6.analyzeconsistency();
    gcf_Q6R6 = gcf;
    % Find all axes objects in the figure
    axesHandles = findobj(gcf_Q6R6, 'Type', 'axes');
    titles = cell(1, length(axesHandles));
    % Get the title of each subplot
    for i = 1:length(axesHandles)
       titleString = get(axesHandles(i), 'Title');
       titles{i}=titleString.String;
    end
    titles{1} = strcat(titles{1}, '  [Q = 8*10^{-6} m^2 & R = 8*10^{-6} m^2]');
    titles{2} = strcat(titles{2}, '  [Q = 8*10^{-6} m^2 & R = 8*10^{-6} m^2]');
    
    axesHandles(1).Title.String = titles{1};
    axesHandles(2).Title.String = titles{2};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q = 8*10^-7 and R = 8*10^-6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q7R6.t = time';
    Q7R6.x = Q7R6.x(indexFirst_Q7R6:indexLast_Q7R6);
    Q7R6.P = Q7R6.P(indexFirst_Q7R6:indexLast_Q7R6);
    Q7R6.R = Q7R6.R(indexFirst_Q7R6:indexLast_Q7R6);
    Q7R6.y = Q7R6.y(indexFirst_Q7R6:indexLast_Q7R6);
    Q7R6.nu = Q7R6.nu(indexFirst_Q7R6:indexLast_Q7R6);
    Q7R6.S = Q7R6.S(indexFirst_Q7R6:indexLast_Q7R6);
    
    Q7R6.analyzeconsistency();
    gcf_Q7R6 = gcf;
    % Find all axes objects in the figure
    axesHandles = findobj(gcf_Q7R6, 'Type', 'axes');
    titles = cell(1, length(axesHandles));
    % Get the title of each subplot
    for i = 1:length(axesHandles)
       titleString = get(axesHandles(i), 'Title');
       titles{i}=titleString.String;
    end
    titles{1} = strcat(titles{1}, '  [Q = 8*10^{-7} m^2 & R = 8*10^{-6} m^2]');
    titles{2} = strcat(titles{2}, '  [Q = 8*10^{-7} m^2 & R = 8*10^{-6} m^2]');
    
    axesHandles(1).Title.String = titles{1};
    axesHandles(2).Title.String = titles{2};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q = 8*10^-8 and R = 8*10^-6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q8R6.t = time';
    Q8R6.x = Q8R6.x(indexFirst_Q8R6:indexLast_Q8R6);
    Q8R6.P = Q8R6.P(indexFirst_Q8R6:indexLast_Q8R6);
    Q8R6.R = Q8R6.R(indexFirst_Q8R6:indexLast_Q8R6);
    Q8R6.y = Q8R6.y(indexFirst_Q8R6:indexLast_Q8R6);
    Q8R6.nu = Q8R6.nu(indexFirst_Q8R6:indexLast_Q8R6);
    Q8R6.S = Q8R6.S(indexFirst_Q8R6:indexLast_Q8R6);
    
    Q8R6.analyzeconsistency();
    gcf_Q8R6 = gcf;
    % Find all axes objects in the figure
    axesHandles = findobj(gcf_Q8R6, 'Type', 'axes');
    titles = cell(1, length(axesHandles));
    % Get the title of each subplot
    for i = 1:length(axesHandles)
       titleString = get(axesHandles(i), 'Title');
       titles{i}=titleString.String;
    end
    titles{1} = strcat(titles{1}, '  [Q = 8*10^{-8} m^2 & R = 8*10^{-6} m^2]');
    titles{2} = strcat(titles{2}, '  [Q = 8*10^{-8} m^2 & R = 8*10^{-6} m^2]');
    
    axesHandles(1).Title.String = titles{1};
    axesHandles(2).Title.String = titles{2};
    
    gcf_Q4R6.Position = [350 160 800 550];
    gcf_Q5R6.Position = [350 160 800 550];
    gcf_Q6R6.Position = [350 160 800 550];
    gcf_Q7R6.Position = [350 160 800 550];
    gcf_Q8R6.Position = [350 160 800 550];
end
