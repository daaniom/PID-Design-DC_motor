currentDir = pwd();
x_length=400;
x_before_rise = 50;
relativePath ='Rho5.csv'; 
csvfile = fullfile(currentDir, relativePath);
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels
threshold = 0.001;
indexFirst = find(abs(data(:, 8)) > threshold, 1, 'first')-x_before_rise;
indexLast = indexFirst+x_length;
x3_meas_Rho5 = data(indexFirst :indexLast,8);
V3_meas_Rho5 = data(indexFirst :indexLast,8+2);
N_step = length(x3_meas_Rho5);
t_step_5 = [0:N_step-1]'*0.01;

relativePath ='Rho10.csv'; 
csvfile = fullfile(currentDir, relativePath);
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels
threshold = 0.001;
indexFirst = find(abs(data(:, 8)) > threshold, 1, 'first')-x_before_rise;
indexLast = indexFirst+x_length;
x3_meas_Rho10 = data(indexFirst :indexLast,8);
V3_meas_Rho10 = data(indexFirst :indexLast,8+2);

relativePath ='Rho20.csv'; 
csvfile = fullfile(currentDir, relativePath);
labels = strsplit(fileread(csvfile), '\n'); % Split file in lines
labels = strsplit(labels{:, 2}, ', '); % Split and fetch the labels (they are in line 2 of every reco
data = dlmread(csvfile, ',', 2, 0); % Data follows the labels
threshold = 0.001;
indexFirst = find(abs(data(:, 8)) > threshold, 1, 'first')-x_before_rise;
indexLast = indexFirst+x_length;
x3_meas_Rho20 = data(indexFirst :indexLast,8);
V3_meas_Rho20 = data(indexFirst :indexLast,8+2);




figure()

plot(t_step_5, x3_meas_Rho5, 'LineWidth', 1.5);
hold on
plot(t_step_5, x3_meas_Rho10, 'LineWidth', 1.5);
plot(t_step_5, x3_meas_Rho20, 'LineWidth', 1.5);

hold off

legend('\rho = 5', '\rho = 10', '\rho = 20' ,'FontSize', 12,'Location', 'southeast');
xlabel('Time (s)');
ylabel('Position (m)');
saveas(gcf, 'Distance_stepPendu.png');


figure()
plot(t_step_5, V3_meas_Rho5, 'LineWidth', 1.5);
hold on
plot(t_step_5, V3_meas_Rho10, 'LineWidth', 1.5);
plot(t_step_5, V3_meas_Rho20, 'LineWidth', 1.5);

hold off

legend('\rho = 5', '\rho = 10', '\rho = 20' ,'FontSize', 12,'Location', 'northeast');
xlabel('Time (s)');
ylabel('Voltage (V)');
saveas(gcf, 'Voltages_stepPendu.png');