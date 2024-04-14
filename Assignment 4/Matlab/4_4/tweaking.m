currentDir = pwd();
x_length=400;
x_before_rise = 50;
relativePath ='Rho17.csv'; 
csvfile = fullfile(currentDir, relativePath);
labels = strsplit(fileread(csvfile), '\n'); 
labels = strsplit(labels{:, 2}, ', '); 
data = dlmread(csvfile, ',', 2, 0); 
threshold = 0.001;
indexFirst = find(abs(data(:, 8)) > threshold, 1, 'first')-x_before_rise;
indexLast = indexFirst+x_length;
x3_meas_Rho5 = data(indexFirst :indexLast,8);
V3_meas_Rho5 = data(indexFirst :indexLast,8+2);
N_step = length(x3_meas_Rho5);
t_step_5 = [0:N_step-1]'*0.01;

relativePath ='Rho18.csv'; 
csvfile = fullfile(currentDir, relativePath);
labels = strsplit(fileread(csvfile), '\n'); 
labels = strsplit(labels{:, 2}, ', '); 
data = dlmread(csvfile, ',', 2, 0); 
threshold = 0.001;
indexFirst = find(abs(data(:, 8)) > threshold, 1, 'first')-x_before_rise;
indexLast = indexFirst+x_length;
x3_meas_Rho10 = data(indexFirst :indexLast,8);
V3_meas_Rho10 = data(indexFirst :indexLast,8+2);

relativePath ='Rho20.csv'; 
csvfile = fullfile(currentDir, relativePath);
labels = strsplit(fileread(csvfile), '\n');
labels = strsplit(labels{:, 2}, ', ');
data = dlmread(csvfile, ',', 2, 0); 
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

legend('\rho = 17', '\rho = 18', '\rho = 20' ,'FontSize', 12,'Location', 'southeast');
xlabel('Time (s)');
ylabel('Position (m)');
saveas(gcf, 'Distance_stepPendu.png');

