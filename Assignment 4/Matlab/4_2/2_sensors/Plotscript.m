




Ts = 0.01;
N = length(obj);
t = [0:N-1]*Ts;

fontSize=12;


figSize = [100, 100,800, 800];


figure('Position', figSize);

subplot(3, 1, 1);
plotstates(obj, 1);
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Position (m)', 'FontSize', fontSize);
ylim([-0.1, 0.6]);
legend('State \zeta_1 with 95% confidence interval', 'FontSize', 12); % Adjust the FontSize as needed



subplot(3, 1, 2);
plotstates(obj, 2);
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Theta (rad)', 'FontSize', fontSize);
ylim([-0.5, 0.5]);
legend('State \zeta_2 with 95% confidence interval', 'FontSize', 12); % Adjust the FontSize as needed




subplot(3, 1, 3);
plotstates(obj, 3);
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Pendulum tip speed (m/s)', 'FontSize', fontSize);
ylim([-1.3, 1.3]);
legend('State \zeta_3 with 95% confidence interval', 'FontSize', 12); % Adjust the FontSize as needed
saveas(gcf, 'Linear_Kalman_2sensors.png');



figSize = [100, 100,800, 800*2/3];

figure('Position', figSize);
subplot(2, 1, 1);
plotstates(obj, 1);
hold on
plotmeasurements(obj, 1);
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Position (m)', 'FontSize', fontSize);
hold off
ylim([-0.1, 0.6]);
legend('State \zeta_1 with 95% confidence interval','Measured state \zeta_1 with 95% confidence interval', 'FontSize', 12); % Adjust the FontSize as needed



subplot(2, 1, 2);
plotstates(obj, 2);
hold on
plotmeasurements(obj, 2);

hold off
ylim([-0.5, 0.5]);
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Theta (rad)', 'FontSize', fontSize);
legend('State \zeta_2 with 95% confidence interval','Measured state \zeta_2 with 95% confidence interval', 'FontSize', 12); % Adjust the FontSize as needed
saveas(gcf, 'Linear_Kalman_2sensors_measured.png');


