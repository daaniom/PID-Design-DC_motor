relativePath = 'PendulumSwing.csv';

currentDir = pwd;


csvfile = fullfile(currentDir, relativePath);

labels = strsplit(fileread(csvfile), '\n');
labels = strsplit(labels{:, 2}, ', '); 
data = dlmread(csvfile, ',', 2, 0); 

output = data(:, 13);
output = output(1:end-250);
samplingFrequency = 100;
timeInSeconds = (0:(length(output)-1)) / samplingFrequency;


plot(timeInSeconds, output, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 14);
ylabel('\theta (rad/s)', 'FontSize', 14);

grid on;
grid minor;
startTime = 1.58;
endTime = 10.59;
xCoords = [startTime, endTime, endTime, startTime];
yCoords = [2.6, 2.6, -2.6, -2.6];
patch('XData', xCoords, 'YData', yCoords, 'EdgeColor', 'black', 'LineStyle', '-', 'FaceColor', [0 0 0], 'FaceAlpha', 0.0001,'LineWidth',2);

startIndex = find(timeInSeconds >= startTime, 1);
endIndex = find(timeInSeconds >= endTime, 1);
selectedSignal = output(startIndex:endIndex);
n = length(selectedSignal);
fftResult = fft(selectedSignal);
frequencies = (0:n-1) * (samplingFrequency/n);
[~, peakIndex] = max(abs(fftResult));


dominantFrequency = frequencies(peakIndex);

fprintf('Dominant Frequency: %.2f Hz\n', dominantFrequency*2*pi);


[peakValues, peakIndices] = findpeaks(selectedSignal, 'MinPeakHeight', 0);


positivePeakTimes = timeInSeconds(startIndex + peakIndices - 1);
positivePeakOutputs = selectedSignal(peakIndices);


disp('Positive Peak Times:');
disp(positivePeakTimes');
disp('Positive Peak Outputs:');
disp(positivePeakOutputs');

positivePeakTimes = positivePeakTimes(:);
positivePeakOutputs = positivePeakOutputs(:);

fitModel = fit(positivePeakTimes, positivePeakOutputs, 'exp1');

fittedValues = feval(fitModel, timeInSeconds(startIndex:endIndex));


% hold on
% plot(timeInSeconds(startIndex:endIndex), fittedValues, 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
% plot(timeInSeconds(startIndex:endIndex), -fittedValues, 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
% xlabel('Time (s)');
% ylabel('\theta (rad/s)');
% grid on;
% grid minor;
coefficients = coeffvalues(fitModel);

% The damping coefficient is the second coefficient
dampingCoefficient = coefficients(2);

fprintf('Damping Coefficient: %.4f\n', dampingCoefficient);

natf = 2*pi*dominantFrequency/sqrt(1-dampingCoefficient^2)
l_calc = 9.81/natf^2