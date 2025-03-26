close all
clear
clc
% Read CSV file
csvData = readtable('A5_Recordings/VowelSegments.csv', 'Delimiter', ',');

% Define search range in milliseconds
minLagMs = 3.33; % Max F0 (300 Hz)
maxLagMs = 13.3; % Min F0 (75 Hz)

% Convert time-domain constraints to frequency domain
minFreq = 1000 / maxLagMs; % Minimum frequency corresponding to maximum lag (Hz)
maxFreq = 1000 / minLagMs; % Maximum frequency corresponding to minimum lag (Hz)

words = {}; % List of unique words
wordsF0_ACF = containers.Map(); % Store ACF F0 values for each word
wordsF0_PSD = containers.Map(); % Store PSD F0 values for each word

% Iterate through each row and process the audio files
for i = 1:height(csvData)

    filename = csvData.Filename{i};
    if ~contains(filename, '200008')
        continue; % Skip this iteration and move to the next i
    end

    parts = split(filename, '_');
    if length(parts) < 3
        continue; % Skip malformatted filenames
    end
    
    word = parts{2};
    
    % Add word to list of unique words if it's not already there
    if ~ismember(word, words)
        words{end+1} = word;
        wordsF0_ACF(word) = [];
        wordsF0_PSD(word) = [];
    end

    startIdx = csvData.StartIdx(i);
    stopIdx = csvData.StopIdx(i);
    
    % Read the audio file
    [audioData, fs] = audioread(filename);
    
    % Convert indices to sample indices and slice audio
    startSample = max(1, startIdx);
    stopSample = min(length(audioData), stopIdx);
    slicedAudio = audioData(startSample:stopSample, :);
    
    % Convert to mono if stereo
    if size(slicedAudio, 2) > 1
        slicedAudio = mean(slicedAudio, 2); % Convert to mono by averaging channels
    end
    
    % Create a figure with 2 subplots
    figure('Position', [100, 100, 1000, 800]);
    
    subplot(2,1,1);
    
    % Compute autocorrelation
    r = xcorr(slicedAudio, 'unbiased');
    r = r(floor(length(r)/2):end); % Keep only positive lags
    
    % Normalize so that r(0) = 1
    r = r / r(1);
    
    % Convert lag indices to time in milliseconds
    lagTimeMs = (0:length(r)-1) * 1000 / fs;
    
    % Find the indices corresponding to search range
    minLagIdx = find(lagTimeMs >= minLagMs, 1);
    maxLagIdx = find(lagTimeMs <= maxLagMs, 1, 'last');
    
    [acfPeakVal, peakIdx] = max(r(minLagIdx:maxLagIdx));
    lastPeakIdx = peakIdx + minLagIdx - 1;
    periodMs_ACF = lagTimeMs(lastPeakIdx);
    f0_ACF = 1000 / periodMs_ACF;

    
    % Plot the autocorrelation function with peak marker
    plot(lagTimeMs, r, 'b-');
    hold on;
    
    % Add search range markers
    xline(minLagMs, 'g--', 'Min Search');
    xline(maxLagMs, 'g--', 'Max Search');
    
    % Add peak marker
    plot(lagTimeMs(lastPeakIdx), r(lastPeakIdx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(lagTimeMs(lastPeakIdx), r(lastPeakIdx), sprintf('  Last Peak: %.2f ms\n  F0: %.2f Hz', periodMs_ACF, f0_ACF), 'VerticalAlignment', 'bottom');
    
    % Add labels and title
    title(sprintf('Method 1: Autocorrelation (Last Peak) - %s', filename));
    xlabel('Lag Time (ms)');
    ylabel('Normalized ACF');
    xlim([0, min(50, max(lagTimeMs))]);  % Limit x-axis for better visibility
    grid on;
    
    % Add a legend
    legend('ACF', 'Min Search (3.5ms)', 'Max Search (11.0ms)', 'Last Peak');
    
    subplot(2,1,2);
    
    % Use Welch's method with default arguments
    % This uses 8 segments with 50% overlap by default
    segment_length = 256; % Use the entire signal
    nfft = 512; % Large FFT size
    [pxx, f] = pwelch(slicedAudio, segment_length, 128, nfft, fs);
    
    % Find the indices corresponding to search range in frequency domain
    minFreqIdx = find(f >= minFreq, 1);
    maxFreqIdx = find(f <= maxFreq, 1, 'last');
    
    % Fallback to max if no peaks found
    [psdPeakVal, peakIdx] = max(pxx(minFreqIdx:maxFreqIdx));
    firstPeakIdx = peakIdx + minFreqIdx - 1;
    f0_PSD = f(firstPeakIdx);
       
    % Store F0 values for this word
    wordsF0_ACF(word) = [wordsF0_ACF(word), f0_ACF];
    wordsF0_PSD(word) = [wordsF0_PSD(word), f0_PSD];
    
    % Plot the PSD with peak marker
    plot(f, 10*log10(pxx), 'b-');
    hold on;
    
    % Add search range markers
    xline(minFreq, 'g--', 'Min Freq');
    xline(maxFreq, 'g--', 'Max Freq');
    
    % Add peak marker
    plot(f(firstPeakIdx), 10*log10(pxx(firstPeakIdx)), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(f(firstPeakIdx), 10*log10(pxx(firstPeakIdx)), sprintf('  First Peak: %.2f Hz', f0_PSD), 'VerticalAlignment', 'bottom');
    
    % Add labels and title
    title(sprintf('Method 2: Power Spectral Density (First Peak) - %s', filename));
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    xlim([0, 500]);  % Limit x-axis for better visibility
    grid on;
    
    % Add a legend
    legend('PSD', sprintf('Min Freq (%.1f Hz)', minFreq), sprintf('Max Freq (%.1f Hz)', maxFreq), 'First Peak');
    
    % Adjust figure layout
    sgtitle(sprintf('Fundamental Frequency Estimation (Highpass filtered) - %s', filename), 'FontSize', 14);
end

% Create a table to store the results
resultTable = table('Size', [length(words), 5], 'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Word', 'ACF_Mean', 'ACF_Std', 'PSD_Mean', 'PSD_Std'});

for i = 1:length(words)
    word = words{i};
    acf_values = wordsF0_ACF(word);
    psd_values = wordsF0_PSD(word);
    
    % Calculate mean and std for both methods
    acf_mean = mean(acf_values);
    acf_std = std(acf_values);
    psd_mean = mean(psd_values);
    psd_std = std(psd_values);
    
    % Store in the table
    resultTable.Word(i) = word;
    resultTable.ACF_Mean(i) = acf_mean;
    resultTable.ACF_Std(i) = acf_std;
    resultTable.PSD_Mean(i) = psd_mean;
    resultTable.PSD_Std(i) = psd_std;
    
    fprintf('%-10s,%-12.2f,%-12.2f,%-12.2f,%-12.2f\n', word, acf_mean, acf_std, psd_mean, psd_std);
end