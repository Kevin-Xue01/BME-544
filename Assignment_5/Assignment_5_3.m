close all
clear
clc

% Load the results from the saved CSV file
results = readtable('FormantEstimates.csv');

% Filter for user "GD391P" and words "Head" and "Had"
filtered_results = results(strcmp(results.UserKey, 'GD391P') & ...
                           (strcmp(results.Word, 'head') | strcmp(results.Word, 'had')), :);

% Compute the average F1 and F2 for each word
mean_F1_Had = mean(filtered_results.F1(strcmp(filtered_results.Word, 'had')));
mean_F1_Head = mean(filtered_results.F1(strcmp(filtered_results.Word, 'head')));
mean_F2_Had = mean(filtered_results.F2(strcmp(filtered_results.Word, 'had')));
mean_F2_Head = mean(filtered_results.F2(strcmp(filtered_results.Word, 'head')));

% Compute the differences
F1_difference = mean_F1_Had - mean_F1_Head;
F2_difference = mean_F2_Had - mean_F2_Head;

% Display results
fprintf('Average F1 of "had" = %.2f Hz\n', mean_F1_Had);
fprintf('Average F1 of "head" = %.2f Hz\n', mean_F1_Head);
fprintf('F1 Difference (had - head) = %.2f Hz\n\n', F1_difference);

fprintf('Average F2 of "had" = %.2f Hz\n', mean_F2_Had);
fprintf('Average F2 of "head" = %.2f Hz\n', mean_F2_Head);
fprintf('F2 Difference (had - head) = %.2f Hz\n', F2_difference);

%%
% Define sampling frequency
fs = 10000; % Example sampling rate, adjust as needed

% Define F1 and F2 (in Hz)
F1 = mean_F1_Head;  % Example formant 1 frequency
F2 = mean_F2_Head; % Example formant 2 frequency

% Define frequency shifts (∆F1, ∆F2)
delta_F1 = abs(F1_difference);  % Example shift for F1
delta_F2 = abs(F2_difference);  % Example shift for F2

% Convert frequencies to normalized angular frequencies
theta_F1 = 2 * pi * F1 / fs;
theta_F2 = 2 * pi * F2 / fs;
theta_F1_shifted = 2 * pi * (F1 + delta_F1) / fs;
theta_F2_shifted = 2 * pi * (F2 + delta_F2) / fs;

% Define zeros (unit circle)
z = [exp(1j * theta_F1), exp(-1j * theta_F1), exp(1j * theta_F2), exp(-1j * theta_F2)];

% Define poles (radius 0.98)
r = 0.98; % Pole radius
p = [r * exp(1j * theta_F1_shifted), r * exp(-1j * theta_F1_shifted), ...
     r * exp(1j * theta_F2_shifted), r * exp(-1j * theta_F2_shifted)];

% Convert zeros and poles to transfer function coefficients
[b, a] = zp2tf(z', p', 1);  % Only two outputs needed

% Plot pole-zero diagram
figure;
zplane(b, a);
title('Pole-Zero Plot of Designed IIR Filter');
grid on;

% Compute and plot frequency response
figure;
freqz(b, a, 1024, fs);
title('Magnitude and Phase Response of the IIR Filter');

%%
% Read CSV file
csvData = readtable('A5_Recordings/VowelSegments.csv', 'Delimiter', ',');

% Define directory for filtered output files
outputDir = 'Filtered_Recordings';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Iterate through each row and process the audio files
for i = 1:height(csvData)

    % Extract filename from CSV
    filename = csvData.Filename{i};

    % Split filename to extract user key, word, and trial number
    parts = split(filename, '_');
    if length(parts) < 3
        continue; % Skip malformatted filenames
    end
    
    user_key = parts{1};
    word = parts{2};
    trial = str2double(erase(parts{3}, '.wav')); % Convert trial number to numeric

    % Filter condition: Only process "head" spoken by user "GD391P"
    if ~strcmp(user_key, 'GD391P') || ~strcmp(word, 'head')
        continue;
    end

    % Read the corresponding audio file
    [inputAudio, fs] = audioread(filename);

    % Convert stereo to mono if needed
    if size(inputAudio, 2) > 1
        inputAudio = mean(inputAudio, 2);
    end

    % Apply the IIR filter
    filteredAudio = filter(b, a, inputAudio);

    % Normalize to avoid clipping (rescale to [-1,1])
    filteredAudio = filteredAudio / max(abs(filteredAudio));

    % Generate output filename
    outputFilename = sprintf('%s/%s_%s_%d_filtered.wav', outputDir, user_key, word, trial);

    % Save filtered audio
    audiowrite(outputFilename, filteredAudio, fs);

    % Compute PSD for original audio using Burg's method
    model_orders = [14, 22, 16];
    [psdOriginal, freqOriginal] = pburg(inputAudio, model_orders(trial), 256, fs);  % Using 24th order AR model
    psdOriginal = 10*log10(psdOriginal);  % Convert to dB

    % Compute PSD for filtered audio using Burg's method
    [psdFiltered, freqFiltered] = pburg(filteredAudio, model_orders(trial), 256, fs);  % Using 24th order AR model
    psdFiltered = 10*log10(psdFiltered);  % Convert to dB

    % Plot PSD comparison
    figure;
    hold on;
    plot(freqOriginal, psdOriginal, 'b', 'LineWidth', 1.5);  % PSD of original
    plot(freqFiltered, psdFiltered, 'r', 'LineWidth', 1.5);  % PSD of filtered
    hold off;

    title(sprintf('PSD Comparison: %s %s %d', user_key, word, trial));
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend('Original PSD', 'Filtered PSD');
    grid on;

    % Display progress
    fprintf('Processed and saved: %s\n', outputFilename);
end

disp('Filtering complete. Processed files are saved in the "Filtered_Recordings" directory.');