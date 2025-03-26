close all
clear
clc

% Load MP3 file
[mp3_audio, mp3_fs] = audioread('Assignment_5_4_train_whistle_sliced.mp3');  % Replace with your file name

% Resample MP3 audio to 10 kHz
mp3_resampled = resample(mp3_audio, 10000, mp3_fs);  % Resample to 10 kHz

% Load WAV file
[wav_audio, wav_fs] = audioread('Assignment_5_4_original.wav');  % Replace with your file name

% Resample WAV audio to 10 kHz
wav_resampled = resample(wav_audio, 10000, wav_fs);  % Resample to 10 kHz

mp3_resampled = mp3_resampled(1:length(wav_resampled));

% Generate white noise with the same length as the MP3 file
white_noise = randn(length(wav_resampled), 1);

% Normalize white noise to the range [-1, 1]
white_noise = white_noise / max(abs(white_noise));

% Plot the resampled signals (optional)
figure;
subplot(3,1,1);
plot(mp3_resampled);
title('Resampled MP3 Audio at 10 kHz');

subplot(3,1,2);
plot(wav_resampled);
title('Resampled WAV Audio at 10 kHz');

subplot(3,1,3);
plot(white_noise);
title('Normalized White Noise Resampled at 10 kHz');

% Save the white noise sample as a WAV file
audiowrite('Assignment_5_4_white_noise.wav', white_noise, 10000);
%%
% Parameters for windowing and segmentation
window_length = 256;
overlap_length = window_length / 2;
window = hamming(window_length);  % Hamming window
model_order = 30;  % Initial LPC model order in the range of 20–30

% Number of segments
num_segments = floor((length(wav_resampled) - overlap_length) / (window_length - overlap_length));

% Preallocate array to store LPC coefficients for each segment
lpc_coeffs = zeros(num_segments, model_order);

% Segmenting the signal and performing LPC analysis on each segment
for i = 1:num_segments
    % Define the segment start and end indices
    start_idx = (i - 1) * (window_length - overlap_length) + 1;
    end_idx = start_idx + window_length - 1;
    
    % Extract the current segment and apply the window
    segment = wav_resampled(start_idx:end_idx) .* window;
    
    % Estimate LPC coefficients using the 'lpc' function
    [a, ~] = lpc(segment, model_order);  % 'a' are the LPC coefficients
    
    % Store the LPC coefficients
    lpc_coeffs(i, :) = a(2:end);  % a(2:end) skips the first coefficient (the gain term)
end

% Preallocate array for the filtered noise signal
filtered_noise = zeros(length(white_noise), 1);

% Filter each segment of the noise signal using the LPC coefficients
for i = 1:num_segments
    % Define the segment start and end indices for the noise signal
    start_idx = (i - 1) * (window_length - overlap_length) + 1;
    end_idx = start_idx + window_length - 1;
    
    % Extract the current segment of the white noise
    noise_segment = white_noise(start_idx:end_idx) .* window;
    
    % Get the LPC coefficients for the current segment
    a = [1, lpc_coeffs(i, :)];  % Add the '1' for the first coefficient (b = 1)
    
    % Filter the noise segment using the LPC coefficients
    filtered_segment = filter(1, a, noise_segment);  % Apply the LPC filter
    
    % Store the filtered segment in the output signal
    filtered_noise(start_idx:end_idx) = filtered_noise(start_idx:end_idx) + filtered_segment;
end

filtered_noise = filtered_noise / max(abs(filtered_noise));
% soundsc(filtered_noise, 10000);  % Play at 10 kHz sample rate
audiowrite('Assignment_5_4_vocoded_speech.wav', filtered_noise, 10000);

%%
% Parameters for windowing and segmentation
window_length = 256;
overlap_length = window_length / 2;
window = hamming(window_length);  % Hamming window
model_order = 30;  % Initial LPC model order in the range of 20–30

% Number of segments
num_segments = floor((length(mp3_resampled) - overlap_length) / (window_length - overlap_length));

% Preallocate array to store LPC coefficients for each segment
lpc_coeffs = zeros(num_segments, model_order);

% Segmenting the signal and performing LPC analysis on each segment
for i = 1:num_segments
    % Define the segment start and end indices
    start_idx = (i - 1) * (window_length - overlap_length) + 1;
    end_idx = start_idx + window_length - 1;
    
    % Extract the current segment and apply the window
    segment = mp3_resampled(start_idx:end_idx) .* window;
    
    % Estimate LPC coefficients using the 'lpc' function
    [a, ~] = lpc(segment, model_order);  % 'a' are the LPC coefficients
    
    % Store the LPC coefficients
    lpc_coeffs(i, :) = a(2:end);  % a(2:end) skips the first coefficient (the gain term)
end

% Preallocate array for the filtered noise signal
filtered_noise = zeros(length(white_noise), 1);

% Filter each segment of the noise signal using the LPC coefficients
for i = 1:num_segments
    % Define the segment start and end indices for the noise signal
    start_idx = (i - 1) * (window_length - overlap_length) + 1;
    end_idx = start_idx + window_length - 1;
    
    % Extract the current segment of the white noise
    noise_segment = white_noise(start_idx:end_idx) .* window;
    
    % Get the LPC coefficients for the current segment
    a = [1, lpc_coeffs(i, :)];  % Add the '1' for the first coefficient (b = 1)
    
    % Filter the noise segment using the LPC coefficients
    filtered_segment = filter(1, a, noise_segment);  % Apply the LPC filter
    
    % Store the filtered segment in the output signal
    filtered_noise(start_idx:end_idx) = filtered_noise(start_idx:end_idx) + filtered_segment;
end

filtered_noise = filtered_noise / max(abs(filtered_noise));
soundsc(filtered_noise, 10000);  % Play at 10 kHz sample rate
audiowrite('Assignment_5_4_vocoded_broadband.wav', filtered_noise, 10000);