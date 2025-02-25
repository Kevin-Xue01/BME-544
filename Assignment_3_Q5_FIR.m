rp = 1;           % Passband ripple in dB 
rs = 40;          % Stopband ripple in dB
fs = 1000;        % Sampling frequency
f = [0.3 0.5];    % Cutoff frequencies
a = [0 1];        % Desired amplitudes

dev = [10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1)]; 
[n,fo,ao,w] = firpmord(f,a,dev,fs);
b = firpm(n,fo,ao,w);

% Load and apply the filter to your data
data = readtable('BaselineWander.txt', 'Delimiter', '\t');
signal_mv = data.Sig;  % Extract signal column
time = data.Time;      % Extract time column

% Apply high-pass filtering
filtered_signal = filter(b, 1, signal_mv);

% Plot original vs. filtered signal
figure;
subplot(2,1,1);
plot(time, signal_mv);
title('Original Signal');
xlabel('Time [s]');
ylabel('Signal [mV]');

subplot(2,1,2);
plot(time, filtered_signal);
title('FIR Filtered Signal');
xlabel('Time [s]');
ylabel('Signal [mV]');

% Verify filter frequency response
figure;
freqz(b,1,fs,fs);
title('High-Pass Filter Frequency Response');

% Compute the group delay
[gd, w] = grpdelay(b,1,fs,fs);

% Plot the magnitude response
figure;
[h, w] = freqz(b,1,fs,fs);
subplot(2,1,1);
plot(w, 20*log10(abs(h))); % Convert magnitude to dB
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Lowpass Filter Frequency Response');
xlim([f(1)*fs*0.8, f(2)*fs*1.2]); % Zoom in on the cutoff
ylim([-rs-5, rp+5]); 

% Plot the group delay
subplot(2,1,2);
plot(w, gd, 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Group Delay (samples)');
title('Group Delay of the Filter');
xlim([f(1)*fs*0.8, f(2)*fs*1.2]); % Zoom in on the cutoff
hold off