% Sampling frequency
fs = 1000;
nyquist_freq = fs / 2;

% Filter specifications
Wp = 0.5 / nyquist_freq; % Passband normalized frequency
Ws = 0.3 / nyquist_freq; % Stopband normalized frequency
Rp = 1;  % Maximum passband ripple (dB)
Rs = 40; % Minimum stopband attenuation (dB)

% Design Elliptic High-Pass Filter
[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[z, p, k] = ellip(n, Rp, Rs, Wp, "high");
sos = zp2sos(z, p, k);

% Visualize frequency response
figure;
freqz(sos, 1024, fs);
title('Elliptic High-Pass Filter Response');

% Load data
data = readtable('BaselineWander.txt', 'Delimiter', '\t');

% Apply the filter to the signal
filtered_signal = sosfilt(sos, data.Sig);

% Plot the original signal
figure;
plot(data.Time, data.Sig, 'b', 'DisplayName', 'Original Signal');
xlabel('Time [s]');
ylabel('Signal [mV]');
title('Original Signal');
grid on;
hold off;

% Plot the filtered signal
figure;
plot(data.Time, filtered_signal, 'r', 'DisplayName', 'IIR Filtered Signal');
xlabel('Time [s]');
ylabel('Signal [mV]');
title('IIR Filtered Signal');
grid on;
hold off;

% Compute the group delay
[gd, f] = grpdelay(sos, 1024, fs);  % Group delay of the filter

% Plot the group delay
figure;
plot(f, gd);
xlabel('Frequency (Hz)');
ylabel('Group Delay (samples)');
title('Group Delay of the Elliptic High-Pass Filter');
grid on;
hold off;