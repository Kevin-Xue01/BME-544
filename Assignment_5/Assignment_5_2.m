close all
clear
clc

% Read CSV file
csvData = readtable('A5_Recordings/VowelSegments.csv', 'Delimiter', ',');

results = table('Size', [0 5], ...
                'VariableTypes', {'string', 'string', 'double', 'double', 'double'}, ...
                'VariableNames', {'UserKey', 'Word', 'Trial', 'F1', 'F2'});

% Convergence parameters
initial_order = 10; % Start with a low model order
max_order = 50; % Upper bound for the model order
tolerance = 50; % Hz tolerance for convergence

% Iterate through each row and process the audio files
for i = 1:height(csvData)

    filename = csvData.Filename{i};

    parts = split(filename, '_');
    if length(parts) < 3
        disp(filename)
        continue; % Skip malformatted filenames
    end
    
    user_key = parts{1};
    word = parts{2};
    trial = str2double(erase(parts{3}, '.wav')); % Convert trial number to numeric

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

    % Iteratively find the best model order based on peak convergence
    prev_formants = [];
    for p = initial_order:max_order
        if mod(p, 2) == 1
            continue
        end
        % Estimate AR coefficients using Burg's method
        a = arburg(slicedAudio, p);

        % Compute the Power Spectral Density using Burg's method
        [psd, freq] = pburg(slicedAudio, p, 256, fs);

        % Find peaks in PSD to estimate formants
        [peak_vals, peak_locs] = findpeaks(10*log10(psd), freq, 'MinPeakDistance', 100);
        
        % Keep only the first two formants
        if length(peak_locs) >= 2
            formants = sort(peak_locs(1:2));
        else
            continue; % Skip if not enough peaks found
        end

        % Check for convergence
        if ~isempty(prev_formants) && all(abs(formants - prev_formants) < tolerance)
            % Convergence achieved, break the loop
            break;
        end

        prev_formants = formants;
    end

    % Store final F1 and F2 after convergence
    if exist('formants', 'var') && length(formants) >= 2
        F1 = formants(1);
        F2 = formants(2);
    else
        F1 = NaN;
        F2 = NaN;
    end

    % Append results to the table
    results = [results; {user_key, word, trial, F1, F2}];
end

writetable(results, 'FormantEstimates.csv');
disp('Results saved to FormantEstimates.csv');