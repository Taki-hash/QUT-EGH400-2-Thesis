function [THD_ratio, THD_dB, THD_per] = computeTHD(y, Fs, fund_freq, num_harmonics)
% computeTHD - Manually computes Total Harmonic Distortion (THD)
%
% Inputs:
%   y             - Output signal (time-domain)
%   Fs            - Sampling frequency [Hz]
%   fund_freq     - Fundamental frequency of the input signal [Hz]
%   num_harmonics - Number of harmonics to include (e.g., 5)
%
% Outputs:
%   THD_ratio - THD as a linear ratio (e.g., 0.02 for 2%)
%   THD_dB    - THD in decibels
%   THD_per   - THD as a percentage

    N = length(y);
    Y = abs(fft(y));
    Y = Y(1:floor(N/2));   % Keep positive frequencies only
    f = (0:floor(N/2)-1) * Fs / N;

    % Find index of fundamental frequency
    [~, idx_fund] = min(abs(f - fund_freq));
    V1 = Y(idx_fund);  % Magnitude of fundamental

    % Sum magnitudes of selected harmonics
    Vh_sq_sum = 0;
    for h = 2:(num_harmonics+1)
        harmonic_freq = h * fund_freq;
        [~, idx_h] = min(abs(f - harmonic_freq));
        if idx_h <= length(Y)
            Vh_sq_sum = Vh_sq_sum + Y(idx_h)^2;
        end
    end

    if V1 < 1e-12  % or a suitable noise threshold
        THD_ratio = NaN;
        THD_dB = NaN;
        THD_per = NaN;
    else
        THD_ratio = sqrt(Vh_sq_sum) / V1;
        THD_dB = 20 * log10(THD_ratio);
        THD_per = 100 * THD_ratio;
    end


end
