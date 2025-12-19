function y_filtered = aslLowPass(Y, time_span, fraction)
%--------------------------------------------------------------------------
%
%   Function to apply a low-pass Butterworth filter to input data
%
%   INPUT:  input data (1D array)                       'Y'             (double)
%           total time span in milliseconds             'time_span'     (double)
%           cutoff as a fraction of sampling frequency  'fraction'      (double)
%   OUTPUT: filtered data                               'y_filtered'    (double)
%
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

    % Define the sampling frequency
    Fs = length(Y) / (time_span / 1000); % Sampling frequency in Hz

    % Define the cutoff frequency
    Fc = fraction * Fs;

    % Ensure the cutoff frequency is within a reasonable range
    Fc = min(Fc, Fs / 2 - 1); % Ensure Fc is less than Nyquist frequency

    % Normalize the cutoff frequency
    Wn = Fc / (Fs / 2);

    % Set design parameters for the Butterworth filter
    rp = 3;   % Passband ripple (dB)
    rs = 40;  % Stopband attenuation (dB)

    if Wn==0 || Wn<0
        Wn=eps;
    end

     Wstop = Wn + 0.1; % Stopband edge frequency (normalized)


    % Calculate the initial filter order
    [n, Wn] = buttord(Wn, Wstop, rp, rs);

    % Calculate the maximum allowable filter order based on the data length
    max_order = floor(length(Y) / 3) - 1;

    % Ensure the filter order does not exceed the maximum allowable value
    n = min(n, max_order);

    % If the maximum allowable order is less than 1, set it to 1
    if n < 1
        n = 1;  % Set to minimum order of 1
    end

    % Design the Butterworth low-pass filter with the adjusted order
    [b, a] = butter(n, Wn, 'low');

    % Apply the filter to the data
    y_filtered = filtfilt(b, a, Y);
    %y_filtered = y_filtered-min(y_filtered);

end