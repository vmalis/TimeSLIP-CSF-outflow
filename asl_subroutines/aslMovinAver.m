function y_filtered = aslMovinAver(x, y)
%--------------------------------------------------------------------------
% Moving average filter with automatic window size estimation
%
% INPUT:  x - time points (1D array)                             (double)
%         y - signal values at time points (1D array)            (double)
%
% OUTPUT: y_filtered - smoothed signal
%--------------------------------------------------------------------------

    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);

    % Estimate average time step
    avg_time_step = mean(diff(x));

    % Estimate window size based on total duration and sampling rate
    total_duration = x(end) - x(1);
    num_points = length(x);
    approx_num_peaks = 2;  % Assuming about 2 peaks in the signal

    % Automatically estimate the window size
    window_size = floor(num_points / (total_duration * approx_num_peaks / avg_time_step));

    % Ensure window size is valid and odd (odd for better centering)
    window_size = max(3, window_size);  % Minimum window size of 3
    if mod(window_size, 2) == 0
        window_size = window_size + 1;  % Convert to odd if necessary
    end

    % Apply moving average using movmean
    y_filtered = movmean(y, window_size);
    %y_filtered=y_filtered-min(y_filtered);

end