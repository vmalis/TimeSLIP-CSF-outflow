function [cycle_num, max_phase, percent_in_phase, max_time, percent_in_rr, max_SIR] ...
    = plot_ecg_wave(BPM, SIR, TI, time_offset_peak, axes_handle, param, ttl)
%--------------------------------------------------------------------------
%
%   Function to simulate and plot a full ECG wave with additional SIR data
%   and output details about the maximum SIR value.
%
%   Outputs:
%       cycle_num        - Number of the cardiac cycle where the max SIR occurs.
%       max_phase        - Phase ("Systolic" or "Diastolic") where the max SIR occurs.
%       percent_in_phase - Percentage of the phase duration where max SIR occurs.
%       max_time         - Time of the maximum SIR in milliseconds.
%       percent_in_rr    - Percentage within the entire RR interval where max SIR occurs.
%
%--------------------------------------------------------------------------


% Ensure the parent figure remains invisible
parent_fig = ancestor(axes_handle, 'figure'); % Get the parent figure
was_visible = strcmp(get(parent_fig, 'Visible'), 'on'); % Check visibility
if ~was_visible
    set(parent_fig, 'HandleVisibility', 'on'); % Allow updates to the figure
end


% Derived Parameters
heart_rate_hz = BPM / 60;                  % Heart rate in Hz
rr_interval = 1000 / heart_rate_hz;       % RR interval in milliseconds
sampling_rate = 1000;                     % Sampling rate (samples per second)

% Add legend
% Determine the correct interpreter based on param.latex
if param.latex
    legend_interpreter = 'latex';
else
    legend_interpreter = 'none';
end

% Determine total time from maximum TI value
total_time_ms = max(TI)+100;
t = 0:(1/sampling_rate)*1000:total_time_ms; % Time vector in milliseconds

% ECG Wave Components
p_wave = @(t, tp) 0.1 * exp(-0.5 * ((t - tp) / 20).^2); % P wave
q_wave = @(t, tq) -0.2 * exp(-0.5 * ((t - tq) / 10).^2); % Q wave
r_wave = @(t, tr) 1.0 * exp(-0.5 * ((t - tr) / 5).^2);  % R wave
s_wave = @(t, ts) -0.5 * exp(-0.5 * ((t - ts) / 10).^2); % S wave
t_wave = @(t, tt) 0.3 * exp(-0.5 * ((t - tt) / 30).^2); % T wave

% Initialize the ECG signal
ecg_wave = zeros(size(t));

% Generate ECG wave by adding components for each cycle
for peak_time = time_offset_peak:rr_interval:total_time_ms
    % Add P, Q, R, S, and T waves centered around the R peak
    ecg_wave = ecg_wave + p_wave(t, peak_time - 100); % P wave (100 ms before R peak)
    ecg_wave = ecg_wave + q_wave(t, peak_time - 30);  % Q wave (30 ms before R peak)
    ecg_wave = ecg_wave + r_wave(t, peak_time);       % R wave
    ecg_wave = ecg_wave + s_wave(t, peak_time + 20);  % S wave (20 ms after R peak)
    ecg_wave = ecg_wave + t_wave(t, peak_time + 100); % T wave centered at 100 ms after R peak
end

% Normalize ECG Wave for right Y-axis plotting
ecg_wave = ecg_wave / max(ecg_wave);

% Configure plot appearance
set(axes_handle, 'Color', 'w');
axes_handle.XColor = 'k';
axes_handle.YColor = 'k';
h1 = [];
h2 = [];
h3 = [];
h4 = [];

% Initialize variables for output
cycle_num = 0;
max_phase = "";
percent_in_phase = 0;
max_time = 0;
percent_in_rr = 0;

% Plot the ECG Wave
hold(axes_handle, 'on'); % Set the target axes for plottings

% Identify systolic and diastolic parts
systolic_mask = false(size(t));
diastolic_mask = false(size(t));

% Find the maximum SIR value and corresponding time
[max_SIR, max_idx] = max(SIR);
max_time = TI(max_idx);

for cycle = 1:length(time_offset_peak:rr_interval:total_time_ms)
    peak_time = time_offset_peak + (cycle - 1) * rr_interval;
    
    % Find the end of the T wave dynamically
    t_wave_end_idx = find(t > (peak_time + 130) & ecg_wave < 0.01, 1); % Locate when T wave ends
    if isempty(t_wave_end_idx)
        % Fallback: Approximation for T wave end if not clearly found
        t_wave_end_idx = find(t > (peak_time + 120), 1);
    end
    if isempty(t_wave_end_idx)
        continue; % Skip if T wave end cannot be found at all
    end
    
    t_wave_end = t(t_wave_end_idx); % Get the time for the end of T wave

    % Define systolic part: From R peak to the end of T wave
    systolic_start = peak_time;
    systolic_end = t_wave_end;
    systolic_mask_segment = (t >= systolic_start & t <= systolic_end);

    % Define diastolic part: The rest of the RR interval
    diastolic_start = systolic_end;
    diastolic_end = peak_time + rr_interval;
    diastolic_mask_segment = (t > diastolic_start & t <= diastolic_end);

    % Check if max SIR occurs in this cycle
    if systolic_start <= max_time && max_time <= systolic_end
        max_phase = "Systolic";
        cycle_num = cycle;
        percent_in_phase = ((max_time - systolic_start) / (systolic_end - systolic_start)) * 100;
    elseif diastolic_start <= max_time && max_time <= diastolic_end
        max_phase = "Diastolic";
        cycle_num = cycle;
        percent_in_phase = ((max_time - diastolic_start) / (diastolic_end - diastolic_start)) * 100;
    end

    % Calculate percent inside RR interval
    if peak_time <= max_time && max_time <= (peak_time + rr_interval)
        percent_in_rr = ((max_time - peak_time) / rr_interval) * 100;
    end

    % Plot systolic part in green for the current segment
    yyaxis right;
    if isempty(h1)
        h1 = plot(t(systolic_mask_segment), ecg_wave(systolic_mask_segment), '-', 'Color', [0.3, 0.8, 0.3], 'LineWidth', 2, 'DisplayName', 'Systolic');
    else
        plot(t(systolic_mask_segment), ecg_wave(systolic_mask_segment), '-', 'Color', [0.3, 0.8, 0.3], 'LineWidth', 2);
    end

    % Plot diastolic part in blue for the current segment
    yyaxis right;
    if isempty(h2)
        h2 = plot(t(diastolic_mask_segment), ecg_wave(diastolic_mask_segment), '-', 'Color', [0.3, 0.3, 0.8], 'LineWidth', 2, 'DisplayName', 'Diastolic');
    else
        plot(t(diastolic_mask_segment), ecg_wave(diastolic_mask_segment), '-', 'Color', [0.3, 0.3, 0.8], 'LineWidth', 2);
    end
end

% Plot SIR data on the left Y-axis after ECG plotting
yyaxis left;
h3 = plot(TI, SIR, 'r.', 'MarkerSize', 8, 'DisplayName', 'SIR Data');
h4 = plot(TI(max_idx), max_SIR, 'pentagram', 'Color', 'k', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Max SIR');



% Reorder elements to bring SIR plots to the front
set(axes_handle, 'Children', [h4; h3; get(axes_handle, 'Children')]);

% Configure SIR axis label based on param.latex
if param.latex
    ylabel('$ \% \; \mathrm{Signal \; Increase \; Ratio}$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k');
else
    ylabel('% Signal Increase Ratio', 'Interpreter', 'none', 'FontSize', 14, 'Color', 'k');
end

if max_SIR>param.SIRLim
    sir_lim=ceil(max_SIR / 5) * 5;
    ylim([-abs(sir_lim), sir_lim]);
else
    ylim([-param.SIRLim, param.SIRLim]);
end

% Set left Y-axis tick label color to black
yyaxis left;
set(axes_handle, 'YColor', 'k');

% Configure right Y-axis for ECG
yyaxis right;
if param.latex
    ylabel('$ \mathrm{ECG}$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k');
else
    ylabel('ECG', 'Interpreter', 'none', 'FontSize', 14, 'Color', 'k');
end
ylim([-1, 1]);
set(axes_handle, 'YColor', 'k');

% Configure axis labels and title based on param.latex
if param.latex
    xlabel(axes_handle, '$\mathrm{TI}$ time [$\mathrm{ms}$]', 'Interpreter', 'latex', 'FontSize', 10);
    set(axes_handle, 'TickLabelInterpreter', 'latex');
else
    xlabel(axes_handle, 'time [ms]', 'Interpreter', 'none', 'FontSize', 10);
    set(axes_handle, 'TickLabelInterpreter', 'none');
end

xlim([0, total_time_ms]);

grid on;
box on;

% Add annotations for clarity
peak_times = time_offset_peak:rr_interval:total_time_ms;
for peak = peak_times
    xline(peak, '--k', 'LineWidth', 1, 'Label', 'R Peak', 'LabelOrientation', 'horizontal','Interpreter',legend_interpreter);
end

% Configure x-axis ticks to be every 100 ms
x_ticks = 0:100:total_time_ms; % Generate ticks from 0 to the total time in steps of 100
set(gca, 'XTick', x_ticks);

% Add the legend with the correct interpreter
legend([h1, h2, h3, h4], {'Systolic', 'Diastolic', 'SIR Data', 'Max SIR'}, 'Location', 'eastoutside', 'Interpreter', legend_interpreter);

% Determine the title string and interpreter
if param.latex
    % Convert spaces to LaTeX-compatible spacing using \;
    titleFont = ['$\mathbf{' strrep(ttl, ' ', '\;') '}$'];
    titleInterpreter = 'latex';
else
    titleFont = ttl; % Use plain text
    titleInterpreter = 'none';
end

% Set the title with right alignment
titleHandle = title(gca, titleFont, 'Interpreter', titleInterpreter, ...
    'FontSize', 12, 'HorizontalAlignment','right');

% Adjust position to align with the right edge of the axes
titlePosition = get(titleHandle, 'Position');
titlePosition(1) = gca().XLim(2); % Align to the right edge
set(titleHandle, 'Position', titlePosition);

hold off;
end
