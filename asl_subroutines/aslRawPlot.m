function varargout = aslRawPlot(ti, sir_mean, param, ttl, varargin)
%--------------------------------------------------------------------------
%
%   Subroutine to plot ASL Raw data
%
%   INPUT:  ti (time intervals)                                 'double'
%           sir_mean (mean signal increase ratio)               'double'
%           param (parameters for plot)                         'struct'
%           ttl (title for the plot)                            'string'
%           sir_pp (per voxel signal increase ratio, optional)  'double'
%           ax (axes handle, optional)                          'handle'
%
%   OUTPUT: y_filtered for sir_mean y_filtered2 for sir_pp
%
%   The function creates a plot of ASL Raw data either in a specified
%   axes or in a new figure if no axes handle is provided. The plot includes
%   the mean signal increase ratio, optional per voxel signal increase ratio,
%   and a low-pass filtered version of the mean signal increase ratio.
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------
   
    sir_pp=[];
    ax = [];
     % Adjust parameters as needed 0.1 (default) strong filter, 0.5 weak
    filter=0.1;

    % Parse optional parameters
    if ~isempty(varargin)
        for i = 1:length(varargin)
            if all(ishandle(varargin{i})) && all(strcmp(get(varargin{i}, ...
                    'type'), 'axes'))
                ax = varargin{i};
            elseif isnumeric(varargin{i})
                sir_pp = varargin{i};
            end
        end
    end
    
    % Check if sir_pp is provided
    has_sir_pp = ~isempty(sir_pp);
    
    % Check if ax is provided and is a valid axes handle
    if isempty(ax) || ~ishandle(ax) || ~strcmp(get(ax, 'type'), 'axes')
        % Create a new figure if no valid axes handle is provided
        fg = figure('Visible', param.visible);
        set(fg, 'Position', [100, 100, 500, 500]);
        set(fg, 'Color', 'w');
        ax = axes(fg);
    end
    
    x = ti;
    y = sir_mean;

    % Plot the original data and the fitted curve
    hold(ax, 'on');
    % Set the axes background to white
    set(ax, 'Color', 'w');
    y(isnan(y))=0;
    y(isinf(y))=0;
    if has_sir_pp
        y2=sir_pp;
        y2(isnan(y2))=0;
        y2(isinf(y2))=0;
        %y2=y2-min(y2);
        plot(ax, x, y2*100, 'ko', 'LineWidth', 1, 'MarkerSize', 5, ...
            'MarkerFaceColor', 'k', 'LineWidth', 1);
        switch param.filter_lowpass
        case 'movinAver'
            y_filtered_pp = aslMovinAver(ti, y2)';
        otherwise
            y_filtered_pp = aslLowPass(y2, max(x), filter);
        end

        plot(ax, x, y_filtered_pp*100, 'k--o', 'LineWidth', 1, ...
            'MarkerSize', 3, 'MarkerFaceColor', 'w', 'LineWidth', 1);
    end

    %y=y-min(y);
    plot(ax, x, y*100, 'ks','LineWidth', 1, 'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor',  'k', 'MarkerSize', 5, 'LineWidth', 1);
    
    switch param.filter_lowpass
        case 'movinAver'
            y_filtered = aslMovinAver(ti, y)';
        otherwise
            y_filtered = aslLowPass(y, max(x), filter);
    end

    plot(ax, x, y_filtered*100, 'k-s', 'MarkerFaceColor','w', ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 3,'LineWidth', 1);
    
    hold(ax, 'off');
    axis(ax, 'square');
    
    % Determine the interpreter type based on param.latex
    if param.latex
        interpreter = 'latex';
    else
        interpreter = 'tex';
    end
    
    % Determine the y-axis label based on the interpreter
    title(ax, ttl, 'Interpreter', interpreter, 'FontSize', 14, 'Color', 'k');
    if param.latex
        ylabel(ax, '$\% \; \mathrm{Signal\, Increase\, Ratio}$', 'Interpreter', ...
            'latex', 'FontSize', 10);
        xlabel(ax, '$\mathrm{TI}$ time [$\mathrm{ms}$]', 'Interpreter', ...
            interpreter, 'FontSize', 10);
    else
        ylabel(ax, '% Signal Increase Ratio', 'Interpreter', 'none', ...
            'FontSize', 10);
        xlabel(ax, 'TI time [ms]', 'Interpreter', interpreter, 'FontSize',...
            10);
    end
    
    % Set legend with the specified interpreter
    if has_sir_pp
        legend_labels = {'per voxel', 'per voxel low pass','mean', ...
            'mean low pass'};
    else
        legend_labels = {'mean', 'mean low pass'};
    end
    h_legend = legend(ax, legend_labels, 'Interpreter', interpreter, ...
        'FontSize', 10, 'Location', 'northeast');
    legend('boxoff');
    set(h_legend, 'TextColor', 'k', 'Box', 'on', 'Color', 'w', ...
        'EdgeColor', 'k');

    % Set other axis properties
    set(ax, 'TickLabelInterpreter', interpreter);

    % Calculate ymax for setting y-axis limits
    ymax = max(max(110 * [y, y_filtered, has_sir_pp * sir_pp], [], 'all'), param.SIRLim);

    if isfield(param,'TILim')
        xlim(ax, [0, param.TILim]);
    else
        xlim(ax, [0, max(x)*1.01]);
    end

    ylim(ax, [0, ymax]);

    yticks(0:10:max(ylim));
    % Set all axis and box properties
    set(ax, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'TickDir', 'in', ...
        'TickLength', [0.03, 0.015], 'LineWidth', 1);
    
    % Set grid color to gray
    set(ax, 'GridColor', [0.5, 0.5, 0.5]);

    % Set font color for title, tick labels, legend, and axis labels to black
    set(ax, 'FontSize', 10, 'FontWeight', 'normal');

    % Add minor ticks
    set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', 'MinorGridLineStyle', '-');
    
    box(ax, 'on');

    varargout{1}=y_filtered;
    if param.perVoxelCalc
        varargout{2}=y_filtered_pp;
    end

end