function [] = aslFitPlot(measured, GKM, Gamma, GammaF, Gauss, Combined, CombinedComponents, GOF, Ratio, param, ttl, varargin)
%--------------------------------------------------------------------------
%
%   Function to plot ASL fit data for Gaussian and gamma variate models
%
%   INPUT:  
%       measured (n by 2 array)             - x and y data points (double)
%       Gamma (n by 2 array)                - x and y data points for Gamma (double)
%       Gauss (n by 2 array)                - x and y data points for Gauss (double)
%       Combined (n by 2 array)             - x and y data points for Combined model (double)
%       GOF (1 by 3 array)                  - R^2 values for Gamma, Gauss, Combined (double)
%       Ratio (double)                      - Fractions of Gauss and Gamma in Combined model (double)
%       param (struct)                      - Parameter structure (struct)
%       ttl                                 - Title string
%       varargin                            - Optional arguments (axes handle)
%
%--------------------------------------------------------------------------

    % Parse optional axes handle from varargin
    ax = [];
    if ~isempty(varargin)
        for i = 1:length(varargin)
            if ishandle(varargin{i}) && strcmp(get(varargin{i}, 'type'), 'axes')
                ax = varargin{i};
            end
        end
    end

    % Default MATLAB colors
    gkmColor = [0, 0, 0];                       % Black for GKM
    gammaColor = [0, 0.4470, 0.8410];           % light Blue for Gamma
    gammaColorF = [1.0, 0.7, 0.5];              % dark Blue for Gamma fixed T1
    gaussColor = [0.4660, 0.6740, 0.1880];      % Green for Gauss
    combinedColor = [0.8500, 0.3250, 0.0980];   % Red for Combined

    % Determine the maximum R^2 value and corresponding line width
    [~, maxIdx] = max(GOF);
    lineWidths = 2*ones(1, 5);
    %lineWidths(maxIdx) = 2;

    % Check if ax is provided and is a valid axes handle
    if isempty(ax) || ~ishandle(ax) || ~strcmp(get(ax, 'type'), 'axes')
        % Create a new figure if no valid axes handle is provided
        fg = figure('Visible', param.visible);
        set(fg, 'Position', [100, 100, 800, 800]);
        set(fg, 'Color', 'w');
        ax = axes(fg);
    end
    
    hold(ax, 'on');
    set(ax, 'Color', 'w');

    % Flags for which plots to include based on param.fitplot
    gkm_flag = any(strcmp(param.fitplot, 'gkm')); % Gamma model
    gm_flag  = any(strcmp(param.fitplot, 'gm')); % Gamma model
    gmF_flag = any(strcmp(param.fitplot, 'gmf')); % Gamma fixed T1 model
    gs_flag  = any(strcmp(param.fitplot, 'gs')); % Gauss model
    bi_flag  = any(strcmp(param.fitplot, 'bi')); % Combined model
    biC_flag  = any(strcmp(param.fitplot, 'biC')); % Combined model components

    % Initialize plot handles and legend labels
    plotHandles = [];
    legendLabels = {};

    % Plot the measured data with circles (always plotted) 
    p4 = plot(ax, measured(:, 1), 100*measured(:, 2), 'ko', 'MarkerFaceColor',...
        'k', 'MarkerSize', 5);
    plotHandles(end+1) = p4;
    legendLabels{end+1} = 'measured';


    % Plot the GKM fit if specified
    if gkm_flag
        p0 = plot(ax, GKM(:, 1), 100*GKM(:, 2), '-', 'Color', gkmColor,...
            'LineWidth', lineWidths(1));
        plotHandles(end+1) = p0;
        if param.latex
            gkmLabel = strcat('$\mathrm{GKM} \, (R^2=$', sprintf('%.2f', GOF(1)),')');
        else
            gkmLabel = sprintf('GKM (R^2 = %.2f) ', GOF(1));
        end
        legendLabels{end+1} = gkmLabel;
    end


    % Plot the Gamma fit if specified
    if gm_flag
        p1 = plot(ax, Gamma(:, 1), 100*Gamma(:, 2), '-', 'Color', gammaColor,...
            'LineWidth', lineWidths(2));
        plotHandles(end+1) = p1;
        if param.latex
            gammaLabel = strcat('$\mathrm{\Gamma} \, (R^2=$', sprintf('%.2f', GOF(2)),')');
        else
            gammaLabel = sprintf('Gamma (R^2 = %.2f) ', GOF(2));
        end
        legendLabels{end+1} = gammaLabel;
    end

    % Plot the Gamma fit with fixed T1 if specified
    if gmF_flag
        p1f = plot(ax, GammaF(:, 1), 100*GammaF(:, 2), '-', 'Color', gammaColorF,...
            'LineWidth', lineWidths(3));
        plotHandles(end+1) = p1f;
        if param.latex
            gammaFLabel = strcat('$\mathrm{\Gamma}|_{T1_{\mathrm{fixed}}} \, (R^2=$', sprintf('%.2f', GOF(3)),')');
        else
            gammaFLabel = sprintf('Gamma|_{T1 fixed} (R^2 = %.2f) ', GOF(3));
        end
        legendLabels{end+1} = gammaFLabel;
    end



    % Plot the Gauss fit if specified
    if gs_flag
        p2 = plot(ax, Gauss(:, 1), 100*Gauss(:, 2), '-', 'Color', gaussColor,...
            'LineWidth', lineWidths(4));
        plotHandles(end+1) = p2;
        if param.latex
            gaussLabel = strcat('$\mathrm{Gaussian} \, (R^2=$', sprintf('%.2f', GOF(4)),')');
        else
            gaussLabel = sprintf('Gaussian (R^2 = %.2f) ', GOF(4));
        end
        legendLabels{end+1} = gaussLabel;
    end

    % Include spacer if necessary
    includeSpacer = ((gm_flag || gs_flag) && bi_flag);
    if includeSpacer
        sp1 = plot(nan, nan, 'w'); % Invisible spacer
        plotHandles(end+1) = sp1;
        legendLabels{end+1} = ''; % Empty string as spacer
    end

    % Plot the Combined fit if specified
    if bi_flag
        p3 = plot(ax, Combined(:, 1), 100*Combined(:, 2), '-', 'Color',...
            combinedColor, 'LineWidth', lineWidths(5));
        plotHandles(end+1) = p3;
        if param.latex
            bicmp = strcat('$\mathrm{bi-component} \, (R^2=', sprintf('%.2f', GOF(5)),') \; $');
            %frac = strcat('$f_{\mathrm{Gaus.}} = ', sprintf('%.1f', CombinedFractions(1)),'\% \quad f_{\mathrm{GKM}} = ',...
            %    sprintf('%.1f', CombinedFractions(2)),'\% \, $');
            frac = strcat('$\mathrm{SIR}_{\mathrm{Gaus.}} / \mathrm{SIR}_{\mathrm{\Gamma}} = ', sprintf('%.1f', Ratio),'$');
            combinedLabel = sprintf('%s\n%s', bicmp, frac);
        else
            %combinedLabel = strcat(sprintf('bi-component (R^2 = %.2f) ', GOF(3)),...
            %    '\newlinef_{Gaus.}=', sprintf('%.1f', CombinedFractions(1)),'%',...
            %    '  f_{GKM}=', sprintf('%.1f', CombinedFractions(2)),'%');
            combinedLabel = strcat(sprintf('bi-component (R^2 = %.2f) ', GOF(5)),...
                '\newlineSIR_{Gaus.} / SIR_{Gamm.} =', sprintf(' %.1f', Ratio));

        end
        legendLabels{end+1} = combinedLabel;
    end

    % Plot the Gauss fit if specified
    if biC_flag
        p3Gauss = plot(ax, Combined(:, 1), 100*CombinedComponents(:, 1), '--', 'Color', gaussColor,...
            'LineWidth', 0.5);
        plotHandles(end+1) = p3Gauss;
        if param.latex
            BiGaussLabel = strcat('$\mathrm{(Bi) Gaussian}$');
        else
            BiGaussLabel = sprintf('(Bi) Gaussian');
        end
        legendLabels{end+1} = BiGaussLabel;

        p3Gamma = plot(ax, Combined(:, 1), 100*CombinedComponents(:, 2), '--', 'Color', gammaColor,...
            'LineWidth', 0.5);
        plotHandles(end+1) = p3Gamma;
        if param.latex
            BiGammaLabel = '$\mathrm{(Bi) \Gamma-variate}$';
        else
            BiGammaLabel = '(Bi) Gamma-variate';
        end
        legendLabels{end+1} = BiGammaLabel;

    end


    % Determine the interpreter type based on param.latex
    if param.latex
        interpreter = 'latex';
    else
        interpreter = 'tex';
    end

    % Create the legend
    h_legend = legend(ax, plotHandles, legendLabels, 'Location', ...
        'northeast', 'Interpreter', interpreter);
    set(h_legend, 'FontSize', 8, 'TextColor', 'k', 'Box', 'on', 'Color',...
        'w', 'EdgeColor', 'k');

    % Set axis labels and title
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

    title(ax, ttl, 'Interpreter', interpreter, 'FontSize', 14, 'Color', 'k');

    % Collect y-data from plotted datasets
    yData = 100*measured(:, 2); % Measured data is always plotted
    yData = [yData; 100*GKM(:, 2)];
    yData = [yData; 100*Gamma(:, 2)];
    yData = [yData; 100*GammaF(:, 2)];
    yData = [yData; 100*Gauss(:, 2)];
    yData = [yData; 100*Combined(:, 2)];


    % Calculate the maximum y-value across plotted data sets
    maxYValue = max(yData);

    % Set ylim based on the maxYValue
    if maxYValue < param.SIRLim 
        ylim(ax, [0, param.SIRLim ]);
    else
        ylim(ax, [0, maxYValue*1.1]);
    end

    % Collect x-data from plotted datasets
    xData = measured(:, 1); % Measured data is always plotted
    xData = [xData; GKM(:, 1)];
    xData = [xData; Gamma(:, 1)];
    xData = [xData; GammaF(:, 1)];
    xData = [xData; Gauss(:, 1)];
    xData = [xData; Combined(:, 1)];
    % Adjust x-axis limits
    if isfield(param,'TILim')
        xlim(ax, [0, param.TILim]);
    else
        xlim(ax, [0, max(xData)*1.01]);
    end

    % Set other axis properties
    set(ax, 'TickLabelInterpreter', interpreter);
    set(ax, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'TickDir', 'in',...
        'TickLength', [0.03, 0.015], 'LineWidth', 1);

    % Set grid and minor ticks
    set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', 'MinorGridLineStyle', '-');
    box(ax, 'on');
    axis square
    ax.YTickMode     = 'auto';
    ax.YTickLabelMode= 'auto';

    hold(ax, 'off');

    if isempty(varargin)
        % Save the figure as a PDF using exportgraphics
        exportgraphics(gcf, [ttl,'.pdf'], 'ContentType', 'vector');
    end

    % Close the figure if it's not supposed to be visible
    if strcmp(param.visible, 'off')
        close(fg);
    end

end