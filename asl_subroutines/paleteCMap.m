function paleteCMap(A, ROI, TI, slice, param)
%--------------------------------------------------------------------------
%
%   Subroutine to create and save palette maps from ON and OFF images
%
%   INPUT:  A         Struct containing:
%                      - On:     3D array of ON images          ('double')
%                      - Ctrl:   3D array of Control images     ('double')
%                      - M0:     3D array of Mo images          ('double')
%                      - SUBT:   3D array of Subtraction images ('double')
%                      - SIR:    3D array of SIR images         ('double')
%           ROI       Struct containing:
%                      - mask: 2D binary mask for region of interest
%                      - mtop, mbottom, mleft, mright: Indices for cropping
%           TI        Array of inversion times ('double')
%           slice     Integer indicating the slice number ('double')
%           param     Struct containing:
%                      - latex: Boolean indicating whether to use LaTeX
%                      - paleteviewZoom: Boolean for zooming into ROI
%                      - acq_dim: Acquisition dimension ('3d' or '2d')
%                      - filter_median2D: Boolean for applying median filter
%
%   OUTPUT: None (saves palette maps as PDF files)
%__________________________________________________________________________
% VM vmalis@ucsd.edu
%--------------------------------------------------------------------------


% Load the colormap
load('MRICmap2.mat');
latex=param.latex;

% Get screen size
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);


% apply mask
mask = ROI.mask;

% Step 1: Find the global minimum and maximum
allData = [A.On(:); A.Ctrl(:); A.M0(:)];
globalMin = min(allData);
globalMax = max(allData);

% Step 2: Normalize each array using the global minimum and maximum, then apply the mask
ON = ((A.On - globalMin) / (globalMax - globalMin)) .* mask;
CTRL = ((A.Ctrl - globalMin) / (globalMax - globalMin)) .* mask;
M0 = ((A.M0 - globalMin) / (globalMax - globalMin)) .* mask;

SUBT = A.SUBT;
SIR  = 100 * (A.SIR).* mask;

%check if subtraction is normalized
if max(SUBT,[],'all','omitnan')>1
    SUBT=mat2gray(SUBT).* mask;
end

if param.noiseM0
    SUBT=SUBT.*ROI.Noise;
    SIR=SIR.*ROI.Noise;
end


% Loop through each slice in the 3D array
for i = 1:size(ON, 3)

    %crop and zoom to ROI
    if param.paleteview(2)
        Ctrl = CTRL(ROI.mtop:ROI.mbottom, ROI.mleft:ROI.mright, i);
        On   = ON(ROI.mtop:ROI.mbottom, ROI.mleft:ROI.mright, i);
        if max(size(size(M0))) > 2
            m0 = M0(ROI.mtop:ROI.mbottom, ROI.mleft:ROI.mright, i);
        else
            m0 = M0(ROI.mtop:ROI.mbottom, ROI.mleft:ROI.mright);
        end
        Subt = SUBT(ROI.mtop:ROI.mbottom, ROI.mleft:ROI.mright, i);
        Sir  = SIR(ROI.mtop:ROI.mbottom, ROI.mleft:ROI.mright, i);
    else
        Ctrl = CTRL(:, :, i);
        On   = ON(:, :, i);
        if max(size(size(M0))) > 2
            m0 = M0(:, :, i);
        else
            m0 = M0;
        end
        Subt = SUBT(:, :, i);
        Sir  = SIR(:, :, i);
    end


     if param.contrastNorm
            Ctrl = mat2gray(Ctrl);
            On   = mat2gray(On);
            m0   = mat2gray(m0);
     end

    

    % Concatenate processed images into a 3D array
    total = cat(3, On, Ctrl, m0, Subt, Sir);

    % Create a tiled layout for the images
    fig = figure('Visible', 'off');
    t = tiledlayout(1, 5);
    
    % Set figure background to black and resize
    set(fig, 'color', 'k');
    set(fig, 'Position', [0, 0, screenWidth, screenWidth/3]);
    set(fig, 'Visible', param.visible);

    % Display the ON image
    clims = [0 1];
    ax1 = nexttile;
    imagesc(total(:, :, 1), clims);
    colormap(ax1, gray);
    axis image;
    axis off;
    if latex
        text(size(total, 2) / 2, -10, '$\mathrm{On}$', 'Interpreter',...
            'latex', 'FontSize', 14, 'Color', 'w', ...
            'HorizontalAlignment','center');
    else
        text(size(total, 2) / 2, -10, 'On', 'FontSize', 14, 'Color',...
            'w', 'HorizontalAlignment', 'center');
    end

    % Display the CTRL image
    clims = [0 1];
    ax2 = nexttile;
    imagesc(total(:, :, 2), clims);
    colormap(ax2, gray);
    axis image;
    axis off;
    if latex
        text(size(total, 2) / 2, -10, '$\mathrm{Control}$', ...
            'Interpreter', 'latex', 'FontSize', 14, 'Color', 'w', ...
            'HorizontalAlignment', 'center');
    else
        text(size(total, 2) / 2, -10, 'Control', 'FontSize', 14, ...
            'Color', 'w', 'HorizontalAlignment', 'center');
    end

    % Display the Subt image with the custom colormap
    clims = [0 100];
    ax4 = nexttile;
    imagesc(100*total(:, :, 4), clims);
    colormap(ax4, uint8(MRICmap2));
    axis image;
    axis off;
    if latex
        text(size(total, 2) / 2, -10,...
            '$|\mathrm{On}-\mathrm{Control}|$',...
            'Interpreter', 'latex', 'FontSize', 14, 'Color', 'w', ...
            'HorizontalAlignment', 'center');
    else
        text(size(total, 2) / 2, -10, ...
           sprintf('| On - Control |'),...
           'FontSize', 14,'Color', 'w', 'HorizontalAlignment', 'center');
    end

    % Add a colorbar to the last image
    yc = colorbar;
    yc.Color = 'w';
    yc.FontSize = 14;
    yc.Location = 'southoutside';
    if latex
        yc.TickLabelInterpreter = 'latex';
        ylabel(yc, '$\% \, (normalized)$', 'FontSize', 14, 'color', 'w', 'Interpreter',...
            'latex');
    else
        ylabel(yc, '% (normalized)', 'FontSize', 14, 'color', 'w');
    end

    % Display the M0 image
    clims = [0 1];
    ax3 = nexttile;
    imagesc(total(:, :, 3), clims);
    colormap(ax3, gray);
    axis image;
    axis off;
    if latex
        text(size(total, 2) / 2, -10, '$\mathrm{M}_0$', 'Interpreter',...
            'latex', 'FontSize', 14, 'Color', 'w', ...
            'HorizontalAlignment', 'center');
    else
        text(size(total, 2) / 2, -10, 'M_0', 'FontSize', 14, ...
            'Color', 'w', 'HorizontalAlignment', 'center');
    end

    

    % Display the SIR image with the custom colormap
    ax5 = nexttile;
    clims = [0 param.SIRClim];
    imagesc(total(:, :, 5), clims);
    colormap(ax5, uint8(MRICmap2));
    axis image;
    axis off;
    if latex
        text(size(total, 2) / 2, -10, '$\mathrm{Signal \, Increase\, Ratio}$', 'Interpreter',...
            'latex', 'FontSize', 14, 'Color', 'w', ...
            'HorizontalAlignment', 'center');
    else
        text(size(total, 2) / 2, -10, 'Signal Increase Ratio', 'FontSize', 14, 'Color',...
            'w', 'HorizontalAlignment', 'center');
    end

    % Add a colorbar to the last image
    xc = colorbar;
    xc.Color = 'w';
    xc.FontSize = 14;
    xc.Location = 'southoutside';
    if latex
        xc.TickLabelInterpreter = 'latex';
        ylabel(xc, '$\% \, \mathrm{SIR}$', 'FontSize', 14, 'color', 'w', 'Interpreter',...
            'latex');
    else
        ylabel(xc, '% SIR', 'FontSize', 14, 'color', 'w');
    end

    % Adjust tile spacing
    t.TileSpacing = 'none';
    t.Padding = 'loose';

    if param.fixedRO
        % Save the figure as a PDF file
        exportgraphics(gcf,['palette_', 'ECG', sprintf('%04d', TI(i)), 's',...
            sprintf('%02d', slice),'.pdf']);
    else
        % Save the figure as a PDF file
        exportgraphics(gcf,['palette_', 'TI', sprintf('%04d', TI(i)), 's',...
            sprintf('%02d', slice),'.pdf']);
    end

    % Close the figure
    close;
    close all;
end


end
