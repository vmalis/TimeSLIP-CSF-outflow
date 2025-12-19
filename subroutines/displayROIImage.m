function displayROIImage(I, mask, masknames, fileNum, colors, ttl, param, varargin)
%--------------------------------------------------------------------------
% Subroutine to display an image with overlaid ROIs and save as a PDF
% varargin can be used to pass axis handle if available
%
% INPUT:
%    I         - Image data to display
%    mask      - 3D array of binary masks, each slice representing an ROI
%    masknames - Cell array of names for each ROI
%    fileNum   - File number for naming the output PDF
%    colors    - Array of colors for each ROI
%    ttl       - Title for the image
%    param     - Structure with display parameters (e.g., visibility, latex)
%    varargin  - Optional axis handle for plotting
%
% OUTPUT: 
%    None (the function saves the image with overlaid ROIs as a PDF file)
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------
% Check if an axis handle is provided
if isempty(varargin)
    fig_roi = figure('Position', [0, 0, 600, 600], 'Color', 'k', 'Visible', 'on');
    set(fig_roi, 'Visible', param.visible);

    set(fig_roi, 'Color', 'k');
    
    % Display image
    imshow(I, 'InitialMagnification', 'fit');
    
    % Overlay different masks
    l = 0;  % Plot number, needed for proper legends
    m = []; % Plots that should be in legend
    for j = 1:size(mask, 3)
        alphamask(mask(:, :, j), colors(j, :), 0.05, gca);
        
        [B, ~] = bwboundaries(mask(:, :, j));
        hold on
        for k = 1:length(B)
            l = l + 1;
            boundary = B{k};
            h(l) = plot(gca, boundary(:, 2), boundary(:, 1), 'Color', ...
                colors(j, :), 'LineWidth', 1);
            if k == 1
                m = [m, l];
            end
        end
    end

     % Determine the interpreter type based on param.latex
    if isfield(param, 'latex') && param.latex
        interpreter = 'latex';
    else
        interpreter = 'none';
    end
    
    % Add legend
    lgd = legend(h(m), masknames, 'Location', 'eastoutside', 'Color', ...
        'none', 'FontSize', 14,'Interpreter',interpreter);
    lgd.TextColor = 'w';
    lgd.EdgeColor = 'w';
    title(ttl, 'Color', 'w','Interpreter',interpreter);

    % Save image with ROIs as PDF
    slice_id = strcat('s', sprintf('%02d', fileNum));
    export_fig(['ROIs_', slice_id, '.pdf'], '-nocrop', ...
        '-opengl','-silent');
    close all
    delete(fig_roi);



else
    ax = varargin{1};

    % Display the image in the specified axis
    imshow(I, 'InitialMagnification','fit', 'Parent', ax);
    axis(ax, 'tight');
    hold(ax, 'on'); % Ensure hold is on
    
    l = 0;  % Plot number, needed for proper legends
    m = []; % Plots that should be in legend
    for j = 1:size(mask, 3)
        alphamask(mask(:, :, j), colors(j, :), 0.05, ax);
        [B, ~] = bwboundaries(mask(:, :, j));
        for k = 1:length(B)
            l = l + 1;
            boundary = B{k};
            hold(ax, 'on'); % Ensure hold is on
            h(l) = plot(ax, boundary(:, 2), boundary(:, 1), 'Color', colors(j, :), 'LineWidth', 0.5);
            if k == 1
                m = [m, l];
            end
        end
    end
    
    % Determine the interpreter type based on param.latex
    if isfield(param, 'latex') && param.latex
        interpreter = 'latex';
    else
        interpreter = 'none';
    end
    
    
    % Add legend and title with the specified interpreter
    lgd = legend(ax, h(m), masknames, 'Location', 'southeast', 'Color','k',...
        'EdgeColor', 'w','FontSize', 12, 'Interpreter', interpreter);
    lgd.TextColor = 'w';
    lgd.EdgeColor = 'none';

end


end