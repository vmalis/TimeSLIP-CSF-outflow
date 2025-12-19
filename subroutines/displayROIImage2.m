function displayROIImage2(I, mask, masknames, fileNum, colors, ttl, param, varargin)
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

if isempty(varargin)
    fig_roi = figure('Position', [0, 0, 600, 600], 'Color', 'k', 'Visible', 'on');
    set(fig_roi, 'Visible', param.visible);
    set(fig_roi, 'Color', 'k');
    
    %---------------------- ORIGINAL MAIN IMAGE CODE ----------------------%
    imshow(I, 'InitialMagnification', 'fit');
    mainAx = gca;
    
    l = 0;
    m = [];
    for j = 1:size(mask, 3)
        alphamask(mask(:, :, j), colors(j, :), 0.05, mainAx);
        [B, ~] = bwboundaries(mask(:, :, j));
        hold(mainAx, 'on')
        for k = 1:length(B)
            l = l + 1;
            boundary = B{k};
            h(l) = plot(mainAx, boundary(:, 2), boundary(:, 1), 'Color', ...
                colors(j, :), 'LineWidth', 1);
            if k == 1
                m = [m, l];
            end
        end
    end

    if isfield(param, 'latex') && param.latex
        interpreter = 'latex';
    else
        interpreter = 'none';
    end

    lgd = legend(h(m), masknames, 'Location', 'southeastoutside', 'Color', ...
        'none', 'FontSize', 12,'Interpreter',interpreter);
    lgd.TextColor = 'w';
    lgd.EdgeColor = 'w';
    title(ttl, 'Color', 'w','Interpreter',interpreter);
    %-------------------- END ORIGINAL MAIN IMAGE CODE --------------------%

    %================== NEW: ZOOMED-IN ROI INSET + RECT ===================%
    combinedMask = any(mask, 3);
    [rowIdx, colIdx] = find(combinedMask);

    if ~isempty(rowIdx)
        pad  = 10;
        imgH = size(I, 1);
        imgW = size(I, 2);

        rmin = max(min(rowIdx) - pad, 1);
        rmax = min(max(rowIdx) + pad, imgH);
        cmin = max(min(colIdx) - pad, 1);
        cmax = min(max(colIdx) + pad, imgW);

        zoomRect = [cmin, rmin, (cmax - cmin), (rmax - rmin)];

        %%%% EDIT #1 — Make main yellow box 3× thinner
        rectangle(mainAx, 'Position', zoomRect, ...
                  'EdgeColor', [1 1 0], 'LineWidth', 0.5);   % was 1.5

        % Crop region
        if ndims(I) == 2
            I_zoom = I(rmin:rmax, cmin:cmax);
        else
            I_zoom = I(rmin:rmax, cmin:cmax, :);
        end
        mask_zoom = mask(rmin:rmax, cmin:cmax, :);

        insetWidth  = 0.25;
        insetHeight = 0.25;

        % Move inset to TOP-RIGHT corner
        insetX = 1 - insetWidth - 0.07;   % right margin
        insetY = 1 - insetHeight - 0.07;  % top margin

        insetAx = axes('Parent', fig_roi, ...
                       'Position', [insetX, insetY, insetWidth, insetHeight]);

        % Draw zoomed region
        imshow(I_zoom, 'InitialMagnification', 'fit', 'Parent', insetAx);
        hold(insetAx, 'on');

        l2 = 0;
        for j = 1:size(mask_zoom, 3)
            alphamask(mask_zoom(:, :, j), colors(j, :), 0.05, insetAx);
            [Bz, ~] = bwboundaries(mask_zoom(:, :, j));
            hold(insetAx, 'on')
            for k = 1:length(Bz)
                l2 = l2 + 1;
                boundary = Bz{k};
                plot(insetAx, boundary(:, 2), boundary(:, 1), ...
                     'Color', colors(j, :), 'LineWidth', 1);
            end
        end

        axis(insetAx, 'image');
        set(insetAx, 'YDir', 'reverse');
        axis(insetAx, 'off');

        %%%% EDIT #2 — Add same bounding box to inset
        % zoomRect in cropped coordinates becomes [1,1,width,height]
        zoomRectInset = [1, 1, (cmax - cmin), (rmax - rmin)];
        rectangle(insetAx, 'Position', zoomRectInset, ...
                  'EdgeColor', [1 1 0], 'LineWidth', 0.5);
    end
    %================ END NEW ZOOMED-IN ROI INSET + RECT ==================%

    slice_id = strcat('s', sprintf('%02d', fileNum));
    exportgraphics(gcf,['ROIs_', slice_id, '.pdf'])
    %export_fig(['ROIs_', slice_id, '.pdf'], '-nocrop', ...
    %    '-opengl','-silent');
    close all
    delete(fig_roi);



else
    % ===== EXTERNAL AXIS BRANCH (UNCHANGED) =====
    ax = varargin{1};
    imshow(I, 'InitialMagnification','fit', 'Parent', ax);
    axis(ax, 'tight');
    hold(ax, 'on');
    
    l = 0;
    m = [];
    for j = 1:size(mask, 3)
        alphamask(mask(:, :, j), colors(j, :), 0.05, ax);
        [B, ~] = bwboundaries(mask(:, :, j));
        for k = 1:length(B)
            l = l + 1;
            boundary = B{k};
            h(l) = plot(ax, boundary(:, 2), boundary(:, 1), 'Color', colors(j, :), 'LineWidth', 0.5);
            if k == 1
                m = [m, l];
            end
        end
    end
    
    if isfield(param, 'latex') && param.latex
        interpreter = 'latex';
    else
        interpreter = 'none';
    end

    lgd = legend(ax, h(m), masknames, 'Location', 'southeast', ...
        'Color','k','EdgeColor', 'w','FontSize', 12, 'Interpreter', interpreter);
    lgd.TextColor = 'w';
    lgd.EdgeColor = 'none';

end

end