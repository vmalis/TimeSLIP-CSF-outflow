function tslipmontage(B, ROI, TI, slice, param)

%--------------------------------------------------------------------------
%
%   Function to display ON, Ctrl, and M0 images with TI values
%
%   INPUT:  B       - Structure with fields M0, On, Ctrl (3D arrays)
%           ROI     - Structure with fields mtop, mbottom, mleft, mright
%                     (integers for bounding box)
%           TI      - Array of doubles, same size as the 3rd dimension of B.On
%           param   - Structure with field montageALL (2-element boolean array)
%                     montageALL(1): If true, show all images in one montage
%                     montageALL(2): If true, use the bounding box for cropping
%                     contrastNorm: If true, normalize the image
%
%   OUTPUT: A figure displaying the ON, Ctrl, and M0 images with TI values
%           and vertical labels for each row.
%__________________________________________________________________________
% VM vmalis@ucsd.edu
%--------------------------------------------------------------------------


% Step 1: Find the global minimum and maximum
allData = [B.On(:); B.Ctrl(:); B.M0(:)];
globalMin = min(allData);
globalMax = max(allData);

% Step 2: Normalize each array using the global minimum and maximum
B.On = (B.On - globalMin) / (globalMax - globalMin);
B.Ctrl = (B.Ctrl - globalMin) / (globalMax - globalMin);
B.M0 = (B.M0 - globalMin) / (globalMax - globalMin);

% Check if M0 is 3D or 2D and normalize
if ismatrix(B.M0)
    B.M0 = repmat(B.M0, 1, 1, size(B.On, 3)); % Duplicate M0 image for all TIs
end

% Determine if bounding box should be used
useBoundingBox = param.montageALL(2);
if useBoundingBox
    % Extract the bounding box
    mtop = ROI.mtop;
    mbottom = ROI.mbottom;
    mleft = ROI.mleft;
    mright = ROI.mright;
end

% Function to get the image or its bounding box and normalize if needed
if useBoundingBox && param.contrastNorm
    getImage = @(img, idx) mat2gray(img(mtop:mbottom, mleft:mright, idx));
elseif useBoundingBox && ~param.contrastNorm
    getImage = @(img, idx) img(mtop:mbottom, mleft:mright, idx);
elseif ~useBoundingBox && param.contrastNorm
    getImage = @(img, idx) mat2gray(img(:, :, idx));
else
    getImage = @(img, idx) img(:, :, idx);
end

% Resize image to be at least 300 pixels in width and height
resizeImage = @(img) imresize(img, max(300 ./ size(img, 1:2)), 'bicubic');

% Number of images
numImages = size(B.On, 3);

% Process images
B.On_resized = arrayfun(@(i) resizeImage(getImage(B.On, i)), 1:numImages, 'UniformOutput', false);
B.Ctrl_resized = arrayfun(@(i) resizeImage(getImage(B.Ctrl, i)), 1:numImages, 'UniformOutput', false);
B.M0_resized = arrayfun(@(i) resizeImage(getImage(B.M0, i)), 1:numImages, 'UniformOutput', false);

% Get the size of a single image
[imgHeight, imgWidth, ~] = size(B.On_resized{1});

% Calculate figure size based on the number of images and their dimensions
figWidth = imgWidth * numImages;
figHeight = imgHeight * 3;  % 3 rows of images

% Get screen size
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% Adjust figure size to fit within screen size
if figWidth > screenWidth || figHeight > screenHeight
    scaleFactor = min(screenWidth / figWidth, screenHeight / figHeight);
    figWidth = figWidth * scaleFactor;
    figHeight = figHeight * scaleFactor;
end

% Create figure with calculated size
fig = figure('Color', 'k', 'InvertHardcopy', 'off', 'Position', ...
    [0, 0, figWidth, figHeight],'Visible','off');
set(fig, 'Visible', param.visible);

% Create a tiled layout for the images with 3 rows
t = tiledlayout(3, numImages, 'TileSpacing', 'none', 'Padding', 'none');

% Determine if 'latex' interpreter should be used
interpreterType = 'none';
if isfield(param, 'latex') && param.latex
    interpreterType = 'latex';
end

FontSize=ceil(figHeight/20);

if FontSize<3
    FontSize=3;
elseif FontSize>16
    FontSize=16;
end

% Plot ON, Ctrl, and M0 images
for row = 1:3
    for col = 1:numImages
        ax = nexttile(t, col + (row-1) * numImages);
        switch row
            case 1
                imshow(B.On_resized{col});
                axis tight;
                % Turn off ticks
                set(gca, 'XTick', [], 'YTick', []);

                if col == 1
                    ylabel(ax, 'On', 'Color', 'w', 'FontSize', FontSize, 'Rotation', 90, 'HorizontalAlignment', 'center', 'Interpreter', interpreterType);
                end
            case 2
                imshow(B.Ctrl_resized{col});
                axis tight;
                % Turn off ticks
                set(gca, 'XTick', [], 'YTick', []);

                if col == 1
                    ylabel(ax, 'Control', 'Color', 'w', 'FontSize', FontSize, 'Rotation', 90, 'HorizontalAlignment', 'center', 'Interpreter', interpreterType);
                end
            case 3
                imshow(B.M0_resized{col});
                axis tight;
                % Turn off ticks
                set(gca, 'XTick', [], 'YTick', []);

                if col == 1
                    if strcmp(interpreterType, 'latex')
                        ylabel(ax, '$M_0$', 'Color', 'w', 'FontSize', FontSize, 'Rotation', 90, 'HorizontalAlignment', 'center', 'Interpreter', interpreterType);
                    else
                        ylabel(ax, 'M_0', 'Color', 'w', 'FontSize', FontSize, 'Rotation', 90, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
                    end
                end

                if strcmp(interpreterType, 'latex')
                    xlabel(ax, ['$' num2str(TI(col)) '\, \mathrm{ms}$'], 'Color', 'w', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'Interpreter', interpreterType);
                else
                    xlabel(ax, [num2str(TI(col)) ' ms'], 'Color', 'w', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
                end
        end
        axis off;
    end
end

% Set figure background to black and resize
set(fig, 'color', 'k');

% Save the figure as a PDF file
export_fig(['montage_', 's', sprintf('%02d', slice)], '-pdf', '-opengl',...
    '-silent','-nocrop');

% Close the figure
close;
close all;


end