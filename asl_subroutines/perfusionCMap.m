function perfusionCMap(B, ROI, TI, slice, param, type)

%--------------------------------------------------------------------------
%
%   Function to display images with B.SIR overlays and TI values
%
%   INPUT:  B       - Structure with fields On, SIR (3D arrays)
%           ROI     - Structure with fields mask, mtop, mbottom, mleft, mright
%                     (integers for bounding box and mask for transparency)
%           TI      - Array of doubles Inversion Times
%           slice   - double
%           param   - Structure with field montageALL (2-element boolean array)
%                     montageALL(2): If true, use the bounding box for cropping
%                     latex: If true, use 'latex' interpreter for labels
%
%   OUTPUT: A figure displaying the images with B.SIR or B.SUBT overlays, 
%           and vertical labels for each row.
%__________________________________________________________________________
% VM vmalis@ucsd.edu
%--------------------------------------------------------------------------

% not passed directly to keep consistency so adjust here if needed
% By default gray scale image is corresponding Tag image for each TI
% with true intensities

underimage = 'On_01';

switch underimage
    case 'On_all'
        Y = mat2gray(B.On);
    case 'On_01'
        Y = mat2gray(B.On(:,:,1));
        Y = repmat(Y,[1,1,numel(TI)]);
    case 'Ctrl_all'
        Y = mat2gray(B.Ctrl);
    case 'Ctrl_01'
        Y = mat2gray(B.Ctrl(:,:,1));
        Y = repmat(Y,[1,1,numel(TI)]);
    case 'M0_all'
        Y = mat2gray(B.M0(:,:,1));
        Y = repmat(Y,[1,1,numel(TI)]);
    case 'M0_01'
        Y = mat2gray(B.M0(:,:,1));
        Y = repmat(Y,[1,1,numel(TI)]);
end

%-------------------------------------------------------------------------

switch type
    case 'SIR' 
        X = 100 * B.SIR;
        cmapUpper = param.SIRClim;
        % Check if the bounding box should be used
        useBoundingBox = param.sir2Dcmap(2);
        if param.latex
            cblable = '$\% \; \mathrm{SIR}$';
        else
            cblable = '% SIR';
        end
    case 'SUBT'
        % Check if subtraction is normalized
        X = B.SUBT;
        if max(X, [], "all") > 1
            X = 100*mat2gray(X);
        else
            X = 100*X;
        end
        if param.latex
            cblable = '$\% \; |\mathrm{On} - \mathrm{CTR}|$';
        else
            cblable = '% | On - CTR |';
        end

        cmapUpper = param.SIRClim;
        % Check if the bounding box should be used
        useBoundingBox = param.subtractionMap(2);
end

if param.sir2DcmapFull
    mask = ROI.mask;
else
    mask = ROI.maskF;
end


% Function to get the image or its bounding box
getImage = @(img, idx) img(:, :, idx);
if useBoundingBox
    getImage = @(img, idx) img(ROI.mtop:ROI.mbottom, ROI.mleft:ROI.mright, idx);
end

% Resize image to be at least 300 pixels in width and height
resizeImage = @(img) imresize(img, max(300 ./ size(img, 1:2)), 'bicubic');

% Number of images
numImages = size(Y, 3);

% Process images
Y_resized = arrayfun(@(i) resizeImage(getImage(Y, i)), 1:numImages, 'UniformOutput', false);
X_resized = arrayfun(@(i) resizeImage(getImage(X, i)), 1:numImages, 'UniformOutput', false);
mask_resized = resizeImage(getImage(mask, 1));

% Get the size of a single image
[imgHeight, imgWidth, ~] = size(Y_resized{1});

% Calculate figure size based on the number of images and their dimensions
figWidth = imgWidth * (numImages + 1);
figHeight = imgHeight * 1.2;  % 1 row of images with space for labels

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
        [0, 0, figWidth, figHeight],  'Visible', 'off');
set(fig, 'Visible', param.visible);

% Create a tiled layout for the images
t = tiledlayout(1, numImages + 1, 'TileSpacing', 'none', 'Padding', 'none');

% Initialize arrays to store axis handles
ax1Array = gobjects(1, numImages);
ax2Array = gobjects(1, numImages);


FontSize=ceil(figHeight/20);

if FontSize<3
    FontSize=3;
elseif FontSize>16
    FontSize=16;
end



% Plot ON images with SIR overlays
for col = 1:numImages

    % Create a new axis for the underlying grayscale image
    ax1Array(col) = nexttile(t, col);
    imgOn = Y_resized{col};
    imshow(imgOn, 'Parent', ax1Array(col));
    axis tight;
    % Turn off ticks
    set(gca, 'XTick', [], 'YTick', []);
    colormap(ax1Array(col), gray);
    ax1Array(col).CLim = [0 1];
    axis(ax1Array(col), 'off');
    
    % Pause to ensure position is updated correctly
    drawnow;
    
    % Create a new axis for the overlay image at the same position
    ax2Array(col) = axes('Position', get(ax1Array(col), 'Position'), 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
    imgSIR = X_resized{col};
    mask = mask_resized;
    
    % Calculate alpha data for gradual transparency
    alphaData = imgSIR / cmapUpper;
    alphaData(alphaData > 1) = 1;  % Cap the values at 1
    alphaData(~mask) = 0;  % Set alpha to 0 for masked areas
    
    h = imagesc(imgSIR, 'Parent', ax2Array(col));
    set(h, 'AlphaData', alphaData);  % Apply calculated alpha data
    colormap(ax2Array(col), blackToRed(1000));
    clim(ax2Array(col), [0 cmapUpper]);
    axis(ax2Array(col), 'off');


    % Add TI label as x-label
    if param.latex
        xlabel(ax1Array(col), ['$' num2str(TI(col)) '\, \mathrm{ms}$'], 'Color', 'w', 'FontSize', FontSize, 'HorizontalAlignment', 'center','Interpreter','latex');
    else
        xlabel(ax1Array(col), [num2str(TI(col)) ' ms'], 'Color', 'w', 'FontSize', FontSize, 'HorizontalAlignment', 'center');    
    end
end

% Link the axes
linkaxes([ax1Array, ax2Array]);

cbAx = nexttile(t, numImages + 1);
axis tight;
 % Turn off ticks
set(gca, 'XTick', [], 'YTick', []);
% Create an overall colorbar outside the tiled layout
h = imagesc(imgSIR*0);
colormap(gca, blackToRed(1000));
clim(gca, [0 cmapUpper]);
axis(gca, 'off');
axis tight;
 % Turn off ticks
set(gca, 'XTick', [], 'YTick', []);


cb = colorbar;
cb.Color = 'w';

if param.latex
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
end
cb.Label.String = cblable;
% Use the default text interpreter

% Match the color limits of the overlay images
cb.Limits = [0 cmapUpper]; 
% Ensure the colorbar uses the same colormap as the overlay images
colormap(cb, blackToRed(1000)); 
% Place tick labels on the right side
cb.TickDirection = 'out';
% Ensure labels and ticks are on the right side
% interpreter
cb.FontSize = FontSize;
cb.Location = "east";

% Get the current position of the colorbar
cb_pos = cb.Position; % [x, y, width, height]
cb_pos(4)=0.8;
cbHeight = figHeight*cb_pos(4); % Desired new width
cb_pos(3) = cbHeight/8/figWidth; % Set the third element (width)
cb_pos(2)=(1-cb_pos(4))/2;
    if col<10
        cb_pos(1)=cb_pos(1)-1/col/2;
    elseif col>=10 && col<25
        cb_pos(1)=cb_pos(1)-1/col/5;
    end
cb.Position = cb_pos;
cb.YAxisLocation = 'right'; 

% Ensure all positions are updated correctly
drawnow;

% Force the position of ax2Array to match ax1Array
for col = 1:numImages
    ax2Array(col).Position = ax1Array(col).Position;
end

linkprop([ax2Array cbAx],'CLim');

% Set figure background to black
set(gcf, 'color', 'k');

% Save the figure as a PDF file

export_fig([type, '_s', sprintf('%02d', slice)], '-pdf', '-opengl','-silent','-nocrop');


% Close the figure
close;
close all;

end

function customColormap = blackToRed(n)
    if nargin < 1
        n = 256; % Default number of colors
    end
    r = linspace(0, 1, n); % Red channel goes from 0 to 1
    g = zeros(1, n);       % Green channel remains 0
    b = zeros(1, n);       % Blue channel remains 0
    customColormap = [r' g' b'];
end