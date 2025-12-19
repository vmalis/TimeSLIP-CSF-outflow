function interpolatedImage = cropAndResizeImage(Image, ROI, roiNumber)
%--------------------------------------------------------------------------
%
%   Function to crop an image based on ROI and resize it to 300x300 pixels
%
%   INPUT:
%       Image - the input image
%       ROI - structure with fields:
%           mtop - top row for the ROI
%           mbottom - bottom row for the ROI
%           mleft - left column for the ROI
%           mright - right column for the ROI
%
%   OUTPUT:
%       interpolatedImage - the cropped and resized image (300x300 pixels)
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

% ----- Default output size -----
if nargin < 3 
    roiNumber = 1;
end

roiLabel = upper(strtrim(char(ROI.name(roiNumber))));

switch roiLabel
    case {'WM','GM','CSF'}
        [bottom, right, ~] = size(Image);
        top = 0;
        left = 0;
    otherwise
        top = ROI.mtop;
        bottom = ROI.mbottom;
        left = ROI.mleft;
        right = ROI.mright;
end

% Get the size of the image
[imageHeight, imageWidth, ~] = size(Image);

% Calculate the width and height of the ROI
roiWidth = right - left + 1;
roiHeight = bottom - top + 1;

% Determine the side length of the square crop
sideLength = max(roiWidth, roiHeight);

% Adjust the crop to be square
if roiWidth < roiHeight
    % Expand width
    extraWidth = sideLength - roiWidth;
    left = left - floor(extraWidth / 2);
    right = right + ceil(extraWidth / 2);
else
    % Expand height
    extraHeight = sideLength - roiHeight;
    top = top - floor(extraHeight / 2);
    bottom = bottom + ceil(extraHeight / 2);
end

% Ensure the crop is within the image bounds
if left < 1
    left = 1;
    right = left + sideLength - 1;
end
if right > imageWidth
    right = imageWidth;
    left = right - sideLength + 1;
end
if top < 1
    top = 1;
    bottom = top + sideLength - 1;
end
if bottom > imageHeight
    bottom = imageHeight;
    top = bottom - sideLength + 1;
end

% Adjust if still out of bounds
left = max(1, left);
right = min(imageWidth, right);
top = max(1, top);
bottom = min(imageHeight, bottom);

% Perform the cropping
croppedImage = Image(top:bottom, left:right, :);

% Check if the image contains 1s and NaNs
if all(isnan(croppedImage(:)) | croppedImage(:) == 1)
    % Replace NaNs with zeros
    croppedImage(isnan(croppedImage)) = 0;
    
    % Use nearest-neighbor interpolation
    interpolatedImage = imresize(croppedImage, [300 300], 'nearest');
    
    % Restore NaN values
    interpolatedImage(interpolatedImage == 0) = NaN;
else
    % Use default interpolation for other images
    interpolatedImage = imresize(croppedImage, [300 300]);
end

end