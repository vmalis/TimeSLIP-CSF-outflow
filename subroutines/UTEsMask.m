function [Mask] = UTEsMask(Images,intensity_threshold)

%==========================================================================
% Subroutine to creat full axial slice mask
%==========================================================================
%
% INput:  Images                  image volume [size = xres, yres, nslices]
%         intensity_threshold     intensity threshold
%
% OUTput: Mask [size = xres, yres, nslices]
%
%--------------------------------------------------------------------------
% written by Vadim Malis (vmalis@ucsd.edu)
% 04/20 at UCSD
%==========================================================================
thrsh = intensity_threshold;
I = Images;

N = size(I, 3); % Number of slices
Mask = zeros(size(I, 1), size(I, 2), N);

for n = 1:N
    bw = zeros(size(I, 1), size(I, 2));
    bw(I(:, :, n, 1) > thrsh) = 1; % Thresholding
    bw = logical(bw);

    % Extract up to two largest objects
    temp_bw = bwareafilt(bw, 2); % Keep up to 2 largest connected components

    % Calculate areas of the two largest objects
    stats = regionprops(temp_bw, 'Area');
    if numel(stats) == 2
        % Check the size difference between the two objects
        area1 = stats(1).Area;
        area2 = stats(2).Area;
        if abs(area1 - area2) / max(area1, area2) > 0.2
            % If difference > 20%, keep only the largest object
            temp_bw = bwareafilt(temp_bw, 1);
        end
    end

    % Fill holes in the retained objects
    temp_bw = imfill(temp_bw, 'holes');

    % Blur by convolution to smooth edges
    windowSize = 20;
    kernel = ones(windowSize) / windowSize^2;
    blurryImage = conv2(single(temp_bw), kernel, 'same');

    % Store in Mask with rethresholding
    Mask(:, :, n) = blurryImage > 0.6; % Rethreshold
end
end