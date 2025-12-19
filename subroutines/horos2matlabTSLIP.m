function ROI = horos2matlabTSLIP(csv, fname, isTSLIP,param)
%--------------------------------------------------------------------------
% Combined function to convert Horos CSV data to MATLAB structure and
% display an image with overlaid ROIs, saving as a PDF.
%
% INPUT:  csv    - array of structures with CSV file information
%         fname  - output filename for ROI image
%         isTSLIP - flag to indicate if data is TSLIP (boolean)
%
% OUTPUT: structure with:
%             maskVm       - volumetric mask                   (double)
%             maskIm       - image mask                        (double)
%             slice_number - slice number                      (int)
%             area         - area of ROIs                      (double)
%             name         - names of ROIs                     (cell array)
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

% Read and process CSV files
for r = 1:size(csv, 1)
    filename = csv(r).name;
    if ~startsWith(filename, '.') && ~startsWith(filename, '_')
        fid = fopen(filename);
        FC = textscan(fid, '%s', 'delimiter', '\n');
        fclose(fid);
        FC = FC{1};
        FC = strcat(',', FC);
        data_temp = regexp(FC, ',([^,]*)', 'tokens');
        if str2double(data_temp{2, 1}{1, 1}) == 0
            delta = 1;
        else
            delta = 0;
        end
    end
end

% Create zero volumes
folder = folder_list(pwd);
[Im, S] = dicom2struct_canon(folder(1).name, 'data');
cd ..

%slice order flag
checked = false; 

% Split if series was acquired with alternate % Determine the number of slices
num_slices = numel(unique([S.location]));
num_images = size(Im,3);
slices=1:num_images;

image_positions=zeros(num_images,1);
for sliceN=1:num_images
    image_positions(sliceN)=dicomSliceLocation(S(sliceN).header);
end

[unique_elements, ~] = unique(image_positions, 'stable');
is_descending = issorted(unique_elements, 'descend');

maskV = zeros(size(Im, 1), size(Im, 2), num_slices);
% maskV = zeros(size(Im, 1), size(Im, 2), size(Im, 3) /...
%    (isTSLIP * 2 + ~isTSLIP));

% Initialize variables
maskIm = [];
maskVm = [];
maskIm0 = [];
n = size(data_temp, 1);
ROIname = cell(n - 1, 1);
area = zeros(n - 1, 1);

% Process each entry in CSV data
for entry = 2:n
    data = flattenCell(data_temp{entry, 1});
    temp = data{8};
    ROIname{entry - 1} = erase(temp, '"');
    area(entry - 1) = str2double(data{13});
    X = cat(1, str2double(data(19:5:end)), str2double(data(20:5:end)));
    X = cat(2, X, X(:, 1));
    X(isnan(X)) = [];
    if size(X, 1) == 1
        X = reshape(X, [2, size(X, 2) / 2]);
    end
    Y = interppolygon(X', 100);

    if is_descending && ~checked
            slices=flip(slices,2);
            checked = true;
    end
    checked = true;
    

    slice_number = slices(str2double(data{1, 1})+1);

    if slice_number>num_slices
            slice_number=slice_number-num_slices;
    end


    maskI = double(poly2mask(Y(:, 1), Y(:, 2), size(Im, 1), size(Im, 2)));
    %maskI = imdilate(maskI, strel('diamond', 1));
    maskIm0 = cat(3, maskIm0, maskI);
    maskV(:, :, slice_number) = maskI;
    maskIm = cat(3, maskIm, maskI);
    maskVm = cat(4, maskVm, maskV);
    clear maskI;
end

% Set ROI colors and display image with overlaid ROIs
roiColors = distinguishable_colors(size(maskIm, 3), 'k');

%ROI_imageWtitle(mat2gray(Im(:, :, slice_number)), maskIm0, ...
%    ROIname', slice_number, roiColors, fname);

displayROIImage2(mat2gray(Im(:, :, slice_number)), maskIm0, ...
     ROIname', slice_number, roiColors, fname, param);


% Prepare output structure
maskIm(maskIm ~= 0) = 1;
maskVm = sum(maskVm, 4);
maskVm(maskVm ~= 0) = 1;
maskIm(maskIm == 0) = nan;
maskVm(maskVm == 0) = nan;

ROI.maskVm = maskVm;
ROI.maskIm = maskIm;
ROI.slice_number = slice_number;
ROI.area = area;
ROI.name = ROIname;
end
