function location = dicomSliceLocation(dicomHeader)

%--------------------------------------------------------------------------
% dicomSliceLocation - Computes the slice location from a DICOM header.
%
% This function calculates the slice location based on the Image Position
% (Patient) and Image Orientation (Patient) fields from the DICOM header.
%
% INPUT:  
%   dicomHeader - Structure containing DICOM metadata (dicominfo output).
%
% OUTPUT: 
%   location    - Computed slice location (double).
%
% USAGE:
%   location = dicomSliceLocation(dicomHeader);
%
% NOTES:
%   - If the DICOM header contains the Slice Location (0020,1041) field,
%     this function will return that value.
%   - Otherwise, the slice location is computed using the normal vector of
%     the image plane and the image position.
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------


    % Extract Image Position (Patient) [X0, Y0, Z0]
    if isfield(dicomHeader, 'ImagePositionPatient')
        imagePosition = dicomHeader.ImagePositionPatient;
    else
        error('ImagePositionPatient field is missing in DICOM header.');
    end

    % Extract Image Orientation (Patient) [Xr, Yr, Zr, Xc, Yc, Zc]
    if isfield(dicomHeader, 'ImageOrientationPatient')
        imageOrientation = dicomHeader.ImageOrientationPatient;
    else
        error('ImageOrientationPatient field is missing in DICOM header.');
    end


    imageOrientation = dicomHeader.ImageOrientationPatient;

    % Extract row and column direction cosines
    Xr = imageOrientation(1);
    Yr = imageOrientation(2);
    Zr = imageOrientation(3);
    Xc = imageOrientation(4);
    Yc = imageOrientation(5);
    Zc = imageOrientation(6);

    % Compute the slice normal vector using the cross-product
    N = cross([Xr, Yr, Zr], [Xc, Yc, Zc]); % [Nx, Ny, Nz]

    % Compute the location using the dot product
    location = dot(N, imagePosition);
end