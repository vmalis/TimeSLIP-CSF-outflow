function filteredA = arraySDfilter(A, M)
%--------------------------------------------------------------------------
%
%   Subroutine to filter out outliers in a 2D or 3D array based on 
%   standard deviations
%
%   INPUT:  2D or 3D array of doubles and NaNs                    'double'
%           Number of standard deviations for filtering (M)       'double'
%   OUTPUT: Filtered 2D or 3D array                               'double'
%__________________________________________________________________________
% VM vmalis@ucsd.edu
%--------------------------------------------------------------------------

    % Initialize the output array
    filteredA = A;
    
    % Check dimensions and process accordingly
    if ismatrix(A)
        A = reshape(A, size(A, 1), size(A, 2), 1);
    end
    
    for c = 1:size(A, 3)
        slice = A(:, :, c);
        non_nan_values = slice(~isnan(slice));
        mean_value = mean(non_nan_values);
        stdev_value = std(non_nan_values);
        
        % Create a logical mask for values within M standard deviations
        mask = abs(slice - mean_value) <= M * stdev_value;
        
        % Set values outside the range to NaN
        slice(~mask) = NaN;
        
        % Update the output array
        filteredA(:, :, c) = slice;
    end
    
    % Reshape back to 2D if necessary
    if size(filteredA, 3) == 1
        filteredA = reshape(filteredA, size(filteredA, 1), size(filteredA, 2));
    end

end