function metadata = dicomCanonScaling(metadata)
%--------------------------------------------------------------------------
%
%   Subroutine to correct a specific field in DICOM metadata structure
%
%   INPUT:  
%           metadata - DICOM metadata structure                   'struct'
%           fieldName - name of the field to correct               'string'
%
%   OUTPUT: 
%           metadata - corrected metadata structure with:
%               - If the specified field contains ASCII-encoded float data:
%                   The field is converted to 'double'.
%               - If the field already contains 'single' or 'double':
%                   The field remains unchanged.
%
%   DESCRIPTION:
%       This function checks if a specific field in the provided DICOM 
%       metadata contains ASCII-encoded numeric data (as a numeric array).
%       If the data is detected as numeric but not already a float type,
%       the function converts it to a float using 'str2double'. If the field
%       is already 'single' or 'double', no conversion is performed.
%
%   EXAMPLE USAGE:
%       filename = 'your_dicom_file.dcm';
%       metadata = dicominfo(filename);  % Load the DICOM metadata
%       metadata = correctPrivateField(metadata, 'Private_700d_1000');
%       disp(metadata.Private_700d_1000);  % Verify the correct float value
%
%   AUTHOR:
%       vmalis@ucsd.edu
%--------------------------------------------------------------------------
fieldName = 'Private_700d_1000';
    
    % Check if the field exists in the metadata
    if isfield(metadata, fieldName)
        fieldValue = metadata.(fieldName);
        
        % Only attempt correction if the field is not already single or double
        if ~isa(fieldValue, 'single') && ~isa(fieldValue, 'double')
            if isnumeric(fieldValue)
                % Convert the numeric array to a character array
                charArray = char(fieldValue)';  % Interpret as ASCII characters
                
                % Convert the character array to a float value
                floatValue = str2double(charArray);  % Convert string to float
                
                % If conversion was successful, update the field
                if ~isnan(floatValue)
                    metadata.(fieldName) = floatValue;
                else
                    warning('Failed to convert %s to a float. Setting to 1.0.', fieldName);
                    metadata.(fieldName) = 1.0;  % Set default value
                end
            else
                warning('Field %s is not numeric. Setting to 1.0.', fieldName);
                metadata.(fieldName) = 1.0;  % Set default value
            end
        end
    else
        warning('Field %s does not exist in the metadata. Creating and setting to 1.0.', fieldName);
        metadata.(fieldName) = 1.0;
    end
end