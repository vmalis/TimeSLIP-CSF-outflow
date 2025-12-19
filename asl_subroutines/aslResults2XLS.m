function aslResults2XLS(S, filename, param)
%--------------------------------------------------------------------------
%
%   Function to remove the last six fields from a structure, convert the 
%   modified structure to a table, and save it as an Excel file.
%
%   This function processes a given structure by removing its last six fields 
%   (if available), then converts the modified structure to a table format. 
%   The resulting table is saved to an Excel file with a specified filename.
%
%   INPUT:  
%           S           - Input structure with multiple fields       (struct)
%           filename    - Filename for the output Excel file         (string)
%
%   OUTPUT: 
%           None (The function saves the modified structure as a table in an 
%           Excel file).
%
%   Example Usage:
%           S = struct('Field1', {1, 2, 3}, 'Field2', {4, 5, 6}, 'Field3', {7, 8, 9});
%           filename = 'output.xlsx';
%           aslResults2XLS(S, filename);
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

    % Get all field names from the structure
    fieldNames = fieldnames(S);  
    numFields = length(fieldNames);

    % Remove the last 6 fields if there are more than 6 fields
    if numFields > 6 & ~param.fixedRO
        fieldsToRemove = fieldNames(end-9:end);
        S = rmfield(S, fieldsToRemove);
    end

    % Convert the modified structure to a table
    if size(S,1)==1
         T = struct2table(S,'AsArray', true);
    else
        T = struct2table(S);
    end

    % Save the table to an Excel file with the specified filename
    writetable(T, filename);

end