function X=TslipStuct2Array(T,field,scaling)
%--------------------------------------------------------------------------
%
%   Subroutine to convert a field of a structure array to a 3D array and optionally scale the values
%
%   INPUT:  Structure array with various fields                  'struct'
%           Field name as a string                               'string'
%           Flag to apply scaling                                'boolean'
%   OUTPUT: 3D array of the specified field values               'double'
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------
    
    n=numel(T);
    X=zeros([size(T(1).(field)),n]);

    for i=1:n
        X(:,:,i)=T(i).(field);
        if scaling
            X(:,:,i) = X(:,:,i)/T(i).dicomScaling;
        end
    end

end