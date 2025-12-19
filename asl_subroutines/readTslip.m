function [tslipStruct,V,param] = readTslip(path, ROI, param)
%--------------------------------------------------------------------------
%
%   Subroutine to read directory with tSLIP data
%
%   INPUT:  path to directory with DICOM folders                 'string'
%           slice to select                                      'int'
%           flag to average slices                               'boolean'
%   OUTPUT: structure with:
%               ON   images                                      (double)
%               CTRL images (if alternate was acquired)          (double)
%               S    header information                          struct
%               TI   in [ms]                                     double
%           stucture with all images:
%               ON                                               (double)
%               CTRL                                             (double)
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

% Ensure the path is valid
if ~isfolder(path)
    error('Invalid directory: %s', path);
end

% Store the current directory to return to later
originalDir = pwd;

% Change to the specified directory
cd(path);
list3 = dir3();

% Initialize the output structure
tslipStruct = struct('On', {}, 'Ctrl', {}, 'S', {}, 'dicomScaling', {}, ...
                     'TI', {},'On_roi',{},'Ctrl_roi',{});
V = [];
tmpOn=[];
tmpCtrl=[];


for k = 1:length(list3)
    fileName = list3(k).name;

    if strcmp(param.acq_dim,'3d')

        % Define pattern to match the TI value
        pattern = '(\d+)[_ ]?[mM][sS]';
    
        % Use regular expression to find the TI value
        tokens = regexp(fileName, pattern, 'tokens');
    
        if ~isempty(tokens)
            TI = str2double(tokens{1}{1});
        else
            TI = NaN; % Default to NaN if TI value is not found
        end
    
        % Read the DICOM file
        [I, S] = dicom2struct_canon(fileName, 'data');
        I = double(I);
        
    
       
    
        % Split if series was acquired with alternate % Determine the number of slices
        num_slices = numel(unique([S.location]));
        
        if num_slices < size(S, 2)

            On = I(:, :, 1:num_slices);
            Ctrl = I(:, :, num_slices + 1:end);

            tmpOn = cat(4,tmpOn,On);
            tmpCtrl = cat(4,tmpCtrl,Ctrl);
    
            if param.filter_median3D
                On  = medfilt3(On);
                Ctrl = medfilt3(Ctrl);
            end
    
            if param.averageSlices
                On = mean(On, 3);
                Ctrl = mean(Ctrl, 3);
            else
                On = On(:, :, ROI.slice_number);
                Ctrl = Ctrl(:, :, ROI.slice_number);
            end
    
            if param.filter_median2D
                On  = medfilt2(On);
                Ctrl = medfilt2(Ctrl);
            end
    
             tslipStruct(k).On = On;
             tslipStruct(k).Ctrl = Ctrl;
    
        else
    
            On = I;
            tmpOn = cat(4,tmpOn,On);
    
            if param.filter_median3D
                On  = medfilt3(On);
            end
    
            if param.averageSlices
                On = mean(On, 3);
            else
                On = On(:, :, ROI.slice_number);
            end
    
            if param.filter_median2D
                On  = medfilt2(On);
            end
            
            tslipStruct(k).On = On;
    
        end
    
        tslipStruct(k).S = S;
    
        % Extract the DICOM scaling factor
        if isfield(S(1).header, 'Private_700d_1000')
            tslipStruct(k).dicomScaling = S(1).header.Private_700d_1000;
        else
            % Default to NaN if scaling factor is not found
            tslipStruct(k).dicomScaling = NaN; 
        end
    
        tslipStruct(k).TI = TI;
        

    elseif strcmp(param.acq_dim,'2d')

        % Read the DICOM file
        [I, S] = dicom2struct_canon(fileName, 'data');
        I = double(I);

        if strcmp(path,'tslip')
            % Use regular expressions to extract the values
            tokens = regexp(fileName, 'dT=(\d+).*?_S(\d+)ms', 'tokens');

            if ~isempty(tokens)
                numbers = str2double(tokens{1});
                param.TI(2)= numbers(1);
                param.TI(1) = numbers(2);
            else
                param.TI=param.TI0;
            end

            param.TI0
            param.TI
        end

        for m=1:size(I,3)

            On=I(:,:,m);
            tmpOn = cat(4,tmpOn,On);
            
            if param.filter_median2D
                On = medfilt2(On);
            end

            tslipStruct(m).On = On;
            tslipStruct(m).S  = S(m);
           
            % Extract the DICOM scaling factor
            if isfield(S(m).header, 'Private_700d_1000')
                tslipStruct(m).dicomScaling = S(m).header.Private_700d_1000;
            else
                 % Default to NaN if scaling factor is not found
                tslipStruct(m).dicomScaling = NaN;
            end
            
            if param.fixedRO 
                tslipStruct(m).TI = ((size(I,3)-1)*param.TI(2)) - ...
                    (param.TI(1)+param.TI(2)*(m-1));
            else
                tslipStruct(m).TI = param.TI(1)+param.TI(2)*(m-1);
            end
            
        end

    else
        error('Acquisition dimension not specified');
    end

    cd ..
end



% Sort the structure array by TI
[tslipStruct,idx] = sortStruct(tslipStruct, 'TI');

V.On=tmpOn(:,:,:,idx);
if ~isempty(tmpCtrl)
    V.Ctrl=tmpCtrl(:,:,:,idx);
end


if param.dicomScale
    for t=1:size(tmpOn,4)
        V.On(:,:,:,t)=V.On(:,:,:,t)/tslipStruct(t).dicomScaling;
        if ~isempty(tmpCtrl)
            V.Ctrl(:,:,:,t)=V.Ctrl(:,:,:,t)/tslipStruct(t).dicomScaling;
        end
    end
end


% apply ROI and calculate average values per ROI

for i=1:numel(tslipStruct)

    On = repmat(tslipStruct(i).On,[1,1,size(ROI.maskIm,3)]);
   
    if param.filter_sd
        On = arraySDfilter(On.*ROI(1).maskIm,2);
    end

    tslipStruct(i).On_roi = ...
        squeeze(mean(On.*ROI(1).maskIm,[1,2],'omitnan'));

    if isfield(tslipStruct,'Ctrl') && ~isempty(tslipStruct(i).Ctrl)
        Ctrl = repmat(tslipStruct(i).Ctrl,[1,1,size(ROI.maskIm,3)]);
        
        if param.filter_sd
            Ctrl = arraySDfilter(Ctrl.*ROI(1).maskIm,2);
        end

        tslipStruct(i).Ctrl_roi = ...
            squeeze(mean(Ctrl.*ROI(1).maskIm,[1,2],'omitnan'));
    end

end

% Change back to the original directory
cd(originalDir);

end