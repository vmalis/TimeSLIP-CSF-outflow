%% ========================================================================
% Top level code for perfusion anlysis
%  
%==========================================================================
%
%   05/2025 - VM  for Canon (vmalis@ucsd.edu)       v2.0
%
%==========================================================================
%
%  expected folder strucutre:
%
%    ---|Subject
%         |---pre/post/postN folders
%            required|--- roi     [dicom folder] + (Horos roi *.csv file)
%            required|--- tslip   [dicom folderS]
%            optional|--- Control [dicom folderS] if opp. tag (MT cancel)
%            optional|--- M0      [dicom folder] optional
%                 
%     use ExportRoi Plugin from:
%     https://horosproject.org/horos-content/plugins/horos/ExportROIs/                    
%
%% ========================================================================
clear all
clc

tic

start_directory=pwd;

% SUBJECTS list
list = dir3();

if exist('progress.mat','file')
    load('progress.mat')
    if ii>numel(list)
       disp('All data processed!');
    return
    end

else

%% Define parameters
%--------------------------------------------------------------------------
param.acq_type      = 'bright';     % bright (NonSelect+Select); dark (Select)
%param.acq_dim      = '2d';         % 2d (single T-SLIP series with multiple TI)
param.acq_dim       = '3d';         % 3d (every TI separate series)
%param.TI            = [1000, 500];    
param.dicomScale    = true;         % MUST BE ON
param.perVoxelCalc  = true;         % more prone to error
param.T1            = 2000;         % fixed CSF T1 for bi-component fit:
                                    % Blood use 1500, CSF use 3000 
param.registration = 'all';        % rigid body registration to M0:
                                        % extreme (myocardium)
                                        % all     all images to M0
                                        % auto    automatic
                                        % off
%------acquisition with fixed RO-------------------------------------------
param.fixedRO       = false;      % BBIRprep WIP1
%param.bpm           = 70;         % heart rate
%--------------------------------------------------------------------------
param.filter_median3D = true;       % for 3D volume median filter across slices
param.filter_median2D = true;       % additional median for slice of interest
param.filter_sd       = true;       % doubleSD window for per ROI average
param.noiseM0         = false;      % dont normalize to noise
param.absValue        = false;       % subtraction as absolute value
param.nullPointFilter = [false,750];% remove short TI [bool, filter]
param.filter_lowpass  = 'lowPass';      % 'lowPass' 'movinAver' 'off' 
param.fitplot  = {'bi','gm','biC'}; % 'gkm','gm', 'gmf', 'gs', 'bi' list them in cell
%--------------------------------------------------------------------------
% average across slices in 3D volume ONLY FOR LARGE PERFUSION VOLUME
param.averageSlices = false;
%--------------------------------------------------------------------------
% to print TEs on colormap, use "other" for SeriesID
param.dataType = 'regular'; 
% use false for NIH (Arial font) and true (latex) for publications
param.latex = true;
%--------------------------------------------------------------------------
% colormaps
% varibles with two values 1 - plot/dont plot; 2 - zoom to roi / full FOV
param.visible        = false;           % dont show figures, plot silently
param.allslices      = false;           % if 3D colormaps for all slices
param.paleteview     = [false,false];    % does all TIs
param.subtractionMap = [false,false];    % abs(ON-OFF)
param.montageALL     = [false,false];    % grayscale: ON, OFF, Control
param.sir2Dcmap      = [false,true];   % 2D SIR
param.sir2DcmapFull  = false;           % to overlay entire FOV
param.contrastNorm   = true;            % enchance contrast, not true Signal
%--------------------------------------------------------------------------
param.SIRLim         = 100;               %default limit for plots
param.SIRClim        = 100;               %default limit for colormaps
param.TILim          = 5000;              %custom x axis limits
%--------------------------------------------------------------------------
% only for advanced user, keep these as is
param.nonAlternate = false;                    % if control is acquired sep.
param.advancedGKM=false;                       % advanced GKM Fit


% Check parameters
param = tslipParamCheck(param);

if isfield(param,'TI')
    param.TI0=param.TI;    % for 2D 
end

%% PROCESSING
%==========================================================================

% structure with results: appended at each subject
ResultM = [];   
ResultP = [];

ii=1;
jj=1;


end

%************************************************************ LOOP: SUBJECT
for i=ii:numel(list)
    
    cd(list(i,1).name)

    % GROUP/CATEGORY list (such as activity/pre-post/etc)
	list2 = dir3();

    
  
    %**************************************************LOOP: GROUP/CATEGORY
    for j=jj:numel(list2)

        if jj==1
            RM=[];
            RP=[];
        end
    
        % Subject/Study info    
        Subject.ID      = list(i,1).name;
        Subject.study   = list2(j,1).name;
        q=1; % study results index

        cd(list2(j,1).name)

    %----------------------------------------------------- START: read data
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: read ROI
        cd roi
        csv=dir('*.csv');
        fname = [list(i,1).name,'  ',list2(j,1).name,];
        ROI = horos2matlabTSLIP(csv,fname,true,param);

        % Find the indices of all 1s in the array create combined Full
        ROI.maskF = sum(~isnan(ROI.maskIm), 3) > 0;
        [row, col] = find(ROI.maskF == 1);
        % Find the minimum and maximum row and column indices
        n=10; %padding
        ROI.mtop    = max(min(row) - n, 1);
        ROI.mbottom = min(max(row) + n, size(ROI.maskF, 1));
        ROI.mleft   = max(min(col) - n, 1);
        ROI.mright  = min(max(col) + n, size(ROI.maskF, 2));

        cd ..  % out of 'roi'

        % do MICO for MPRAGE or T1w (MPRAGE is priority)
        if isfolder('MPRAGE')
            cd 'MPRAGE'
            MPRAGE2FASE('MPRAGE','FASE2D')
            seriesPath=fullfile(pwd,'Image4Segmentation');
            evalc('R = aslMICO(seriesPath, param, true);');
            if ~isempty(R)
                ROI = mergeASLroi(ROI, R);
            end
    
            cd ..

        elseif ~isfolder('MPRAGE') && isfolder('T1w')
            seriesPath=fullfile(pwd,'T1w');
            evalc('R = aslMICO(seriesPath, param, true);');
            if ~isempty(R)
                ROI = mergeASLroi(ROI, R);
            end



        end


        clearvars csv fname n  row col
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: read ROI

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: read M0
        if exist('M0','dir')
            [M0,VM0]=readTslip('M0',ROI,param);
            if ~isempty(VM0) && strcmp(param.acq_dim,'3d')  %only created for 3d
                 V.M0 = VM0.On;
            end
            if isempty(VM0) && strcmp(param.acq_dim,'2d')  %only created for 3d
                M0=M0(:,:,1);
            end
            % M0 should only have single series
            if isfield(M0,'Ctrl') && isfield(M0,'Ctrl_roi') ...
                && ~isempty([M0.Ctrl]) && ~isempty(Ctrl_roi)

                M0.On     = M0.Ctrl;
                M0.On_roi = M0.Ctrl_roi;
                M0=rmfield(M0,'Ctrl');
                M0=rmfield(M0,'Ctrl_roi');
            end
            clearvars VM0
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: read M0

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: read CONTROL
        if exist('Control','dir')
            [Control,VC]=readTslip('control',ROI,param);
            if ~isempty(VC)&& strcmp(param.acq_dim,'3d')  %only created for 3d
                V.Ctrl=VC.On;
            end
            % control should only have ON series
            if isfield(Control,'Ctrl') && isfield(Control,'Ctrl_roi')
                Control=rmfield(Control,'Ctrl');
                Control=rmfield(Control,'Ctrl_roi');
            end
            clearvars VC
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: read Control

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: read TSLIP
        if exist('tslip','dir')
            [TSlip,VT,param]=readTslip('tslip',ROI,param);
            if ~isempty(VT)&& strcmp(param.acq_dim,'3d')  %only created for 3d
                V.On=VT.On;
                if isfield(VT,'Ctrl')
                    V.Ctrl=VT.Ctrl;
                end
                if ~isfield(V,'M0') && strcmp(param.acq_type,'bright') && ...
                        strcmp(param.acq_dim,'3d') && ~isfield(VT,'Ctrl')
                    V.M0 = V.Ctrl(:,:,:,end);
                    param.nonAlternate=true;
                elseif ~isfield(V,'M0') && strcmp(param.acq_type,'dark') && ...
                        strcmp(param.acq_dim,'3d') && ~isfield(VT,'Ctrl')
                    V.M0 = V.Ctrl;
                    param.nonAlternate=true;
                elseif ~isfield(V,'M0') && strcmp(param.acq_type,'bright') && ...
                        strcmp(param.acq_dim,'3d')
                    V.M0 = VT.Ctrl(:,:,:,end);
                elseif ~isfield(V,'M0') && strcmp(param.acq_type,'dark') && ...
                        strcmp(param.acq_dim,'3d')
                    V.M0 = VT.Ctrl;
                end
            end
            clearvars VT
        else
            error('Data should have tslip folder');
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: read TSLIP

	    
    %------------------------------------------------------- END: read data
    nTIs=numel(TSlip);    


    %------------------------------------------------ START: ararys for SIR
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: denominator
        if exist('M0','var')
            A.M0=TslipStuct2Array(M0,'On',param.dicomScale);
            m0=squeeze(TslipStuct2Array(M0,'On_roi',...
                                            param.dicomScale));
        elseif ~exist('M0','var') && strcmp(param.acq_dim,'3d') ...
            && ~param.nonAlternate
            A.M0=repmat(TslipStuct2Array(TSlip(end),'Ctrl',param.dicomScale),[1,1,nTIs]);
            m0=squeeze(TslipStuct2Array(TSlip(end),'Ctrl_roi',...
                                            param.dicomScale));
        elseif ~exist('M0','var') && strcmp(param.acq_dim,'3d') ...
                && param.nonAlternate
            A.M0=repmat(TslipStuct2Array(Control(end),'On',param.dicomScale),[1,1,nTIs]);
            m0=squeeze(TslipStuct2Array(TSlip(end),'On_roi',...
                                            param.dicomScale));
        elseif ~exist('M0','var') && strcmp(param.acq_dim,'2d')
            A.M0=repmat(TslipStuct2Array(Control(end),'On',param.dicomScale),[1,1,nTIs]);
            m0=squeeze(TslipStuct2Array(TSlip(end),'On_roi',...
                                            param.dicomScale));
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: denominator
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: On
        A.On=TslipStuct2Array(TSlip,'On',param.dicomScale);
        on=squeeze(TslipStuct2Array(TSlip,'On_roi',param.dicomScale));
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: On
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: Control
        if exist('Control','var')
            A.Ctrl=TslipStuct2Array(Control,'On',param.dicomScale);
            ctrl=squeeze(TslipStuct2Array(Control,'On_roi',param.dicomScale));
        else
            A.Ctrl=TslipStuct2Array(TSlip,'Ctrl',param.dicomScale);
            ctrl=squeeze(TslipStuct2Array(TSlip,'Ctrl_roi',...
                                            param.dicomScale));
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: Control

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: check m0
        if size(m0,2)~=size(on,2)
            m0=repmat(m0,[1,size(on,2)]);
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: check m0


    %-------------------------------------------------- END: ararys for SIR

    %%-------------------------------------check if pre-null point to be removed

    % check dims of M0 and m0 cause sometimes they are acquired seperat;y
    % with less dT

    if size(A.M0,3)>1 && size(A.M0,3) ~= size(A.On,3)
             A.M0=A.M0(:,:,1);
    end


    if min(size(m0))>1
        if size(m0, 2) > 1 && size(m0,2) ~= size(on,2)
                m0=m0(:,1);
        end
    else
        if size(m0, 1) > 1 && size(m0,1) ~= size(on,1)
            m0=m0(1);
        end
    end


    if param.nullPointFilter(1)
    
        % Get TI values from the structure and find indices below the threshold
        TI_values = [TSlip.TI];  
        idx = TI_values < param.nullPointFilter(2);  % Logical indexing is more efficient than `find`
    
        
        % Conditional checks before modifying arrays to avoid errors
        if exist('Control','var')
            Control(idx) = [];
        end

        if exist('M0', 'var') && size(M0, 2) > 1 && size(M0, 2)== size(TSlip, 2)
            M0(idx) = [];
        end

        % Remove entries based on identified indices
        TSlip(idx) = [];
    
        if min(size(m0))>1
            if size(m0, 2) > 1 && size(m0,2) == size(on,2)
                m0(:, idx) = [];
            else
                m0=m0(:,1);
            end
        elseif size(ROI.name,1)>1 && size(m0,2) == 1

        else
            if size(m0, 1) > 1 && size(m0,1) == size(on,1)
                m0(idx) = [];
            else
                m0=m0(1);
            end
        end


        if size(size(on))<2
            on(idx) = [];
            ctrl(idx) = [];
        else
            on(:, idx) = [];
            ctrl(:, idx) = [];
        end
        
        if isfield(A, 'M0') && size(A.M0,3)==size(A.On,3)
            A.M0(:, :, idx) = [];
        else
            
        end
    
        % Update fields within structure A
        A.On(:, :, idx) = [];
        A.Ctrl(:, :, idx) = [];

    
        % Check acquisition dimension and update variables accordingly
        if strcmp(param.acq_dim, '3d')
            V.On(:, :, :, idx) = [];
            V.Ctrl(:, :, :, idx) = [];
            
            if isfield(V, 'M0') && size(V.M0, 4) > 1
                V.M0(:, :, :, idx) = [];
            end
        elseif strcmp(param.acq_dim, '2d')
            param.TI(1) = TI_values(1);
        end
    end

    nTIs=max(size([TSlip.TI]));  

    clearvars idx TI_values
    %%-------------------------------------------------------------------------
    if ~strcmpi(param.registration,'off')
        A=ASL_imReg(A, param);
    end






    %-------------------------------------------------- START: get noise
    %masks for colormaps
        if strcmp(param.acq_dim,'3d') && ~isempty(V)  
            ROI.mask=UTEsMask(mat2gray(max(V.Ctrl,[],[3,4])),0.02);
            if param.noiseM0
                ROI.noiseLevel = std(V.On(V.On.*abs(ROI.mask-1) ~= 0), 0,"all", "omitnan");
                temp = ones(size(V.On));
                temp(abs(V.On-V.Ctrl) < ROI.noiseLevel/sqrt(2)) = NaN;
                ROI.Noise=temp.*ROI.mask;
            end
        else
            ROI.mask=UTEsMask(mat2gray(max(A.Ctrl,[],[3,4])),0.02);
            if param.noiseM0
                ROI.noiseLevel = std(A.On(A.On.*abs(ROI.mask-1) ~= 0),0, "all", "omitnan");
                temp = ones(size(A.On)); 
                temp(abs(A.On-A.Ctrl) < ROI.noiseLevel/sqrt(2)) = NaN;
                ROI.Noise=temp.*ROI.mask;
            end
        end

        % mean recalcs with Noise mask
        if param.noiseM0
            on=squeeze(mean(permute(repmat(A.On.*ROI.Noise,...
                [1,1,1,size(ROI.maskIm,3)]),[1,2,4,3]).*repmat...
                (ROI.maskIm,[1,1,1,size(A.On,3)]),[1,2],'omitnan'));
            ctrl=squeeze(mean(permute(repmat(A.Ctrl.*ROI.Noise,...
                [1,1,1,size(ROI.maskIm,3)]),[1,2,4,3]).*repmat...
                (ROI.maskIm,[1,1,1,size(A.Ctrl,3)]),[1,2],'omitnan'));
            m0=squeeze(mean(permute(repmat(A.M0.*ROI.Noise,...
                [1,1,1,size(ROI.maskIm,3)]),[1,2,4,3]).*repmat...
                (ROI.maskIm,[1,1,1,size(A.M0,3)]),[1,2],'omitnan'));

            on(isnan(on))=0;
            ctrl(isnan(ctrl))=0;
            m0(isnan(m0))=0;

            on(isinf(on))=0;
            ctrl(isinf(ctrl))=0;
            m0(isinf(m0))=0;
        end


        clearvars temp
     %-------------------------------------------------- END: get noise


    %------------------------------------------------- START: calculate SIR
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: full array
        if strcmp(param.acq_dim,'3d') && ~isempty(V)
            
            tmpOn   = V.On;
            tmpCtrl = V.Ctrl;
            tmpM0   = V.M0;

            if param.noiseM0
                tmpM0=tmpM0.*ROI.Noise;
            end

            if param.filter_median3D
                tmpOn = arrayfun(@(i) medfilt3(tmpOn(:,:,:,i)), ...
                              1:size(tmpOn,4), 'UniformOutput', false);
                tmpOn = cat(4, tmpOn{:});
                tmpCtrl = arrayfun(@(i) medfilt3(tmpCtrl(:,:,:,i)), ...
                              1:size(tmpCtrl,4), 'UniformOutput', false);
                tmpCtrl = cat(4, tmpCtrl{:});
                tmpM0 = arrayfun(@(i) medfilt3(tmpM0(:,:,:,i)), ...
                              1:size(tmpM0,4), 'UniformOutput', false);
                tmpM0 = cat(4, tmpM0{:});
            end
            
            if param.absValue
                V.SUBT=abs(tmpOn-tmpCtrl);
            elseif ~param.absValue && strcmp(param.acq_type,'bright')
                V.SUBT=tmpOn-tmpCtrl;
                V.SUBT(V.SUBT<0)=0;
            else
                V.SUBT=tmpCtrl-tmpOn;
                V.SUBT(V.SUBT<0)=0;
            end
 
            if param.noiseM0
                V.SUBT=V.SUBT.*ROI.Noise;
            end

            V.SIR=V.SUBT./tmpM0;

            V.SUBT(isinf(V.SUBT))=NaN;
            V.SIR(isinf(V.SIR))=NaN;
            


            clearvars tmp tmpOn tmpCtrl tmpM0
        end


        if param.filter_median2D
            % Preallocate arrays with the same size as A.On
            tmpOn = zeros(size(A.On));
            tmpCtrl = tmpOn;
            
            % Check if A.M0 has multiple slices and preallocate accordingly
            if size(A.M0, 3) > 1
                tmpM0 = zeros(size(A.M0));
            end
            
            % Apply median filtering across all time points (nTIs)
            for t = 1:nTIs
                tmpOn(:, :, t) = medfilt2(A.On(:, :, t));
                tmpCtrl(:, :, t) = medfilt2(A.Ctrl(:, :, t));
                
                if size(A.M0, 3) > 1
                    tmpM0(:, :, t) = medfilt2(A.M0(:, :, t));
                end
            end
            
            % Handle the case where A.M0 has only one slice
            if size(A.M0, 3) == 1
                tmpM0 = medfilt2(A.M0);
            end

        else
                tmpOn=A.On;
                tmpCtrl=A.Ctrl;
                tmpM0=A.M0;
        end

        if param.absValue
            A.SUBT=abs(tmpOn-tmpCtrl);
        elseif ~param.absValue && strcmp(param.acq_type,'bright')
            A.SUBT=tmpOn-tmpCtrl;
            A.SUBT(A.SUBT<0)=0;
        else
            A.SUBT=tmpCtrl-tmpOn;
            A.SUBT(A.SUBT<0)=0;
        end

        if param.noiseM0
            A.SUBT = A.SUBT.*ROI.Noise;
        end
        
        A.SIR=A.SUBT./tmpM0;

        A.SUBT(isinf(A.SUBT))=NaN;
        A.SIR(isinf(A.SIR))=NaN;


        clearvars tmpOn tmpCtrl tmpM0 t
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: full array
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start: roi quantification
        
        if param.absValue
            SIR.mean=abs(on-ctrl)./m0;
        elseif ~param.absValue && strcmp(param.acq_type,'bright')
            temp=(on-ctrl)./m0;
            temp(temp<0)=0;
            SIR.mean=temp;
        else 
            temp=(ctrl-on)./m0;
            temp(temp<0)=0;
            SIR.mean=temp;
        end
        clearvars on ctrl m0 tmp

        if param.perVoxelCalc
            pp=[];
            for m=1:size(ROI.maskIm,3)
                tmp=A.SIR.*repmat(ROI.maskIm(:,:,m),[1,1,nTIs]);
                tmp=(mean(tmp,[1,2],'omitnan'));
                pp=cat(1,pp,tmp);
            end
            SIR.pp=squeeze(pp);
        end

        clearvars tmp pp m
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End: roi quantification
    %--------------------------------------------------- END: calculate SIR

    %----------------------------------------------------------- START: fit

    if param.fixedRO
        if param.perVoxelCalc
            [rM,rP]=aslFixedRO(SIR,[TSlip.TI],Subject,ROI,param,A.M0(:,:,1));
        else
            rM=aslFixedRO(SIR,[TSlip.TI],Subject,ROI,param,A.M0(:,:,1));
        end
    else
        if param.perVoxelCalc
            [rM,rP]=aslFit(SIR,[TSlip.TI],Subject,ROI,param,A.M0(:,:,1));
        else
            rM=aslFit(SIR,[TSlip.TI],Subject,ROI,param,A.M0(:,:,1));
        end

    end
    %------------------------------------------------------------- END: fit
    cd ..
    %----------------------------------------------------- START: colormaps
    mkdir('colormaps')
    cd('colormaps')

        if strcmp(param.acq_dim, '3d') && ~isempty(V) && param.allslices
            tempSUBT = mat2gray(V.SUBT .* ROI.mask);
            slice_indices = 1:size(V.On, 3);
            data_source = V;
        else
            slice_indices = ROI.slice_number;
            data_source = A;
        end

        for s = slice_indices
            if strcmp(param.acq_dim, '3d') && ~isempty(V) && param.allslices
    
                

                B = structfun(@(field) squeeze(field(:,:,s,:)), V,...
                    'UniformOutput', false);
                B.SUBT = squeeze(tempSUBT(:,:,s,:));
            else
                B = A;
            end

            operations = {
            'palette', @paleteCMap, param.paleteview(1);
            'montage', @tslipmontage, param.montageALL(1);
            'subtraction', @(B, ROI, TI, s, param) perfusionCMap(B, ...
                ROI, TI, s, param, 'SUBT'), param.subtractionMap(1);
            'SIR', @(B, ROI, TI, s, param) perfusionCMap(B, ROI, TI, ...
                s, param, 'SIR'), param.sir2Dcmap(1);
            };

            for x = 1:size(operations, 1)
                if operations{x, 3}
                    mkdir(operations{x, 1});
                    cd(operations{x, 1});
                    operations{x, 2}(B, ROI, [TSlip.TI], s, param);
                    cd ..;
                end
            end
        clearvars B x operations slice_indices s data_source
        end
        cd ..  
    %------------------------------------------------------- END: colormaps
       
        if param.perVoxelCalc
            RM=[RM,rM];
            RP=[RP,rP];
            clearvars rM rP
            ii=i;
            jj=j+1;
            save([start_directory,'/progress.mat'],'ResultM','ResultP',...
                'RM','RP','ii','jj','list','list2','param')
        else
            RM=[RM,rM];
            clearvars rM
            ii=i;
            jj=j+1;
            save([start_directory,'/progress.mat'],'ResultM','RM','ii','jj',...
            'list','list2','param')
        end


        cd ..

    end

    %**************************************************LOOP: GROUP/CATEGORY

    if param.perVoxelCalc
        ResultM=[ResultM,RM];
        ResultP=[ResultP,RP];
        clearvars RM RP
        ii=i+1;
        jj=1;
        save([start_directory,'/progress.mat'],'ResultM','ResultP','ii','jj',...
            'list','list2','param')
    else
        ResultM=[ResultM,RM];
        clearvars RM
        ii=i+1;
        jj=1;
        save([start_directory,'/progress.mat'],'ResultM','ii','jj',...
            'list','list2','param')
    end
cd ..

end
%************************************************************ LOOP: SUBJECT  



data_directory=pwd;

folder_results = ['results_', char(datetime('today'), 'yyMMdd')];
mkdir(folder_results);

movePDFs(data_directory, folder_results)

cd(folder_results)

if param.perVoxelCalc
    aslResults2XLS(ResultM, 'ResultMean.xlsx',param)
    aslResults2XLS(ResultP, 'ResultPP.xlsx',param)
    save('results.mat', 'ResultM','ResultP')
else
    aslResults2XLS(ResultM, 'ResultMean.xlsx',param)
    save('results.mat', 'ResultM')
end

toc