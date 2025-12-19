function varargout = aslFit(SIR,ti,Subject,ROI,param,Image)

%--------------------------------------------------------------------------
%
%   Function to perform ASL data fitting and analysis
%
%   This function processes ASL (Arterial Spin Labeling) data to fit 
%   various models (Gamma variate, Gaussian, and combined fits) to the 
%   provided SIR (Signal Intensity Ratio) data. The function outputs the 
%   results for both mean and per-pixel data within specified ROIs.
%
%   INPUT:
%           SIR     - Structure containing ASL signal intensity ratio data.
%                     Fields:
%                       SIR.mean - Mean SIR values (1D array or 2D matrix)
%                       SIR.pp   - Per pixel SIR values (optional)
%           ti      - Inversion times (1D array)                       (double)
%           Subject - Structure containing subject information:
%                       Subject.ID    - Subject ID (string)
%                       Subject.study - Study group (string)
%           ROI     - Structure defining the Region of Interest (ROI):
%                       ROI.maskIm     - ROI mask image (3D matrix)
%                       ROI.name       - ROI names (cell array)
%                       ROI.slice_number - Slice number (integer)
%           param   - Structure of parameters:
%                       param.acq_type      - Acquisition type (string)
%                       param.perVoxelCalc  - Per voxel calculation flag (logical)
%                       param.filter_lowpass- Apply low pass filter (logical)
%                       param.visible       - Plot visibility ('on'/'off')
%                       param.T1            - T1 value (double)
%                       param.latex         - LaTeX formatting flag (logical)
%           Image   - Image data for displaying ROIs (2D matrix)
%
%   OUTPUT:
%           varargout{1} - ResultsM: Structure array containing results for mean data.
%                           Fields:
%                               ID, Group, ROI, TTP, FWHM, PH, AUC, Max, T1, f,
%                               fGauss, fGamma, fCSF, bTTP, bFWHM, bPH, bAUC,
%                               R2, RMSE, bR2, bRMSE, ti, y, x_fit, y_fit,
%                               bx_fit, by_fit.
%           varargout{2} - ResultsP (optional): Structure array containing results for per-pixel data.
%                           Fields are the same as in ResultsM.
%
%   Additional Notes:
%       - The function creates plots for visualizing the fits and saves them as PDFs.
%       - The output PDF filenames include the slice number and ROI name.
%       - The function handles both mean and per-pixel analysis depending on the input parameters.
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------



if strcat(param.acq_type,'bright')
    PerfusionK = 6000/0.976;
else
    PerfusionK = 6000;
end

Rm=struct();            % fit results for mean
Rpp=struct();           % fit results for per pixel
ResultsM = struct();    % results for mean
ResultsP = struct();    % results for per pixel

% deafult max TI [ms] for plots
if max(ti)<4000
    TImax=4000;
else 
    TImax=max(ti);
end

I = cropAndResizeImage(mat2gray(Image), ROI);

mkdir('plots');
cd('plots')

for ii=1: size(ROI.maskIm,3)

    fgH = figure('Position', [100, 100, 800, 800],'Visible','off');
    set(fgH, 'Color', 'w');
    set(fgH, 'Visible', param.visible);
    
    % Create a 2x2 tiled layout
    tiledLayoutHandle = tiledlayout(2, 2);

    % Initialize handles array
    ax = gobjects(1, 4);
    
    if min(size(SIR.mean))==1
            sir_mean = SIR.mean';         
    else
            sir_mean = SIR.mean(ii,:);   
    end


    % show roi in HighRes Square CROP 
    mask = cropAndResizeImage(ROI.maskIm(:,:,ii), ROI);
    mask(isnan(mask))=0;
    mask=logical(mask);

    if isfield(SIR,'pp') && param.perVoxelCalc

        % Set all axes to be square
        for jj = 1:4
            ax(jj) = nexttile(tiledLayoutHandle);
            axis(ax(jj), 'square');
        end

        tiledLayoutHandle.TileSpacing = 'tight';
        tiledLayoutHandle.Padding = 'tight';

        displayROIImage(I,mask,ROI.name{ii},0,[0,1,0],'', param, ax(1))

        if min(size(SIR.pp))==1
            sir_pp = SIR.pp';
        else
            sir_pp = SIR.pp(ii,:);   
        end

        % low pass filter
        [ym,ypp] = aslRawPlot(ti, sir_mean, param, ... 
            '', ax(2), sir_pp);


        ym(isnan(ym) | isinf(ym)) = 0;
        ypp(isnan(ypp) | isinf(ypp)) = 0;

        if ~strcmpi(param.filter_lowpass, 'off')

            % GKM numeric
            if param.advancedGKM
            [Rpp(ii).fitGKM,Rpp(ii).gofGKM,Rpp(ii).xGKM,Rpp(ii).yGKM] = ...
                                aslGKMFit3(ti,ypp,param,TImax);
            else
            [Rpp(ii).fitGKM,Rpp(ii).gofGKM,Rpp(ii).xGKM,Rpp(ii).yGKM] = ...
                                aslGKMFit(ti,ypp,param,TImax);
            end    

            % Gamma
            [Rpp(ii).fitGm,Rpp(ii).gofGm,Rpp(ii).xGm,Rpp(ii).yGm] = ...
                                         aslGammaFitCSF2(ti,ypp,TImax);

             % Gamma fixed
            [Rpp(ii).fitGmF,Rpp(ii).gofGmF,Rpp(ii).xGmF,Rpp(ii).yGmF] = ...
                                 aslGammaFitCSF2(ti,ypp,TImax,'-T1',param.T1);

            % Gauss
            [Rpp(ii).fitGs,Rpp(ii).gofGs,Rpp(ii).xGs,Rpp(ii).yGs] = ...
                aslGaussianFit(ti,ypp,TImax);
            
            % Combined: fit of combination rather than combination of fits
            [Rpp(ii).fitC, Rpp(ii).gofC, Rpp(ii).xGc, Rpp(ii).yGc, ...
                Rpp(ii).fGs, Rpp(ii).fGm,Rpp(ii).fit_bc] = ...
                aslCombinedFit3(ti, ypp, param, TImax);

            Rpp(ii).yGc(Rpp(ii).yGc<0)=0;
            measured = [ti;ypp]';
            Rpp(ii).measured = measured;

            maxPP = max(ypp);

        else

            % GKM numeric
            if param.advancedGKM
                [Rpp(ii).fitGKM,Rpp(ii).gofGKM,Rpp(ii).xGKM,Rpp(ii).yGKM] = ...
                            aslGKMFit3(ti,sir_pp,param,TImax);
            else
                [Rpp(ii).fitGKM,Rpp(ii).gofGKM,Rpp(ii).xGKM,Rpp(ii).yGKM] = ...
                            aslGKMFit(ti,sir_pp,param,TImax);
            end

            % Gamma
            [Rpp(ii).fitGm,Rpp(ii).gofGm,Rpp(ii).xGm,Rpp(ii).yGm] = ...
                                         aslGammaFitCSF2(ti,sir_pp,TImax);

            % Gamma fixed
            [Rpp(ii).fitGmF,Rpp(ii).gofGmF,Rpp(ii).xGmF,Rpp(ii).yGmF] = ...
                                 aslGammaFitCSF2(ti,sir_pp,TImax,'-T1',param.T1);

            % Gauss
            [Rpp(ii).fitGs,Rpp(ii).gofGs,Rpp(ii).xGs,Rpp(ii).yGs] = ...
                aslGaussianFit(ti,sir_pp,TImax);
            % Combined
            [Rpp(ii).fitC, Rpp(ii).gofC, Rpp(ii).xGc, Rpp(ii).yGc, ...
                Rpp(ii).fGs, Rpp(ii).fGm,Rpp(ii).fit_bc] = ...
                    aslCombinedFit3(ti, sir_pp,param, TImax);
            Rpp(ii).yGc(Rpp(ii).yGc<0)=0;
            measured = [ti;sir_pp]';
            Rpp(ii).measured = measured;

            maxPP=max(sir_pp);
         
        end

        GKM      = [Rpp(ii).xGKM',Rpp(ii).yGKM'];
        Gamma    = [Rpp(ii).xGm',Rpp(ii).yGm];
        GammaF   = [Rpp(ii).xGmF',Rpp(ii).yGmF];
        Gauss    = [Rpp(ii).xGs',Rpp(ii).yGs];
        Combined = [Rpp(ii).xGc',Rpp(ii).yGc];
        CombinedComponents = [Rpp(ii).fit_bc.gauss, Rpp(ii).fit_bc.gamma];

        CombinedFracions = 100*[Rpp(ii).fGs;Rpp(ii).fGm];
        Ratio = Rpp(ii).fit_bc.PHGauss/Rpp(ii).fit_bc.PHGamma;
        GOF=[Rpp(ii).gofGKM.rsquare,Rpp(ii).gofGm.rsquare,...
            Rpp(ii).gofGmF.rsquare,Rpp(ii).gofGs.rsquare,Rpp(ii).gofC.rsquare];
        
        % plot fits
        aslFitPlot(measured, GKM,Gamma,GammaF,Gauss, Combined,...
            CombinedComponents, GOF, Ratio, param, 'per pixel',ax(3));

        [TTP, FWHM, PH, AUC] = aslFitParams(Rpp(ii).xGm,Rpp(ii).yGm);
        [TTPF, FWHMF, PHF, AUCF] = aslFitParams(Rpp(ii).xGmF,Rpp(ii).yGmF);
        [TTPG, FWHMG, PHG, AUCG] = aslFitParams(Rpp(ii).xGs,Rpp(ii).yGs);
        [bTTP, bFWHM, bPH, bAUC] = aslFitParams(Rpp(ii).xGc, Rpp(ii).yGc);

        ResultsP(ii).ID                    =   Subject.ID;
        ResultsP(ii).Group                 =   Subject.study;
        ResultsP(ii).ROI                   =   ROI.name{ii};
        ResultsP(ii).Max                   =   maxPP*100;
        % GKM
        ResultsP(ii).gkmT1                 =   Rpp(ii).fitGKM(1);
        ResultsP(ii).Delta_t               =   Rpp(ii).fitGKM(2);
        ResultsP(ii).tau                   =   Rpp(ii).fitGKM(3);
        ResultsP(ii).PerfusionGKM          =   Rpp(ii).fitGKM(4)*PerfusionK;
        % Gamma-variate
        ResultsP(ii).deltaTau              =   Rpp(ii).fitGm.dt;
        ResultsP(ii).TTPgm                 =   TTP;
        ResultsP(ii).FWHMgm                =   FWHM;
        ResultsP(ii).PHgm                  =   PH*100;
        ResultsP(ii).AUCgm                 =   AUC*0.1;
        ResultsP(ii).T1gm                  =   Rpp(ii).fitGm.T1;
        ResultsP(ii).PerfusionGamma        =   Rpp(ii).fitGm.perfC*PerfusionK;
        % Gamma-variate fixed T1
        ResultsP(ii).deltaTauF             =   Rpp(ii).fitGmF.dt;
        ResultsP(ii).TTPgmF                =   TTPF;
        ResultsP(ii).FWHMgmF               =   FWHMF;
        ResultsP(ii).PHgmF                 =   PHF*100;
        ResultsP(ii).AUCgmF                =   AUCF*0.1;
        ResultsP(ii).T1gmF                 =   param.T1;
        ResultsP(ii).PerfusionGammaFixedT1 =   Rpp(ii).fitGmF.perfC*PerfusionK;
        % Gaussian
        ResultsP(ii).deltaT                =  Rpp(ii).fitGs.b-2*Rpp(ii).fitGs.c;
        ResultsP(ii).TTPgs                 =  TTPG;
        ResultsP(ii).FWHMgs                =  FWHMG;
        ResultsP(ii).PHgs                  =  PHG*100;
        ResultsP(ii).AUCgs                 =  AUCG*0.1;
        % bi-component
        ResultsP(ii).bi_deltaTGauss        =   Rpp(ii).fit_bc.deltaTGauss;
        ResultsP(ii).bi_deltaTGamma        =   Rpp(ii).fit_bc.deltaTGamma;  
        ResultsP(ii).fGauss                =   CombinedFracions(1);
        ResultsP(ii).fGamma                =   CombinedFracions(2);
        ResultsP(ii).biRatio               =   Rpp(ii).fit_bc.ratio;
        ResultsP(ii).biPerfusion           =   Rpp(ii).fit_bc.perfC*PerfusionK;
        ResultsP(ii).biTTP                 =   bTTP;
        ResultsP(ii).biTTPGauss            =   Rpp(ii).fit_bc.TTPGauss;
        ResultsP(ii).biTTPGamma            =   Rpp(ii).fit_bc.TTPGamma;
        ResultsP(ii).biFWHM                =   bFWHM;
        ResultsP(ii).biPH                  =   bPH*100;
        ResultsP(ii).biPHGauss             =   Rpp(ii).fit_bc.PHGauss*100;
        ResultsP(ii).biPHGamma             =   Rpp(ii).fit_bc.PHGamma*100;
        ResultsP(ii).biAUC                 =   bAUC*0.1;
        % Goodness of fit
        ResultsP(ii).R2_GKM                =   Rpp(ii).gofGKM.rsquare;
        ResultsP(ii).RMSE_GKM              =   Rpp(ii).gofGKM.rmse;
        ResultsP(ii).SSE_GKM               =   Rpp(ii).gofGKM.ssr;
        ResultsP(ii).R2_GV                 =   Rpp(ii).gofGm.rsquare;
        ResultsP(ii).RMSE_GV               =   Rpp(ii).gofGm.rmse;
        ResultsP(ii).SSE_GV                =   Rpp(ii).gofGm.sse;
        ResultsP(ii).R2_GVfT1              =   Rpp(ii).gofGmF.rsquare;
        ResultsP(ii).RMSE_GVfT1            =   Rpp(ii).gofGmF.rmse;
        ResultsP(ii).SSE_GVfT1             =   Rpp(ii).gofGmF.sse;
        ResultsP(ii).R2_Gauss              =   Rpp(ii).gofGs.rsquare;
        ResultsP(ii).RMSE_Gauss            =   Rpp(ii).gofGs.rmse;
        ResultsP(ii).SSE_Gauss             =   Rpp(ii).gofGs.sse;
        ResultsP(ii).biR2                  =   Rpp(ii).gofC.rsquare;
        ResultsP(ii).biRMSE                =   Rpp(ii).gofC.rmse;
        ResultsP(ii).biSSE                 =   Rpp(ii).gofC.sse;
        % actual values and fits x (TI) and y (SIR)
        ResultsP(ii).ti                    =   ti(:); 
        ResultsP(ii).y                     =   sir_pp(:);
        ResultsP(ii).y_filtered            =   ypp(:);
        ResultsP(ii).x_fit_GKM             =   Rpp(ii).xGKM(:);
        ResultsP(ii).y_fit_GKM             =   Rpp(ii).yGKM(:);
        ResultsP(ii).x_fit_GV              =   Rpp(ii).xGm(:);
        ResultsP(ii).y_fit_GV              =   Rpp(ii).yGm(:);
        ResultsP(ii).x_fit_GVfT1           =   Rpp(ii).xGmF(:);
        ResultsP(ii).y_fit_GVfT1           =   Rpp(ii).yGmF(:);
        ResultsP(ii).x_fit_Gauss           =   Rpp(ii).xGs(:);
        ResultsP(ii).y_fit_Gauss           =   Rpp(ii).yGs(:);
        ResultsP(ii).x_fit_bi              =   Rpp(ii).xGc(:);
        ResultsP(ii).y_fit_bi              =   Rpp(ii).yGc(:);



    else
        % axis handle
        ax(1) = nexttile(tiledLayoutHandle);
        ax(2) = nexttile(tiledLayoutHandle);
        ax(3) = nexttile(tiledLayoutHandle,[1,2]);
        
        displayROIImage(I,mask,ROI.name{ii},0,[0,1,0],'', param, ax(1))

        % low pass filter
        ym = aslRawPlot(ti, sir_mean, param,'', ...
            ax(2));

    end

    if ~strcmpi(param.filter_lowpass, 'off')

        % GKM numeric
        if param.advancedGKM
        [Rm(ii).fitGKM,Rm(ii).gofGKM,Rm(ii).xGKM,Rm(ii).yGKM] = ...
                            aslGKMFit3(ti,ym,param,TImax);
        else
        [Rm(ii).fitGKM,Rm(ii).gofGKM,Rm(ii).xGKM,Rm(ii).yGKM] = ...
                            aslGKMFit(ti,ym,param,TImax);
        end

        % Gamma
        [Rm(ii).fitGm,Rm(ii).gofGm,Rm(ii).xGm,Rm(ii).yGm] = ...
                                     aslGammaFitCSF2(ti,ym,TImax);

         % Gamma fixed
        [Rm(ii).fitGmF,Rm(ii).gofGmF,Rm(ii).xGmF,Rm(ii).yGmF] = ...
                             aslGammaFitCSF2(ti,ym,TImax,'-T1',param.T1);

        % Gauss
        [Rm(ii).fitGs,Rm(ii).gofGs,Rm(ii).xGs,Rm(ii).yGs] = ...
            aslGaussianFit(ti,ym,TImax);
        
        % Combined: fit of combination rather than combination of fits
        [Rm(ii).fitC, Rm(ii).gofC, Rm(ii).xGc, Rm(ii).yGc, ...
            Rm(ii).fGs, Rm(ii).fGm,Rm(ii).fit_bc] = ...
            aslCombinedFit3(ti, ym, param, TImax);

        Rm(ii).yGc(Rm(ii).yGc<0)=0;
        measured = [ti;ym]';
        Rm(ii).measured = measured;
        maxM=max(ym);
    else

       % GKM numeric
        if param.advancedGKM
        [Rm(ii).fitGKM,Rm(ii).gofGKM,Rm(ii).xGKM,Rm(ii).yGKM] = ...
                            aslGKMFit3(ti,sir_mean,param,TImax);
        else
        [Rm(ii).fitGKM,Rm(ii).gofGKM,Rm(ii).xGKM,Rm(ii).yGKM] = ...
                            aslGKMFit(ti,sir_mean,param,TImax);
        end

        % Gamma
        [Rm(ii).fitGm,Rm(ii).gofGm,Rm(ii).xGm,Rm(ii).yGm] = ...
                                     aslGammaFitCSF2(ti,sir_mean,TImax);

         % Gamma fixed
        [Rm(ii).fitGmF,Rm(ii).gofGmF,Rm(ii).xGmF,Rm(ii).yGmF] = ...
                             aslGammaFitCSF2(ti,sir_mean,TImax,'-T1',param.T1);

        % Gauss
        [Rm(ii).fitGs,Rm(ii).gofGs,Rm(ii).xGs,Rm(ii).yGs] = ...
            aslGaussianFit(ti,sir_mean,TImax);
        
        % Combined: fit of combination rather than combination of fits
        [Rm(ii).fitC, Rm(ii).gofC, Rm(ii).xGc, Rm(ii).yGc, ...
            Rm(ii).fGs, Rm(ii).fGm,Rm(ii).fit_bc] = ...
            aslCombinedFit3(ti, sir_mean, param, TImax);

        Rm(ii).yGc(Rm(ii).yGc<0)=0;
        measured = [ti;sir_mean]';
        Rm(ii).measured = measured;
        maxM=max(sir_mean);
    end

    GKM      = [Rm(ii).xGKM',Rm(ii).yGKM'];
    Gamma    = [Rm(ii).xGm',Rm(ii).yGm];
    GammaF   = [Rm(ii).xGmF',Rm(ii).yGmF];
    Gauss    = [Rm(ii).xGs',Rm(ii).yGs];
    Combined = [Rm(ii).xGc',Rm(ii).yGc];
    CombinedComponents = [Rm(ii).fit_bc.gauss, Rm(ii).fit_bc.gamma];
    CombinedFracions = 100*[Rm(ii).fGs;Rm(ii).fGm];
    Ratio = Rm(ii).fit_bc.PHGauss/Rm(ii).fit_bc.PHGamma;
    GOF=[Rm(ii).gofGKM.rsquare,Rm(ii).gofGm.rsquare,...
    Rm(ii).gofGmF.rsquare,Rm(ii).gofGs.rsquare,Rm(ii).gofC.rsquare];
        
    % plot
    if isfield(SIR,'pp') && param.perVoxelCalc
       aslFitPlot(measured, GKM,Gamma,GammaF,Gauss, Combined,...
            CombinedComponents, GOF, Ratio, param, 'average',ax(4));
    else
       aslFitPlot(measured, GKM,Gamma,GammaF,Gauss, Combined,...
            CombinedComponents, GOF, Ratio, param, 'average',ax(3));
    end


    [TTP, FWHM, PH, AUC] = aslFitParams(Rm(ii).xGm,Rm(ii).yGm);
    [TTPF, FWHMF, PHF, AUCF] = aslFitParams(Rm(ii).xGmF,Rm(ii).yGmF);
    [TTPG, FWHMG, PHG, AUCG] = aslFitParams(Rm(ii).xGs,Rm(ii).yGs);
    [bTTP, bFWHM, bPH, bAUC] = aslFitParams(Rm(ii).xGc, Rm(ii).yGc);

        ResultsM(ii).ID                    =   Subject.ID;
        ResultsM(ii).Group                 =   Subject.study;
        ResultsM(ii).ROI                   =   ROI.name{ii};
        ResultsM(ii).Max                   =   maxM*100;
        % GKM
        ResultsM(ii).gkmT1                 =   Rm(ii).fitGKM(1);
        ResultsM(ii).Delta_t               =   Rm(ii).fitGKM(2);
        ResultsM(ii).tau                   =   Rm(ii).fitGKM(3);
        ResultsM(ii).PerfusionGKM          =   Rm(ii).fitGKM(4)*PerfusionK;
        % Gamma-variate
        ResultsM(ii).deltaTau              =   Rm(ii).fitGm.dt;
        ResultsM(ii).TTPgm                 =   TTP;
        ResultsM(ii).FWHMgm                =   FWHM;
        ResultsM(ii).PHgm                  =   PH*100;
        ResultsM(ii).AUCgm                 =   AUC*0.1;
        ResultsM(ii).T1gm                  =   Rm(ii).fitGm.T1;
        ResultsM(ii).PerfusionGamma        =   Rm(ii).fitGm.perfC*PerfusionK;
        % Gamma-variate fixed T1
        ResultsM(ii).deltaTauF             =   Rm(ii).fitGmF.dt;
        ResultsM(ii).TTPgmF                =   TTPF;
        ResultsM(ii).FWHMgmF               =   FWHMF;
        ResultsM(ii).PHgmF                 =   PHF*100;
        ResultsM(ii).AUCgmF                =   AUCF*0.1;
        ResultsM(ii).T1gmF                 =   param.T1;
        ResultsM(ii).PerfusionGammaFixedT1 =   Rm(ii).fitGmF.perfC*PerfusionK;
        % Gaussian
        ResultsM(ii).deltaT                =   Rm(ii).fitGs.b-2*Rm(ii).fitGs.c;
        ResultsM(ii).TTPgs                 =   TTPG;
        ResultsM(ii).FWHMgs                =   FWHMG;
        ResultsM(ii).PHgs                  =   PHG*100;
        ResultsM(ii).AUCgs                 =   AUCG*0.1;
        % bi-component
        ResultsM(ii).bi_deltaTGauss        =   Rm(ii).fit_bc.deltaTGauss;
        ResultsM(ii).bi_deltaTGamma        =   Rm(ii).fit_bc.deltaTGamma;
        ResultsM(ii).fGauss                =   CombinedFracions(1);
        ResultsM(ii).fGamma                =   CombinedFracions(2);
        ResultsM(ii).biRatio               =   Rm(ii).fit_bc.ratio;
        ResultsM(ii).biPerfusion           =   Rm(ii).fit_bc.perfC*PerfusionK;
        ResultsM(ii).biTTP                 =   bTTP;
        ResultsM(ii).biTTPGauss            =   Rm(ii).fit_bc.TTPGauss;
        ResultsM(ii).biTTPGamma            =   Rm(ii).fit_bc.TTPGamma;
        ResultsM(ii).biFWHM                =   bFWHM;
        ResultsM(ii).biPH                  =   bPH*100;
        ResultsM(ii).biPHGauss             =   Rm(ii).fit_bc.PHGauss*100;
        ResultsM(ii).biPHGamma             =   Rm(ii).fit_bc.PHGamma*100;
        ResultsM(ii).biAUC                 =   bAUC*0.1;
        % Goodness of fit
        ResultsM(ii).R2_GKM                =   Rm(ii).gofGKM.rsquare;
        ResultsM(ii).RMSE_GKM              =   Rm(ii).gofGKM.rmse;
        ResultsM(ii).SSE_GKM               =   Rm(ii).gofGKM.ssr;
        ResultsM(ii).R2_GV                 =   Rm(ii).gofGm.rsquare;
        ResultsM(ii).RMSE_GV               =   Rm(ii).gofGm.rmse;
        ResultsM(ii).SSE_GV                =   Rm(ii).gofGm.sse;
        ResultsM(ii).R2_GVfT1              =   Rm(ii).gofGmF.rsquare;
        ResultsM(ii).RMSE_GVfT1            =   Rm(ii).gofGmF.rmse;
        ResultsM(ii).SSE_GVfT1             =   Rm(ii).gofGmF.sse;
        ResultsM(ii).R2_Gauss              =   Rm(ii).gofGs.rsquare;
        ResultsM(ii).RMSE_Gauss            =   Rm(ii).gofGs.rmse;
        ResultsM(ii).SSE_Gauss             =   Rm(ii).gofGs.sse;
        ResultsM(ii).biR2                  =   Rm(ii).gofC.rsquare;
        ResultsM(ii).biRMSE                =   Rm(ii).gofC.rmse;
        ResultsM(ii).biSSE                 =   Rm(ii).gofC.sse;
        % actual values and fits x (TI) and y (SIR)
        ResultsM(ii).ti                    =   ti(:); 
        ResultsM(ii).y                     =   sir_mean(:);
        ResultsM(ii).y_filtered            =   ym(:);
        ResultsM(ii).x_fit_GKM             =   Rm(ii).xGKM(:);
        ResultsM(ii).y_fit_GKM             =   Rm(ii).yGKM(:);
        ResultsM(ii).x_fit_GV              =   Rm(ii).xGm(:);
        ResultsM(ii).y_fit_GV              =   Rm(ii).yGm(:);
        ResultsM(ii).x_fit_GVfT1           =   Rm(ii).xGmF(:);
        ResultsM(ii).y_fit_GVfT1           =   Rm(ii).yGmF(:);
        ResultsM(ii).x_fit_Gauss           =   Rm(ii).xGs(:);
        ResultsM(ii).y_fit_Gauss           =   Rm(ii).yGs(:);
        ResultsM(ii).x_fit_bi              =   Rm(ii).xGc(:);
        ResultsM(ii).y_fit_bi              =   Rm(ii).yGc(:);

    % Generate the filename
    filename = sprintf('results_s%02d_%s.pdf', ROI.slice_number, ROI.name{ii});


    % Add an overall title to the entire tiled layout
    if param.latex
        titleHandle = title(tiledLayoutHandle, strcat('$ \mathrm{',...
            Subject.ID, '} - \mathrm{', Subject.study, '}$'), 'Interpreter',...
            'latex');
    else
        titleHandle = title(tiledLayoutHandle, strcat(Subject.ID, ' - ',...
            Subject.study),'Interpreter','none');
    end
    
    
    % Optionally, you can set additional properties for the title
    titleHandle.FontSize = 20;
    titleHandle.Color = 'k';
    % Export the current figure or specific axes as a PDF
    exportgraphics(gcf, filename, 'ContentType', 'vector','Resolution', 300);

    if strcmp(param.visible, 'off')
        close(fgH);
    end
    delete(fgH)

end
    
varargout{1}=ResultsM;

if param.perVoxelCalc
    varargout{2}=ResultsP;
end

end