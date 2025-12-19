function varargout = aslFixedRO(SIR, ti, Subject, ROI, param, Image)
%--------------------------------------------------------------------------
%
%   Function to process ASL data and generate plots
%   This function generates results and plots for ASL data analysis.
%   It uses ROI and SIR data to generate results for mean and per-pixel data.
%
%--------------------------------------------------------------------------

% Adjust parameters as needed 0.1 (default) strong filter, 0.5 weak
filter=0.1;

ResultsM = struct();    % results for mean
ResultsP = struct();    % results for per pixel

% Default max TI [ms] for plots
TImax = max(ti);

% Crop and resize the image based on ROI
I = cropAndResizeImage(mat2gray(Image), ROI);

mkdir('plots');
cd('plots')

for ii = 1:size(ROI.maskIm, 3)

    fgH = figure('Position', [100, 100, 1200, 400], 'Visible', 'off');
    set(fgH, 'Color', 'w');
    set(fgH, 'Visible', param.visible);
    
    % Create a 2x5 tiled layout
    tiledLayoutHandle = tiledlayout(2, 5);

    % Initialize handles array
    ax = gobjects(1, 3);
    
    if min(size(SIR.pp)) == 1
        sir_mean = SIR.mean';
    else
        sir_mean = SIR.mean(ii, :);   
    end

    % Apply low pass filter
    if strcmpi(param.filter_lowpass,'lowPass')
        ym = aslLowPass(sir_mean,max(ti),filter);
    else
        ym = sir_mean;
    end

    % Show ROI in HighRes Square Crop
    mask = cropAndResizeImage(ROI.maskIm(:,:,ii), ROI);
    mask(isnan(mask))=0;
    mask=logical(mask);

    % Define main axes layout
    ax(1) = nexttile(tiledLayoutHandle,[2,2]);
    axis(ax(1), 'square');

    tiledLayoutHandle.TileSpacing = 'tight';
    tiledLayoutHandle.Padding = 'tight';
    displayROIImage(I, mask, ROI.name{ii}, 0, [0, 1, 0], '', param, ax(1));

    % Check for per voxel calculation
    if isfield(SIR, 'pp') && param.perVoxelCalc
        
        if min(size(SIR.pp)) == 1
            sir_pp = SIR.pp';
        else
            sir_pp = SIR.pp(ii, :);   
        end

        % Apply low pass filter
        if strcmpi(param.filter_lowpass,'lowPass')
            ypp = aslLowPass(sir_pp,max(ti),filter);
        else
            ypp = sir_pp;
        end

        ax(2) = nexttile(tiledLayoutHandle, [1, 3]); % Per-pixel plot in the first row
        % Plot ECG wave for per-pixel data
        set(fgH, 'Visible', param.visible);
        [cycle_num, phase, percent_in_phase, max_time, p_rr, max_SIR] = ...
            plot_ecg_wave(param.bpm, ypp * 100, ti, 0, ax(2), param,'per pixel');


        ResultsP(ii).ID                 = Subject.ID;
        ResultsP(ii).Group              = Subject.study;
        ResultsP(ii).ROI                = ROI.name{ii};
        ResultsP(ii).max_SIR            = max_SIR;
        ResultsP(ii).cycle_num          = cycle_num;
        ResultsP(ii).phase              = phase;
        ResultsP(ii).percent_in_phase   = percent_in_phase;
        ResultsP(ii).max_time           = max_time;
        ResultsP(ii).percent_rr         = p_rr;

        ax(3) = nexttile(tiledLayoutHandle, [1, 3]); % Mean plot occupies the first row

    else

        ax(3) = nexttile(tiledLayoutHandle, [2, 3]); % Mean plot occupies the first row

    end

    % Plot ECG wave for mean data
    [cycle_num, phase, percent_in_phase, max_time, p_rr, max_SIR] = ...
            plot_ecg_wave(param.bpm, ym * 100, ti, 0, ax(3), param,'mean');


    ResultsM(ii).ID                 = Subject.ID;
    ResultsM(ii).Group              = Subject.study;
    ResultsM(ii).ROI                = ROI.name{ii};
    ResultsM(ii).max_SIR            = max_SIR;
    ResultsM(ii).cycle_num          = cycle_num;
    ResultsM(ii).phase              = phase;
    ResultsM(ii).percent_in_phase   = percent_in_phase;
    ResultsM(ii).max_time           = max_time;
    ResultsM(ii).percent_rr         = p_rr;
   

    % Generate the filename
    filename = sprintf('results_s%02d_%s.pdf', ROI.slice_number, ROI.name{ii});

    % Add an overall title to the entire tiled layout
    if param.latex
        titleHandle = title(tiledLayoutHandle, strcat('$ \mathrm{', Subject.ID, '} - \mathrm{', Subject.study,'\, HR=', num2str(param.bpm), 'bpm', '}$'), 'Interpreter', 'latex');
    else
        titleHandle = title(tiledLayoutHandle, strcat(Subject.ID, ' - ', Subject.study, ' HR=',num2str(param.bpm), 'bpm'), 'Interpreter', 'none');
    end

    % Optionally, you can set additional properties for the title
    titleHandle.FontSize = 20;
    titleHandle.Color = 'k';

    % Export the current figure or specific axes as a PDF
    exportgraphics(fgH, filename, 'ContentType', 'vector', 'Resolution', 300);

    if strcmp(param.visible, 'off')
        close(fgH);
    end
    delete(fgH);

end

varargout{1} = ResultsM;

if param.perVoxelCalc
    varargout{2} = ResultsP;
end

end