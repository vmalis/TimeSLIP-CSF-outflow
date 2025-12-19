function A_out = ASL_imReg(A, param)
%--------------------------------------------------------------------------
%
%   ASL_imReg
%
%   Subroutine to perform rigid registration between ASL M0, ON, and CTRL
%   image volumes. Supports both 2D and 3D acquisitions with alternate and
%   non-alternate ASL labeling schemes. 
%
%   INPUT:  A       Struct containing:
%                      - On:   3D array of ON images ('double')
%                      - Ctrl: 3D array of Control images ('double')
%                      - M0:   3D array of M0 images ('double')
%           param   Struct containing:
%                      - acq_dim: '2d' or '3d' (acquisition type)
%                      - nonAlternate: boolean (true = non-alternate labeling)
%
%   OUTPUT: A_out   Struct with registered On, Ctrl (M0 unchanged)
%                   Also exports registration overlays as PDF files
%
%   NOTES:
%     - Rigid registration: 3 degrees of freedom (rotation + x/y translation)
%     - Similarity metric: Mattes Mutual Information
%     - Optimizer: One Plus One Evolutionary Optimizer
%     - Interpolation: Linear
%     - Overlay outputs: purple (before) vs green (after), black background
%
%__________________________________________________________________________
%  VM vmalis@ucsd.edu
%--------------------------------------------------------------------------

A_out = A;

if ~isfield(param, 'registration')
    error('param.registration must be specified as either "auto" or "all".');
end

if strcmpi(param.registration, 'all')
    fprintf('Registration = "all" — Register ALL ON and CTRL slices independently to M0\n');
    
    if strcmp(param.acq_dim,"2d")
        numSlices = size(A.On, 3);
    else
        numSlices = size(A.M0, 3);
    end
    
    A_out.On   = zeros(size(A.On),   'like', A.On);
    A_out.Ctrl = zeros(size(A.Ctrl), 'like', A.Ctrl);
    target_size = size(A.M0(:,:,1));
    
    for k = 1:numSlices
        tformON   = estimate_rigid(A.On(:,:,k),   A.M0(:,:,1));
        tformCTRL = estimate_rigid(A.Ctrl(:,:,k), A.M0(:,:,1));
        
        A_out.On(:,:,k)   = imwarp(A.On(:,:,k),   tformON,   'OutputView', imref2d(target_size), 'Interp', 'linear', 'FillValues', 0);
        A_out.Ctrl(:,:,k) = imwarp(A.Ctrl(:,:,k), tformCTRL, 'OutputView', imref2d(target_size), 'Interp', 'linear', 'FillValues', 0);
    end

elseif strcmpi(param.registration, 'extreme')
    fprintf('Registration = "extreme" — Non-rigid slice-by-slice registration to M0 using imregdemons\n');
    figure; imshow(mat2gray(A.On(:,:,1))); title('Draw box around the heart (moving region)');
    h = drawrectangle('Color', 'r');
    wait(h);
    
    % Get bounding box coordinates
    rect = round(h.Position);  % [x, y, width, height]
    c1 = rect(1); r1 = rect(2);
    c2 = c1 + rect(3) - 1;
    r2 = r1 + rect(4) - 1;
    bbox = struct('r1', r1, 'r2', r2, 'c1', c1, 'c2', c2);
    close(gcf);

    fixed_crop = mat2gray(A.M0(bbox.r1:bbox.r2, bbox.c1:bbox.c2, 1));  % Normalized for demons
    
    numSlices = size(A.On, 3);
    A_out = A;
    
    for k = 1:numSlices
        % --- Get original crops (preserve intensity)
        orig_on_crop   = A.On(bbox.r1:bbox.r2, bbox.c1:bbox.c2, k);
        orig_ctrl_crop = A.Ctrl(bbox.r1:bbox.r2, bbox.c1:bbox.c2, k);
    
        % --- Also get normalized versions for demons registration
        moving_on_norm   = mat2gray(orig_on_crop);
        moving_ctrl_norm = mat2gray(orig_ctrl_crop);
    
        % --- Run demons registration on normalized images
        [dispON, ~] = imregdemons(moving_on_norm, fixed_crop, 800, ...
            'AccumulatedFieldSmoothing', 0.9, 'PyramidLevels', 3,'DisplayWaitbar', false);
        [dispCTRL, ~] = imregdemons(moving_ctrl_norm, fixed_crop, 800, ...
            'AccumulatedFieldSmoothing', 0.9, 'PyramidLevels', 3,'DisplayWaitbar', false);
    
        % --- Warp the ORIGINAL unnormalized image using same deformation
        A_out.On(bbox.r1:bbox.r2, bbox.c1:bbox.c2, k) = imwarp(orig_on_crop, dispON, 'linear');
        A_out.Ctrl(bbox.r1:bbox.r2, bbox.c1:bbox.c2, k) = imwarp(orig_ctrl_crop, dispCTRL, 'linear');
    
    end


elseif strcmpi(param.registration, 'auto')
    if strcmpi(param.acq_dim, '2d')
        fprintf('Registration = "auto", Case: 2D — Single transform for full volumes\n');
        
        tformON   = estimate_rigid(A.On(:,:,end),   A.M0(:,:,1));
        tformCTRL = estimate_rigid(A.Ctrl(:,:,end), A.M0(:,:,1));
        
        A_out.On   = apply_tform_3d(A.On,   tformON,   size(A.M0(:,:,1)));
        A_out.Ctrl = apply_tform_3d(A.Ctrl, tformCTRL, size(A.M0(:,:,1)));

    elseif strcmpi(param.acq_dim, '3d') && param.nonAlternate
        fprintf('Registration = "auto", Case: 3D — Independent slice-by-slice registration to M0\n');
        
        numSlices = size(A.M0, 3);
        A_out.On   = zeros(size(A.On),   'like', A.On);
        A_out.Ctrl = zeros(size(A.Ctrl), 'like', A.Ctrl);
        target_size = size(A.M0(:,:,1));
        
        for k = 1:numSlices
            tformON   = estimate_rigid(A.On(:,:,k),   A.M0(:,:,1));
            tformCTRL = estimate_rigid(A.Ctrl(:,:,k), A.M0(:,:,1));
            
            A_out.On(:,:,k)   = imwarp(A.On(:,:,k),   tformON,   'OutputView', imref2d(target_size), 'Interp', 'linear', 'FillValues', 0);
            A_out.Ctrl(:,:,k) = imwarp(A.Ctrl(:,:,k), tformCTRL, 'OutputView', imref2d(target_size), 'Interp', 'linear', 'FillValues', 0);
        end

    elseif strcmpi(param.acq_dim, '3d') && ~param.nonAlternate
        fprintf('Registration = "auto", Case: 3D — Register Ctrl slice, apply to same On slice\n');
        
        numSlices = size(A.M0, 3);
        A_out.On   = zeros(size(A.On),   'like', A.On);
        A_out.Ctrl = zeros(size(A.Ctrl), 'like', A.Ctrl);
        target_size = size(A.M0(:,:,1));
        
        for k = 1:numSlices
            tformCTRL = estimate_rigid(A.Ctrl(:,:,k), A.M0(:,:,1));
            
            A_out.Ctrl(:,:,k) = imwarp(A.Ctrl(:,:,k), tformCTRL, 'OutputView', imref2d(target_size), 'Interp', 'linear', 'FillValues', 0);
            A_out.On(:,:,k)   = imwarp(A.On(:,:,k),   tformCTRL, 'OutputView', imref2d(target_size), 'Interp', 'linear', 'FillValues', 0);
        end

    else
        error('Unsupported parameter combination for "auto" registration.');
    end

else
    error('Invalid value for param.registration. Must be "auto" or "all".');
end

% Export overlays
export_overlay_pdf(A.On,   A_out.On,   'On_Registration_Result.pdf');
export_overlay_pdf(A.Ctrl, A_out.Ctrl, 'Ctrl_Registration_Result.pdf');

end

% -------------------------------------------------------
% Helper: Estimate rigid transformation
function tform = estimate_rigid(moving, fixed)
    [optimizer, metric] = imregconfig('multimodal');
    tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
end

% -------------------------------------------------------
% Helper: Apply transformation to 3D array
function vol_out = apply_tform_3d(vol, tform, target_size)
    numSlices = size(vol, 3);
    vol_out = zeros([target_size, numSlices], 'like', vol);
    outputView = imref2d(target_size);
    for k = 1:numSlices
        vol_out(:,:,k) = imwarp(vol(:,:,k), tform, 'OutputView', outputView, 'Interp', 'linear', 'FillValues', 0);
    end
end

% -------------------------------------------------------
% Helper: Export overlays to PDF
function export_overlay_pdf(before_vol, after_vol, filename)
    fig = figure('Units', 'normalized', 'OuterPosition', [0 0 1 0.3], 'Color', 'k'); % Black figure background
    t = tiledlayout(1, size(before_vol,3), 'TileSpacing', 'none', 'Padding', 'none');
    
    for k = 1:size(before_vol, 3)
        ax = nexttile;
        
        before = mat2gray(before_vol(:,:,k));
        after = mat2gray(after_vol(:,:,k));
        
        % Create RGB overlay
        overlay = cat(3, ...
            before, ... % Red channel: BEFORE (Original)
            after,  ... % Green channel: AFTER (Registered)
            before);    % Blue channel: BEFORE (Original)

        imshow(overlay);
        title(sprintf('Slice %d', k), 'Color', 'w');
        ax.XColor = 'none'; 
        ax.YColor = 'none';
        ax.Color = 'k'; % Black axes background
    end
    
    % Export using export_fig
    if exist('export_fig', 'file') ~= 2
        error('export_fig not found. Please download it from the MATLAB File Exchange.');
    end

    export_fig(filename, '-pdf', '-opengl', '-silent'); % <- Your requested options
    
    close(fig);
end