function param = tslipParamCheck(param)
%--------------------------------------------------------------------------
%
%   Subroutine to check tSLIP parameters
%
%   INPUT:  
%       param - structure with fields:
%           acq_type                  - 'string'  ('bright', 'dark')
%           acq_dim                   - 'string'  ('2d', '3d')
%           TI                        - 'array'   (2 positive integers, only if acq_dim = '2d')
%           dicomScale                - 'boolean'
%           perVoxelCalc              - 'boolean'
%           fitplot                   -  cell array of {'gm', 'gs', 'bi'}
%           filter_median3D           - 'boolean'
%           filter_median2D           - 'boolean'
%           filter_sd                 - 'boolean'
%           filter_lowpass            - 'boolean'
%           averageSlices             - 'boolean'
%           dataType                  - 'string'  ('regular', 'other')
%           latex                     - 'boolean'
%           visible                   - 'boolean'
%           allslices                 - 'boolean'
%           subtractionMap            - 'logical' (2-element array)
%           paleteview                - 'logical' (2-element array)
%           sir2Dcmap                 - 'logical' (2-element array)
%           sir2DcmapFull             - 'boolean'
%           montageALL                - 'logical' (2-element array)
%           contrastNorm              - 'boolean'
%           fixedRO                   - 'boolean' (default: false)
%           bpm                       - 'integer' (30-200, only checked if fixedRO is true)
%           SIRLim                    - 'integer' (> 0)
%           SIRClim                   - 'integer' (> 0)
%           advancedGKM               - 'boolean'
%
%   OUTPUT: 
%       Validated and possibly modified param structure
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

    % Check acq_type
    valid_acq_types = {'bright', 'dark'};
    if ~isfield(param, 'acq_type') || ~ismember(param.acq_type, valid_acq_types)
        error('param.acq_type must be either ''bright'' or ''dark''.');
    end
    
    % Check acq_dim
    valid_acq_dims = {'2d', '3d'};
    if ~isfield(param, 'acq_dim') || ~ismember(param.acq_dim, valid_acq_dims)
        error('param.acq_dim must be either ''2d'' or ''3d''.');
    end
    
    % Check TI based on acq_dim
    if strcmp(param.acq_dim, '2d')
        if ~isfield(param, 'TI') || ~isnumeric(param.TI) || ...
                numel(param.TI) ~= 2 || any(param.TI < 0) || ...
                    any(mod(param.TI, 1) ~= 0)
            error('param.TI must be an array with 2 positive integers if acq_dim is ''2d''.');
        end
        param.filter_median3D = false;
        param.averageSlices = false;
    elseif strcmp(param.acq_dim, '3d')
        if isfield(param, 'TI')
            error('param.TI should not be provided if acq_dim is ''3d''.');
        end
    end

    % Check fitplot
    valid_fitplot_options = {'gkm', 'gm','gmf','gs', 'bi', 'biC'};
    if ~isfield(param, 'fitplot') || ~iscell(param.fitplot) || ...
            ~all(ismember(param.fitplot, valid_fitplot_options))
        error('param.fitplot must be a cell array containing any combination of ''gkm'',''gm'',''gmf'',''gs'' , ''bi'' &  ''biC''.');
    end


    valid_filter_options = {'lowPass' 'movinAver' 'off' };
    if ~isfield(param, 'fitplot') || ~iscell(param.fitplot) || ...
            ~all(ismember(param.filter_lowpass, valid_filter_options))
        error('param.filter_lowpass can be: ''lowPass'',''movinAver'' or ''off''.');
    end

    % Check boolean parameters
    bool_params = {'dicomScale', 'perVoxelCalc', 'filter_median3D', ...
        'filter_median2D', 'filter_sd', 'averageSlices', ...
        'latex', 'allslices', 'contrastNorm', ...
        'sir2DcmapFull', 'visible', 'fixedRO','advancedGKM'};
    
    for i = 1:length(bool_params)
        param_name = bool_params{i};
        if ~isfield(param, param_name) || ~islogical(param.(param_name))
            error('param.%s must be a boolean value.', param_name);
        end
    end

    % Check 2-element logical array parameters
    array_params = {'subtractionMap', 'paleteview', 'sir2Dcmap', 'montageALL'};
    for i = 1:length(array_params)
        param_name = array_params{i};
        if ~isfield(param, param_name) || ~islogical(param.(param_name)) || ...
                numel(param.(param_name)) ~= 2
            error('param.%s must be a 2-element logical array.', param_name);
        end
    end

    % Check dataType
    valid_data_types = {'regular', 'other'};
    if ~isfield(param, 'dataType') || ~ismember(param.dataType, valid_data_types)
        error('param.dataType must be either ''regular'' or ''other''.');
    end

    % Check bpm if fixedRO is true
    if param.fixedRO
        if ~isfield(param, 'bpm') || ~isnumeric(param.bpm) || ...
                ~isscalar(param.bpm) || param.bpm < 30 || param.bpm > 200
            error('param.bpm must be a scalar integer between 30 and 200 when param.fixedRO is true.');
        end
    end

    % Check SIRLim
    if ~isfield(param, 'SIRLim') || ~isnumeric(param.SIRLim) || ...
        param.SIRLim <= 0
    error('param.SIRLim must be a positive numeric value.');
    end

    % Check SIRClim
    if ~isfield(param, 'SIRClim') || ~isnumeric(param.SIRClim) || ...
            ~isscalar(param.SIRClim) || param.SIRClim <= 0 || ...
            mod(param.SIRClim, 1) ~= 0
        error('param.SIRClim must be a positive integer greater than 0.');
    end
end