function [fit_result, gof, x_fine, y_fine] = aslGammaFitCSF(x, y, varargin)
%--------------------------------------------------------------------------
%
%   Function to fit data to a gamma variate model for ASL
%
%   INPUT:  x values (1D array)                         'x'             (double)
%           y values (1D array)                         'y'             (double)
%           maximum x value for fine grid (optional)    'maxX'          (double)
%           T1 value (optional, specified by flag -T1)  'T1'            (double)
%
%   OUTPUT: fit result                                  'fit_result'    (fit object)
%           goodness of fit                             'gof'           (struct)
%           fine grid x values                          'x_fine'        (double)
%           fitted y values on fine grid                'y_fine'        (double)
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

    % Default values
    maxX = max(x);
    T1_fixed = false;
    T1_default = 4000; % Default value for T1 if not provided
    min_hShift = 400; % Minimum horizontal shift to enforce non-zero start
    hShift_initial = 400 + (900 - 400) * rand; % Random initial guess for hShift between 300 and 700

    % Parse optional inputs
    i = 1;
    while i <= length(varargin)
        if isnumeric(varargin{i})
            if i == 1 % Assume the first numeric input is maxX
                maxX = varargin{i};
            end
        elseif ischar(varargin{i}) && strcmp(varargin{i}, '-T1') && i+1 <= length(varargin)
            T1 = varargin{i+1};
            T1_fixed = true;
            i = i + 1; % Skip the next value as it is part of '-T1'
        end
        i = i + 1;
    end

    % Add zero points for the fit
    x = cat(2, 0, x);
    y = cat(2, 0, y);

    % Define the fitting function for Y(TI)
    gamma_variate = @(perfC, T1, hShift, TI) ...
        perfC .* max(0, (TI - hShift)) .* exp(-max(0, (TI - hShift)) ./ T1);

    % Calculate upper bound for perfC
    e = exp(1); % Euler's number
    max_y = max(y); % Maximum of the input data

    % Estimate initial hShift based on data (optional refinement)
    non_zero_idx = find(y > 0, 1, 'first'); % Find first non-zero y-value
    if ~isempty(non_zero_idx) && non_zero_idx > 1
        hShift_initial = min(max(x(non_zero_idx-1), min_hShift), 700); % Use x-value before rise, bounded
    end

    if T1_fixed
        % T1 is fixed
        perfC_upper = max_y * e / T1; % Upper bound for perfC
        fit_type = fittype(@(perfC, hShift, TI) gamma_variate(perfC, T1, hShift, TI), ...
            'coefficients', {'perfC', 'hShift'}, ...
            'independent', 'TI');
        initial_guess = [min(max_y, perfC_upper), hShift_initial]; % Initial guesses for perfC and hShift
        lower_bound = [0, min_hShift]; % Lower bounds for perfC and hShift
        upper_bound = [perfC_upper, 800]; % Upper bounds for perfC and hShift
        options = fitoptions('Method', 'NonlinearLeastSquares', ...
            'StartPoint', initial_guess, ...
            'Lower', lower_bound, ...
            'Upper', upper_bound);
    else
        % T1 is a variable parameter
        perfC_upper = max_y * e / T1_default; % Initial estimate for perfC bound
        fit_type = fittype(@(perfC, T1, hShift, TI) gamma_variate(perfC, T1, hShift, TI), ...
            'coefficients', {'perfC', 'T1', 'hShift'}, ...
            'independent', 'TI');
        initial_guesses = [min(max_y, perfC_upper), T1_default, hShift_initial]; % Initial guesses
        lower_bounds = [0, 1, min_hShift]; % Lower bounds for perfC, T1, and hShift
        upper_bounds = [perfC_upper * 2, Inf, 800]; % Relaxed upper bound for perfC
        options = fitoptions('Method', 'NonlinearLeastSquares', ...
            'StartPoint', initial_guesses, ...
            'Lower', lower_bounds, ...
            'Upper', upper_bounds);
    end

    % ======== FIT WITH FALLBACK TO ZERO-PARAMETER MODEL ========
    try
        % Perform the fit (original behaviour)
        [fit_result, gof] = fit(x', y', fit_type, options);
        
        % Generate fine grid for fitted values
        x_fine = linspace(min(x), maxX, 1000); % Use maxX for fine grid generation
        y_fine = feval(fit_result, x_fine);
        
    catch ME
        % If fit fails (e.g., bounds too close), return a "null" gamma fit
        warning('aslGammaFitCSF:FitFailed', ...
                'Fit failed (%s). Returning zero-parameter fit.', ME.message);

        % Fine grid consistent with original
        x_fine = linspace(min(x), maxX, 1000);
        y_fine = zeros(size(x_fine));
        y_fine = y_fine';
        % Fake fit_result with required fields
        fit_result = struct();
        fit_result.perfC  = 0;
        fit_result.hShift = 0;
        if T1_fixed
            fit_result.T1 = T1;   % keep fixed T1 value if provided
        else
            fit_result.T1 = 0;
        end

        % Build gof struct with same fields as MATLAB fit()
        nData    = numel(x);              % number of data points
        if T1_fixed
            nParams = 2;                  % perfC, hShift
        else
            nParams = 3;                  % perfC, T1, hShift
        end
        dfe = max(nData - nParams, 0);    % degrees of freedom (non-negative)

        gof = struct();
        gof.sse        = NaN;       % no meaningful SSE from a failed fit
        gof.rsquare    = 0;         % no explained variance
        gof.dfe        = dfe;
        gof.adjrsquare = 0;
        gof.rmse       = NaN;

        return;
    end
    % ===========================================================

    % Check if fitted values exceed max(y) and warn
    if max(y_fine) > max_y
        warning('Fitted values exceed maximum input data (%.2f > %.2f). Consider adjusting constraints.', max(y_fine), max_y);
    end

    % Check if hShift is too close to zero and warn
    if T1_fixed
        hShift_fitted = fit_result.hShift;
    else
        hShift_fitted = fit_result.hShift;
    end
    if hShift_fitted < min_hShift * 1.1
        warning('Fitted hShift (%.2f) is close to minimum bound (%.2f). Consider adjusting initial guess or bounds.', hShift_fitted, min_hShift);
    end

end