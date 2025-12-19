function [fit_result, gof, x_fine, y_fine, fracGauss, fracGamma, fit_bc] = aslCombinedFit2(x, y, param, varargin)
    % Add zero points for the fit
    x = cat(2, 0, x);
    y = cat(2, 0, y);

        % Add zero points for the fit
        x = cat(2, 0, x);
        y = cat(2, 0, y);
        
        % Ensure at least 6 data points
        if numel(x) < 6
            x = zeros(1, 6);
            y = zeros(1, 6);
        end

    % --- Existing fitting code below this line ---

    % Default values
    maxX = max(x);

    % Parse optional inputs
    i = 1;
    while i <= length(varargin)
        if isnumeric(varargin{i})
            maxX = varargin{i};
        end
        i = i + 1;
    end

    % Extract T1 value from param structure
    T1 = param.T1;

    % Define the fitting functions for Gaussian and Gamma
    gaussian = @(a, b, c, t) a * exp(-((t - b).^2 / (2 * c^2)));  % Positive Gaussian
    gamma_variate = @(perfC, shift, TI) perfC .* (TI - shift) .* exp(-(TI - shift) ./ T1) .* ((TI - shift) > 0); % Flexible Gamma Variate with shift

    % Define the combined fitting function with constraints
    combinedModel = @(fracGauss, a, b, c, perfC, shift, t) ...
        fracGauss * gaussian(a, b, c, t) + (1 - fracGauss) * gamma_variate(perfC, shift, t);

    % Set initial guesses and bounds
    initial_guesses_combined = [0.5, max(y), mean(x), std(x), max(y),0];
    lower_bounds_combined = [0, 0, min(x), 0, 0, 0];
    upper_bounds_combined = [1, max(y), max(x), Inf, Inf, max(x)];

    % Define the fit type
    fit_type_combined = fittype(@(fracGauss, a, b, c, perfC, shift, t) combinedModel(fracGauss, a, b, c, perfC, shift, t), ...
        'coefficients', {'fracGauss', 'a', 'b', 'c', 'perfC', 'shift'}, ...
        'independent', 't');

    % Set the fit options
    options_combined = fitoptions('Method', 'NonlinearLeastSquares', ...
        'StartPoint', initial_guesses_combined, ...
        'Lower', lower_bounds_combined, ...
        'Upper', upper_bounds_combined);

    % Perform the fit
    [fit_result, gof] = fit(x', y', fit_type_combined, options_combined);

    % Extract the fraction of Gaussian component from the fit result
    fracGauss = fit_result.fracGauss;
    fracGamma = 1 - fracGauss;

    % Generate fine grid for fitted values
    x_fine = linspace(min(x), maxX, 1000);
    y_fine = feval(fit_result, x_fine);

    % Remove negative values
    y_fine(y_fine < 0) = 0;

    % Extract the fit parameters for Gaussian and Gamma
    a = fit_result.a;
    b = fit_result.b;
    c = fit_result.c;
    perfC = fit_result.perfC;
    shift = fit_result.shift;

    % Compute Gaussian and Gamma components separately
    gaussian_y = gaussian(a, b, c, x_fine);
    gamma_y = gamma_variate(perfC, shift, x_fine);

    % Find x-coordinate of maximum for Gaussian
    [max_gaussian_y, max_gaussian_idx] = max(gaussian_y);
    x_max_gaussian = x_fine(max_gaussian_idx);
    fit_bc.PHGauss = max_gaussian_y * fracGauss;
    fit_bc.TTPGauss = x_max_gaussian;

    % Find x-coordinate of maximum for Gamma
    [max_gamma_y, max_gamma_idx] = max(gamma_y);
    x_max_gamma = x_fine(max_gamma_idx);
    fit_bc.PHGamma = max_gamma_y * fracGamma;
    fit_bc.TTPGamma = x_max_gamma;

    % Ratio of Gaussian to Gamma
    fit_bc.ratio = (max_gaussian_y * fracGauss) / (max_gamma_y * fracGamma);
end