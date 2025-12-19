function [fit_result, gof, x_fine, y_fine] = aslGaussianFit(x, y, maxX)
%--------------------------------------------------------------------------
%
%   Function to fit data to a Gaussian model for ASL
%
%   INPUT:  x values (1D array)                         'x'             (double)
%           y values (1D array)                         'y'             (double)
%           maximum x value for fine grid               'maxX'          (double)
%   OUTPUT: fit result                                  'fit_result'    (fit object)
%           goodness of fit                             'gof'           (struct)
%           fine grid x values                          'x_fine'        (double)
%           fitted y values on fine grid                'y_fine'        (double)
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

 % Append new element to the beginning of x and y
    x = cat(2, 0, x);
    y = cat(2, 0, y);

    % Constrained Gaussian function (starting from 0,0 and positive y)
    gaussian = @(a, b, c, t) a * (exp(-((t - b).^2 / (2 * c^2))) - exp(-((0 - b).^2 / (2 * c^2))));

    % Set the initial guesses for the parameters
    initial_guesses_gaussian = [max(y), mean(x), std(x)];
    % Set bounds for the parameters
    lower_bounds_gaussian = [0, min(x), 0];
    upper_bounds_gaussian = [max(y)*1.5, max(x), Inf];
    
    % Constrain the Gaussian to pass through (0,0)
    fit_type_gaussian = fittype(@(a, b, c, t) gaussian(a, b, c, t), ...
        'coefficients', {'a', 'b', 'c'}, ...
        'independent', 't');
    
    % Perform the fit with bounds
    options_gaussian = fitoptions('Method', 'NonlinearLeastSquares', ...
        'StartPoint', initial_guesses_gaussian, ...
        'Lower', lower_bounds_gaussian, ...
        'Upper', upper_bounds_gaussian);
    [fit_result, gof] = fit(x', y', fit_type_gaussian, options_gaussian);

    % Generate fine grid for fitted values
    x_fine = linspace(min(x), maxX, 1000);
    y_fine = feval(fit_result, x_fine);
    
    % Remove negative values
    % positive_indices = y_fine >= 0;
    % x_fine = x_fine(positive_indices);
    % y_fine = y_fine(positive_indices);

    negative_indices = y_fine < 0;
    y_fine(negative_indices) = 0;

end