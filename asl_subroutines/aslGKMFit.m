function [paramsGKM, gof, x_fine, y_fit_fine] = aslGKMFit(timeInfo, sir, param, varargin)
%--------------------------------------------------------------------------
%   Function to fit ASL data to the GKM (General Kinetic Model) using 
%   numerical fitting with lsqcurvefit and flexible parameter handling 
%   through varargin parsing.
%
%   INPUT: 
%       timeInfo         - Time vector (1D array)                  (double)
%       sir              - Signal intensity ratio (1D array)       (double)
%       param            - Struct with default parameter values    (struct)
%       varargin         - Optional inputs                         (mixed)
%
%   OUTPUT:
%       paramsGKM        - Fitted parameters                       (double)
%       gof              - Goodness of fit (R^2, RMSE, SSR)        (struct)
%       x_fine           - Fine time grid                          (double)
%       y_fit_fine       - Fitted values on the fine time grid     (double)
%
%   Optional Inputs (varargin parsing):
%       maxX: Override maximum X value for fine grid generation    (double)
%       '-T1delta': Define a range for T1 bounds                   (double)
%__________________________________________________________________________
% vmalis@ucsd.edu
%--------------------------------------------------------------------------

 % Add zero points for the fit
    timeInfo = cat(2, 0, timeInfo);
    sir = cat(2, 0, sir);

    % Default values for optional inputs
    maxX = max(timeInfo);  % Default maxX value for the fine time grid
    T1 = param.T1;         % Default T1 from param
    T1delta = 0;           % Default T1delta (no bounds modification)

    timeInfo = cat(2, timeInfo,maxX+200);
    sir = cat(2, sir,0);


    % Parse optional inputs
    i = 1;
    while i <= length(varargin)
        if isnumeric(varargin{i})
            maxX = varargin{i};  % Override maxX if provided
        elseif ischar(varargin{i}) && strcmp(varargin{i}, '-T1delta') && i+1 <= length(varargin)
            T1delta = varargin{i+1};  % Set T1delta if provided
            i = i + 1;  % Skip the next value as it is part of '-T1delta'
        end
        i = i + 1;
    end

    % Define bounds for T1 based on T1delta, ensuring T1_lower >= 0
    T1_lower = max(0, T1 - T1delta);
    T1_upper = T1 + T1delta;

    % Define the perfusion [ml/min/ml]
    f = 100 / 60000;

    % Fixed pre-arrival delay and labeling duration
    [~,idx]=max(sir);
    Delta_t = max(timeInfo(idx)-(400 + (500 - 200) * rand),20);
    tau = 200 + (1000 - 200) * rand;

    % Initial guesses for the fit parameters
    initial_guess = [T1, Delta_t, tau, f];

    % Define bounds for the fit parameters
    % Define bounds for the fit parameters
    lb = [T1_lower, 100, 200, 5e-6];  
    ub = [T1_upper, 600, 800, 1e-4];   

    % Optimization options (default)
    %options = optimoptions('lsqcurvefit', 'Display', 'off');

    options = optimoptions('lsqcurvefit', ...
    'MaxIterations', 1000, ... % Increase from default (default is 400)
    'MaxFunctionEvaluations', 1000, ... % Allow more function calls
    'TolFun', 1e-12, ... % Reduce function tolerance for better convergence
    'TolX', 1e-12, ...
    'Algorithm', 'interior-point',...
    'Display', 'off'); % Reduce step tolerance)

    % Perform the fit using lsqcurvefit

    
    try
        [paramsGKM, ~] = lsqcurvefit(@gkm_numeric, initial_guess, timeInfo, sir, lb, ub, options);
    catch
        warning('lsqcurvefit failed — returning zeros for paramsGKM and output');
        paramsGKM = zeros(1, numel(initial_guess));
        gof = struct('rmse', NaN, 'rsquare', NaN, 'ssr', NaN);
        x_fine = linspace(min(timeInfo), maxX, 1000);
        y_fit_fine = zeros(size(x_fine));
    return;
    end
    
    % Calculate fitted values on the given time grid
    y_fit = gkm_numeric(paramsGKM, timeInfo);
    
    % Calculate residuals and Goodness of Fit
    residuals = sir - y_fit;
    RMSE = sqrt(mean(residuals.^2));
    SSR = sum(residuals.^2); % Sum of squared residuals
    SST = sum((sir - mean(sir)).^2); % Total sum of squares
    R_squared = 1 - (SSR / SST);
    
    % Store goodness of fit metrics in a struct
    gof = struct('rmse', RMSE, 'rsquare', R_squared, 'ssr', SSR);
    
    % Generate fine grid and calculate fitted values
    x_fine = linspace(min(timeInfo), maxX, 1000);  % Use maxX for fine grid
    y_fit_fine = gkm_numeric(paramsGKM, x_fine);

end