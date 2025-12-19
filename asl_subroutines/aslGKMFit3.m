function [paramsGKM, gof, x_fine, y_fit_fine] = aslGKMFit3(timeInfo, sir, param, varargin)
%--------------------------------------------------------------------------
%   Function to fit ASL data to the GKM (General Kinetic Model) using
%   global optimization (GlobalSearch) to find the best parameters 
%   that maximize R².
%
%   INPUT: 
%       timeInfo         - Time vector (1D array)                  (double)
%       sir              - Signal intensity ratio (1D array)       (double)
%       param            - Struct with default parameter values    (struct)
%       varargin         - Optional inputs                         (mixed)
%
%   OUTPUT:
%       paramsGKM        - Best fitted parameters                  (double)
%       gof              - Goodness of fit (R², RMSE, SSR)         (struct)
%       x_fine           - Fine time grid                          (double)
%       y_fit_fine       - Fitted values on the fine time grid     (double)
%__________________________________________________________________________
% vmalis@ucsd.edu
%--------------------------------------------------------------------------

    % Add zero points for the fit
    timeInfo = cat(2, 0, timeInfo);
    sir = cat(2, 0, sir);

    % Default values for optional inputs
    maxX = max(timeInfo);  
    T1 = param.T1;         
    T1delta = 0;           

    % Parse optional inputs
    i = 1;
    while i <= length(varargin)
        if isnumeric(varargin{i})
            maxX = varargin{i};  
        elseif ischar(varargin{i}) && strcmp(varargin{i}, '-T1delta') && i+1 <= length(varargin)
            T1delta = varargin{i+1};  
            i = i + 1;  
        end
        i = i + 1;
    end

    % Define bounds for T1 based on T1delta, ensuring T1_lower >= 0
    T1_lower = max(0, T1 - T1delta);
    T1_upper = T1 + T1delta;

    % Define the perfusion [ml/min/ml]
    f = 100*rand(1,1) / 60000;

    % Fixed pre-arrival delay and labeling duration
    [~,idx] = max(sir);
    Delta_t = max(timeInfo(idx) - (400 + (500 - 200) * rand), 20);
    tau = 700 + (500 - 300) * rand;

    % Initial guess for the fit parameters
    initial_guess = [T1, Delta_t, tau, f*0.6];

    % Define bounds for the fit parameters
    lb = [T1_lower, 100, 500, 5e-6];  
    ub = [T1_upper, 600, 1200, 2e-4];   

    % Optimization options for fmincon
    options = optimoptions('fmincon', ...
        'MaxIterations', 5000, ...
        'MaxFunctionEvaluations', 5000, ...
        'TolFun', 1e-16, ...
        'TolX', 1e-16, ...
        'Algorithm', 'interior-point', ...
        'Display', 'iter-detailed');

    % Define GlobalSearch solver
    gs = GlobalSearch('Display', 'final');

    % Define the optimization problem
    problem = createOptimProblem('fmincon', ...
        'x0', initial_guess, ...
        'objective', @(p) objective_rsquare(p, timeInfo, sir), ...
        'lb', lb, ...
        'ub', ub, ...
        'options', options);

    % Run GlobalSearch optimization
    paramsGKM = run(gs, problem);

    % Calculate fitted values on the given time grid
    y_fit = gkm_numeric(paramsGKM, timeInfo);

    % Calculate residuals and Goodness of Fit
    residuals = sir - y_fit;
    RMSE = sqrt(mean(residuals.^2));
    SSR = sum(residuals.^2);
    SST = sum((sir - mean(sir)).^2);
    R_squared = 1 - (SSR / SST);

    % Store goodness of fit metrics in a struct
    gof = struct('rmse', RMSE, 'rsquare', R_squared, 'ssr', SSR);

    % Generate fine grid and calculate fitted values
    x_fine = linspace(min(timeInfo), maxX, 1000);
    y_fit_fine = gkm_numeric(paramsGKM, x_fine);
end


function cost = objective_rsquare(params, timeInfo, sir)
    try
        y_fit = gkm_numeric(params, timeInfo);
        
        % If gkm_numeric outputs NaNs or wrong size, catch early
        if any(isnan(y_fit)) || length(y_fit) ~= length(sir)
            cost = 1e6;  % Big penalty
            return
        end

        residuals = sir - y_fit;
        SSR = sum(residuals.^2);
        SST = sum((sir - mean(sir)).^2);
        
        if SST == 0  % Prevent divide by zero
            cost = 1e6;
        else
            R_squared = 1 - (SSR / SST);
            cost = abs(1 - R_squared);
        end

    catch
        % If anything at all goes wrong, return large cost
        cost = 1e6;
    end
end