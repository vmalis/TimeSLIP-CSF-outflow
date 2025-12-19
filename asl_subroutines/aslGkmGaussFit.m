function y_combined = aslGkmGaussFit(params, timeInfo)
%--------------------------------------------------------------------------
%
%   Function to fit ASL data to a combined numeric model and Gaussian model
%   7 parameters are passed and ony 3 are fixed
%
%
%   INPUT:  params - Parameter vector with fields:
%           frac_numeric - Fractional weight for the numeric model
%           T1 - T1 value for the numeric model
%           Delta_t - Time delay parameter for the numeric model
%           tau - Time constant parameter shared between models
%           f - Flow-related parameter for the numeric model
%           a - Amplitude parameter for the Gaussian model
%           b - Center (mean) parameter for the Gaussian model
%           timeInfo - Time values (1D array)       'timeInfo'    (double)
%
%   OUTPUT: y_combined - Combined model output values (1D array)    'y_combined'  (double)
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------


    % Unpack the parameters
    frac_numeric = params(1);
    frac_gaussian = 2-frac_numeric;
    
    % Parameters for the numeric function
    T1 = params(2);             %comes as fixed
    Delta_t = params(3);
    tau = params(4);
    f = params(5);
    
    % Parameters for the Gaussian function
    a = params(6);              %comes fixed
    b = params(7);              %comes fixed
    c = tau;


    % Compute the numeric function component
    y_numeric = gkm_numeric([T1, Delta_t, tau, f], timeInfo);
    
    % Compute the Gaussian function component
    gaussian = @(a, b, c, t) a * (exp(-((t - b).^2 / (2 * c^2))) - exp(-((0 - b).^2 / (2 * c^2))));
    y_gaussian = gaussian(a, b, c, timeInfo);
    
    % Combine the two functions
    y_combined = frac_numeric * y_numeric + frac_gaussian * y_gaussian;

end