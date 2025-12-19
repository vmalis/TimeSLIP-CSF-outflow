function [fit_result, gof, x_fine, y_fine] = aslGammaFitCSF2(x, y, varargin)
%--------------------------------------------------------------------------
%   Smooth Gamma-variate fit for ASL inflow curves
%   Matches style of bi-component model (soft onset)
%
%   INPUT:
%       x       time array (ms)
%       y       ASL signal
%       maxX    (optional) fine-grid limit
%       -T1 X   (optional) fix T1 to value X
%
%   OUTPUT (same as before):
%       fit_result
%       gof
%       x_fine
%       y_fine
%--------------------------------------------------------------------------

    % ----------------------- Parse inputs -----------------------
    maxX = max(x);
    T1_fixed = false;
    T1_default = 4000;

    i = 1;
    while i <= length(varargin)
        if isnumeric(varargin{i})
            maxX = varargin{i};
        elseif ischar(varargin{i}) && strcmp(varargin{i}, '-T1')
            T1 = varargin{i+1};
            T1_fixed = true;
            i = i + 1;
        end
        i = i + 1;
    end

    if ~T1_fixed
        T1 = T1_default;
    end

    % ----------------------- Pre-processing -----------------------
    x = [0, x];
    y = [0, y];
    % x = cat(2, x, x(end)+T1);
    % y = cat(2, y, 0);

    max_y = max(y);

    % ------------------------ Soft step ---------------------------
    softstep = @(z,eps) 1 ./ (1 + exp(-z./eps));
    epsS = 50;   % smoothing width

    % ---------------------- Gamma variate -------------------------
    gamma_smooth = @(perfC, T1, dt, TI) ...
        perfC .* (TI - dt) .* exp(-(TI - dt)./T1) .* softstep(TI - dt, epsS);

    % -------------------------- Bounds ----------------------------
    % perfC ≤ max_y
    % dtGamma in [0,1000]
    % T1 in [1,8000]

    if T1_fixed
        % --------- FIXED T1 ----------
        fit_type = fittype(@(perfC, dt, TI) ...
            gamma_smooth(perfC, T1, dt, TI), ...
            'coefficients', {'perfC','dt'}, ...
            'independent','TI');

        start = [0.5*max_y, 300];

        lb = [0,    0];
        ub = [max_y, 1000];

        opts = fitoptions('Method','NonlinearLeastSquares',...
                          'StartPoint', start,...
                          'Lower', lb,...
                          'Upper', ub);

    else
        % --------- FREE T1 ----------
        fit_type = fittype(@(perfC, T1, dt, TI) ...
            gamma_smooth(perfC, T1, dt, TI), ...
            'coefficients', {'perfC','T1','dt'}, ...
            'independent','TI');

        start = [0.5*max_y,  T1_default/100, 300];

        lb = [0,   200,    0];
        ub = [max_y, 8000, 1000];

        opts = fitoptions('Method','NonlinearLeastSquares',...
                          'StartPoint', start,...
                          'Lower', lb,...
                          'Upper', ub);
    end

    % ------------------------- FIT -------------------------------
    [fit_result, gof] = fit(x', y', fit_type, opts);

    % ---------------------- Fine grid eval ------------------------
    x_fine = linspace(min(x), maxX, 1000);

    if T1_fixed
        y_fine = gamma_smooth(fit_result.perfC, T1, fit_result.dt, x_fine);
    else
        y_fine = gamma_smooth(fit_result.perfC, fit_result.T1, fit_result.dt, x_fine);
    end

    y_fine = max(y_fine, 0);   % no negative values
    y_fine=y_fine';
end