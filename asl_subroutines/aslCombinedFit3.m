function [fit_result, gof, x_fine, y_fine, fracGaus, fracGamma, fit_bc] = aslCombinedFit3(x, y, param, varargin)
% Bi-component ASL fit: Gamma variate + Gaussian bolus
% Constraints:
%   dtG >= 0
%   dtGamma in [dtG, dtG + tau]   (via dtGamma = dtG + h*tau, 0<=h<=1)
%   perfC <= A                    (via perfC = A*k, 0<=k<=1)
%
% Extra rules:
%   1) If last measured point has y=0 AND max(x) < 2*T1:
%         treat as "no perfusion" and force Gamma component to zero.
%   2) Combined model must improve R^2 over Gaussian-only by at least minDeltaR2,
%         otherwise force Gamma component to zero.
%   3) Gaussian tail is scaled so it never goes negative (post-peak).
%   4) Gamma peak must be reasonably close to measured data (not far outside),
%         otherwise force Gamma component to zero.
%   5) Final fitted curve is clipped to be non-negative.

    % ---------- original data copy ----------
    T1     = param.T1;   % fixed T1 in the Gamma component
    x_orig = x(:).';
    y_orig = y(:).';

    if isempty(x_orig) || isempty(y_orig)
        x_orig = 0;
        y_orig = 0;
    end

    % max of the ORIGINAL data (kept for compatibility, even if unused)
    max_y_raw = max(y_orig); %#ok<NASGU>

    % ---------- no-perfusion check on ORIGINAL data ----------
    [xLast, idxLast] = max(x_orig);   % max TI in original data
    yLast            = y_orig(idxLast);
    noPerf           = (yLast == 0) && (xLast < 2*T1);

    % ---------- basic pre-processing (padding for stability) ----------
    % prepend (0,0)
    x = cat(2, 0, x_orig);
    y = cat(2, 0, y_orig);

    if numel(x) < 6
        x = zeros(1,6);
        y = zeros(1,6);
    end

    maxX = max(x);
    if nargin > 3 && isnumeric(varargin{1})
        maxX = varargin{1};
    end

    % ---------- component models ----------
    softstep = @(z,eps) 1 ./ (1 + exp(-z./eps));   % smoothed Heaviside

    epsG  = 50;  % ms, transition width for Gaussian
    epsGa = 50;  % ms, transition width for Gamma

    % Gamma variate with soft entrance at dtGamma
    gamma_variate = @(perfC, dtGamma, TI) ...
        perfC .* (TI - dtGamma) .* exp(-(TI - dtGamma) ./ T1) ...
          .* softstep(TI - dtGamma, epsGa);

    % ------------------------------------------------------------------
    % Re-parameterization for constraints in COMBINED model:
    %
    %  Fitted coefficients:  A, dtG, tau, k, h
    %
    %     perfC   = A * k           with 0 <= k <= 1  => perfC <= A
    %     dtGamma = dtG + h * tau   with 0 <= h <= 1  => dtG <= dtGamma <= dtG + tau
    % ------------------------------------------------------------------
    combinedModel = @(A, dtG, tau, k, h, TI) ...
        gaussian_bc(A, dtG, tau, TI, epsG, softstep) + ...
        gamma_variate(A .* k, dtG + h .* tau, TI);

    % ---------- initial guesses ----------
    maxy = max(y);
    A0   = maxy;
    dtG0 = max(0, mean(x)/4);
    tau0 = std(x);
    if ~isfinite(tau0) || tau0 <= 0
        tau0 = (max(x) - min(x)) / 6;
    end
    k0   = 0.5;
    h0   = 0.5;

    % COMBINED start/bounds: [A, dtG, tau, k, h]
    startC = [A0, dtG0, tau0, k0, h0];
    lowerC = [0,    0,    0,   0,  0];
    upperC = [2*maxy, maxX, 3000, 1, 1];

    % ---------- GAUSSIAN-ONLY model (for R^2 comparison) ----------
    ftG = fittype(@(A,dtG,tau,TI) ...
                  gaussian_bc(A,dtG,tau,TI,epsG,softstep), ...
                  'coefficients', {'A','dtG','tau'}, ...
                  'independent',  'TI');

    startG = [A0, dtG0, tau0];
    lowerG = lowerC(1:3);
    upperG = upperC(1:3);

    optsG = fitoptions('Method','NonlinearLeastSquares', ...
                       'StartPoint', startG, ...
                       'Lower',      lowerG, ...
                       'Upper',      upperG);

    [fitG, gofG] = fit(x', y', ftG, optsG);

    % ---------- COMBINED model fit ----------
    ftC = fittype(@(A,dtG,tau,k,h,TI) ...
                    combinedModel(A,dtG,tau,k,h,TI), ...
                  'coefficients', {'A','dtG','tau','k','h'}, ...
                  'independent',  'TI');

    optsC = fitoptions('Method','NonlinearLeastSquares', ...
                       'StartPoint', startC, ...
                       'Lower',      lowerC, ...
                       'Upper',      upperC);

    [fitC, gofC] = fit(x', y', ftC, optsC);

    % By default, return the COMBINED fit object and its gof
    fit_result = fitC;
    gof        = gofC;

    % ---------- evaluate on fine grid ----------
    x_fine = linspace(min(x), maxX, 1000);

    % Extract combined parameters
    A   = fitC.A;
    dtG = fitC.dtG;
    tau = fitC.tau;
    k   = fitC.k;
    h   = fitC.h;

    % Derived parameters for Gamma
    perfC   = A .* k;
    dtGamma = dtG + h .* tau;

    % Component curves (from combined model, with Gaussian tail scaled inside gaussian_bc)
    gauss_y = gaussian_bc(A, dtG, tau, x_fine, epsG, softstep);
    gamma_y = gamma_variate(perfC, dtGamma, x_fine);

    % ---------- model selection / perfusion suppression ----------
    % Gate 0: heuristic "no perfusion"
    useGamma = true;
    if noPerf
        useGamma = false;
    end

    % Gate 1: R^2 improvement threshold (ΔR² must exceed minDeltaR2)
    if useGamma
        deltaR2    = gofC.rsquare - gofG.rsquare;
        minDeltaR2 = 0.01;   % tune if needed
        if deltaR2 < minDeltaR2
            useGamma = false;
        end
    end

    % Gate 2: Gamma peak must be reasonably close to measured data
    if useGamma
        minXdata = min(x_orig);
        maxXdata = max(x_orig);
        [PHGa_pre, idxGa_pre] = max(gamma_y);
        TTPGamma_pre = x_fine(idxGa_pre);
        margin = 0.2 * (maxXdata - minXdata);  % allow slight extrapolation

        if PHGa_pre <= 0 || isnan(TTPGamma_pre) || TTPGamma_pre > maxXdata + margin
            useGamma = false;
        end
    end

    if ~useGamma
        % Kill Gamma component
        perfC   = 0;
        gamma_y = zeros(size(x_fine));
        dtGamma = NaN;  % marker
    end

    % ---------- combine and cleanup ----------
    y_fine  = gauss_y + gamma_y;
    y_fine(y_fine < 0) = 0;
    y_fine  = y_fine.';   % column vector

    % ---------- summary metrics ----------
    [PHG,  idxG]  = max(gauss_y);

    if any(gamma_y > 0)
        [PHGa, idxGa] = max(gamma_y);
        TTPGamma      = x_fine(idxGa);
    else
        PHGa     = 0;
        TTPGamma = NaN;
    end

    fit_bc.PHGauss     = PHG;
    fit_bc.TTPGauss    = x_fine(idxG);
    fit_bc.PHGamma     = PHGa;
    fit_bc.TTPGamma    = TTPGamma;
    fit_bc.deltaTGauss = dtG;
    fit_bc.deltaTGamma = dtGamma;
    fit_bc.tau         = tau;
    fit_bc.A           = A;
    fit_bc.perfC       = perfC;
    fit_bc.gauss       = gauss_y';
    fit_bc.gamma       = gamma_y';

    if PHGa > 0
        fit_bc.ratio = PHG / PHGa;
    else
        fit_bc.ratio = NaN;  % no gamma peak
    end

    % ---------- fractions ----------
    if perfC > 0
        fracGaus  = A / (A + perfC);
        fracGamma = perfC / (A + perfC);
    else
        fracGaus  = 1;
        fracGamma = 0;
    end

end

% ================= local helper =================
function yG = gaussian_bc(A, dtG, tau, TI, epsG, softstep)
% Gaussian bolus component with:
%   - soft entrance at dtG
%   - post-peak tail scaled so it never goes negative

    % Base "Eq.2 style" Gaussian with soft entrance
    base = A .* ( exp( -((TI - dtG - 2*tau).^2) ./ (2*tau.^2) ) - exp(-2) ) ...
             .* softstep(TI - dtG, epsG);

    % --- Tail scaling to avoid negatives after the peak ---
    [~, idxPeak] = max(base);
    tail    = base(idxPeak:end);
    minTail = min(tail);

    if minTail < 0
        peakVal = tail(1);   % value at the peak

        if peakVal > 0 && peakVal > minTail
            % Linear remap on the tail:
            %   tail_peak stays at peakVal
            %   minTail -> 0
            a = peakVal / (peakVal - minTail);
            b = -peakVal * minTail / (peakVal - minTail);
            tail = a * tail + b;
        else
            % Fallback: just clamp if something degenerate happens
            tail(tail < 0) = 0;
        end

        base(idxPeak:end) = tail;
    end

    % Final safety: no negatives at all
    base(base < 0) = 0;

    yG = base;
end