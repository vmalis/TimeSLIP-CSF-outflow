function [TTP, FWHM, peakHeight, AUC] = aslFitParams( x_fine, y_fine)
    % Calculate the peak height and the corresponding time (TTP)
    [peakHeight, peakIndex] = max(y_fine);
    TTP = x_fine(peakIndex);
    
    % Calculate Full Width at Half Maximum (FWHM)
    halfMax = peakHeight / 2;
    indicesAboveHalfMax = find(y_fine >= halfMax);
    if length(indicesAboveHalfMax) >= 2
        FWHM = x_fine(indicesAboveHalfMax(end)) - x_fine(indicesAboveHalfMax(1));
    else
        FWHM = '-'; % If FWHM can't be determined
    end
    
    % Calculate Area Under the Curve (AUC) using numerical integration
    AUC = trapz(x_fine, y_fine);
end