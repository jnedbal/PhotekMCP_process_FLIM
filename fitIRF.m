
%% Create a waitbar
hwaitbar = waitbar(0, '', ...
                   'Name', 'Fitting IRF with exGaussian function...');
% Make it a big higher
hwaitbar.Position(4) = hwaitbar.Position(4) + 10;



%% Find the position of the peak in each pixel
fit.h = zeros(inputSize([1, 2]));
fit.mu = fit.h;
fit.sigma = fit.h;
fit.tau = fit.h;
fit.offset = fit.h;
fit.exitFlag = fit.h;
peak.Pos = nan(inputSize([1,2]));
peak.Max = peak.Pos;
peak.risingHM = peak.Pos;
peak.fallingHM = peak.Pos;

tic
numberPixels = double(numberPixels); %#ok<NODEF>
                                     % loaded by extractHistogramData 

%% Go through the image pixel by pixel
for i = 1 : numberPixels
    % update the waitbar
    if mod(i, 100) == 0
        if ishandle(hwaitbar)
            % Estimated time to completion
            eta = numberPixels / i * toc / 86400;
            if isinf(eta); eta = 0; end
            % Update the waitbar
            waitbar(i / numberPixels, ...
                    hwaitbar, ...
                    sprintf('%d of %d\n%s (%s ETA)', ...
                            i, ...
                            numberPixels, ...
                            datestr(toc / 86400, 'HH:MM:SS'), ...
                            datestr(eta, 'HH:MM:SS')));
        end
    end
    % Extract the histogram for pixel
    pixIndex = double(1 : (lastCalBin(i) - firstCalBin(i) + 1))';
    %pixHist = double(inputHist(firstCalBin(i) : lastCalBin(i), i));
    pixHist = double(XYZimage(pixIndex, i));
    % Find the maximum in each pixel
    [h0, mu0] = max(pixHist);
    sigma0 = 3;
    tau0 = 8;
    offset0 = median(pixHist);
    % These are the initial guesses for the fit
    param0 = [h0 * tau0 * exp(1), mu0, sigma0, tau0, offset0];
    % Run the fit
    [fit.h(i), ...
        fit.mu(i), ...
        fit.sigma(i), ...
        fit.tau(i), ...
        fit.offset(i), ...
        fit.exitFlag(i)] = exgfit(param0);
    %fprintf('%d\n', i);
end

% Reshape outputs into matrices of the same size as the image sensor
%h = reshape(h, inputSize([1, 2]));
%mu = reshape(mu, inputSize([1, 2]));
%sigma = reshape(sigma, inputSize([1, 2]));
%tau = reshape(tau, inputSize([1, 2]));
%exitFlag = reshape(exitFlag, inputSize([1, 2]));

%% Calculate the well fitted bins
% These are ones that give results similar to the median
fit.goodfit = (fit.sigma < 1.5 * median(fit.sigma(:))) & ...
              (fit.sigma > 0.5 * median(fit.sigma(:))) & ...
              (~isoutlier(fit.offset)) & ...
              (fit.tau < 1.5 * median(fit.tau(:))) & ...
              (fit.tau > 0.5 * median(fit.tau(:)));
% The first and last row are alway poor on my SPAD, never use them
fit.goodfit([1, end], :) = false;

%% Calculate the fit parameters
tic
% Set options for the minimization function
options = optimoptions('fmincon', 'Display', 'none');
for i = find(fit.goodfit)'
    if mod(i, 100) == 0
        if ishandle(hwaitbar)
            % Estimated time to completion
            eta = numberPixels / i * toc / 86400;
            if isinf(eta); eta = 0; end
            % Update the waitbar
            waitbar(i / numberPixels, ...
                    hwaitbar, ...
                    sprintf('%d of %d\n%s (%s ETA)', ...
                            i, ...
                            numberPixels, ...
                            datestr(toc / 86400, 'HH:MM:SS'), ...
                            datestr(eta, 'HH:MM:SS')));
        end
    end
    % Find the position of the peak, and the the maximum of the peak
    % Minimization function
    mFun = @(x) (-(exGauss(x, ...
                           fit.h(i), ...
                           fit.mu(i), ...
                           fit.sigma(i), ...
                           fit.tau(i), ...
                           fit.offset(i))));
    % Run minimization routine to find the peak position
    peak.Pos(i) = fminsearch(mFun, fit.mu(i));
    % Evaluate the exGauss function to get the value of the peak
    peak.Max(i) = exGauss(peak.Pos(i), ...
                          fit.h(i), ...
                          fit.mu(i), ...
                          fit.sigma(i), ...
                          fit.tau(i), ...
                          fit.offset(i));
    % Create a minimization function for half-maximum
    mFun = @(x) abs(peak.Max(i) / 2 - ...
                    exGauss(x, fit.h(i), fit.mu(i), fit.sigma(i), ...
                               fit.tau(i), fit.offset(i)) + ...
                    fit.offset(i) / 2);
    % Run minimization routine to find the left half-maximum
    peak.risingHM(i) = fmincon(mFun, ...
                               peak.Pos(i) - fit.sigma(i), ...
                               [], [], [], [], ...
                               -Inf, ...
                               peak.Pos(i), ...
                               [], options);
    % Run minimization routine to find the right half-maximum
    peak.fallingHM(i) = fmincon(mFun, ...
                                peak.Pos(i) + fit.tau(i), ...
                                [], [], [], [], ...
                                peak.Pos(i), ...
                                Inf, ...
                                [], options);
end
% Calculate the full width at half-maximum
peak.FWHM = peak.fallingHM - peak.risingHM;

%% Run interpolations for the data that was not properly fitted
% interpolate the FWHM for points that are not fitted well
% Create a matrix of points
[X, Y] = meshgrid(1 : size(peak.FWHM, 2), 1 : size(peak.FWHM, 1));
Xv = X(~fit.goodfit);
Yv = Y(~fit.goodfit);
X = X(fit.goodfit);
Y = Y(fit.goodfit);

% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the full-width half maximum for the missing pixels
F = scatteredInterpolant(X, Y, peak.FWHM(fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peak.FWHMinterp = peak.FWHM;
peak.FWHMinterp(~fit.goodfit) = F(Xv, Yv);



% Interpolate the peak position for points that are not fitted well
% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the peak position for the missing pixels
F = scatteredInterpolant(X, Y, peak.Pos(fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peak.PosInterp = peak.Pos;
peak.PosInterp(~fit.goodfit) = F(Xv, Yv);



% Interpolate the rising edge half-maximum for pixels not fitted well
% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the peak position for the missing pixels
F = scatteredInterpolant(X, Y, peak.risingHM(fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peak.risingHMinterp = peak.risingHM;
peak.risingHMinterp(~fit.goodfit) = F(Xv, Yv);



% Interpolate the falling edge half-maximum for pixels not fitted well
% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the peak position for the missing pixels
F = scatteredInterpolant(X, Y, peak.fallingHM(fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peak.fallingHMinterp = peak.fallingHM;
peak.fallingHMinterp(~fit.goodfit) = F(Xv, Yv);

% Close the waitbar
delete(hwaitbar)

% Assign results to the main structure for future use
correction.IRF.fit = fit;
correction.IRF.peak = peak;




%% Run a routine to correct the IRF offset
correction.IRF.corrected = resampleHistogramPar(correction.IRF.raw);