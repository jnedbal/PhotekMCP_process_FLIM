function peakPos = IRFpeakFind(XYZimage, sensor, input)
global pixHist
global pixIndex

% Distance from center
%Cdist = sqrt((X - sensorMeas.Center(1)) .^ 2 + (Y - sensorMeas.Center(2)) .^ 2);
% Sensor gate
%sensorGate = Cdist < sensorMeas.Radius(3);
sensorGate = flipud(logical(sensor));
% Run through the pixels of the sensor
sensorPix = find(sensorGate)';
% Number of sensor pixels
totalPix = numel(sensorPix);

%% Run the lifetime fit for pixels with more that 10 photons
%  reshape XYZ image into an array where the first dimension is the time
%  and the other dimensions are all the pixels
TPimage = reshape(XYZimage, ...
                  size(XYZimage, 1) * size(XYZimage, 2), ...
                  size(XYZimage, 3))';
clear('XYZimage')
switch input.type
    case 'shift'
        % Extract the IRFs from the sensor only
        IRFtraces = dip_image(TPimage(:, sensorGate(:)), 'double');
        % Calculate the average IRF across the sensor
        avgTrace = squeeze(mean(IRFtraces, [], 1));
        % Create matrix of shifts
        peakPos = NaN(size(sensorGate));
    otherwise
        % Extract the IRFs from the sensor only
        IRFtraces = double(TPimage(:, sensorGate(:)));
        % Calculate the average IRF across the sensor
        avgTrace = mean(IRFtraces, 2);

        % Find the maximum in each pixel
        pixHist = avgTrace;
        [h0, mu0] = max(pixHist);

        % Create an index of the bins
        pixIndex = (1 : size(IRFtraces, 1))';
end
clear('TPimage')

switch input.type
    case 'gauss'
        %% Find the position of the peak in each pixel based on the fit of
        %  the IRF trace in each pixel and localilzing the maximum
        %  Model sigma
        sigma0 = 6;
        % These are the initial guesses for the fit
        param0 = [h0, mu0, sigma0];

        % Run the fit of the average IRF
        [h0, mu0, sigma0] = gfit(param0);
        % Create an average sigma vector that is going to evolve as more
        % IRFs are going to be fitted. To start, it just takes the sigma of
        % the average IRF trace
        sigmaAvg = ones(100, 1) * sigma0;
        % Create another variable that is going to change with each fit
        sigma01 = sigma0;
        % Create a struct with parameters of the fit
        Fit.h = NaN(size(sensorGate));  Fit.h(sensorGate) = h0;
        peakPos = Fit.h;            peakPos(sensorGate) = mu0;
        Fit.sigma = Fit.h;              Fit.sigma(sensorGate) = sigma0;
        % Create starting fit parameter vectors
        param0 = [h0, mu0, sigma0];
        param01 = param0;
    case 'exGauss'
        %% Find the position of the peak in each pixel based on the fit of
        %  the IRF trace in each pixel and localilzing the maximum
        %  Model sigma
        sigma0 = 3;
        tau0 = 8;
        offset0 = median(pixHist);
        % These are the initial guesses for the fit
        param0 = [h0 * tau0 * exp(1), mu0, sigma0, tau0, offset0];

        % Run the fit of the average IRF
        [h0, mu0, sigma0, ~, offset0] = exgfit(param0);
        tau0 = 1;     % The decay is better fit by a shorter tau
        % Create average sigma and tau vectors that are going to evolve as
        % more IRFs are going to be fitted. To start, it just takes the
        % sigma and tau of the average IRF trace
        sigmaAvg = ones(100, 1) * sigma0;
        tauAvg = ones(100, 1) * tau0;
        % Create extra variables that are going to change with each fit
        sigma01 = sigma0;
        tau01 = tau0;
        % Create a struct with parameters of the fit
        Fit.h = NaN(size(sensorGate));  Fit.h(sensorGate) = h0;
        Fit.mu = Fit.h;                 Fit.mu(sensorGate) = mu0;
        Fit.sigma = Fit.h;              Fit.sigma(sensorGate) = sigma0;
        Fit.tau = Fit.h;                Fit.tau(sensorGate) = tau0;
        Fit.offset = Fit.h;             Fit.offset(sensorGate) = offset0;
        peakPos = Fit.h;
        Fit.exitFlag = logical(size(Fit.h));
        % Create starting fit parameter vectors
        param0 = [h0 * tau0 * exp(1), mu0, sigma0, tau0, offset0];
        param01 = param0;
end

%% Go through all the pixels and find the shift between the average IRF and
%  the pixel IRF
% Create a waitbar that's larger than normal
h = waitbar(0, 'Correcting IRF shifts...');
h.Position(4) = h.Position(4) * 1.3;
% Record the time of the analysis start
tic;
switch input.type
    case 'shift'
        for i = 1 : totalPix
            % Update the waitbar every 100 pixels
            updateWaitbar(h, i, totalPix);

            % Get the current pixel index
            pixIn = sensorPix(i);
            % Calculate the shift of the IRF in the pixel compared to the
            % average IRF
            peakPos(pixIn) = findshift(avgTrace, ...
                                           squeeze(IRFtraces(i - 1, :)),...
                                           'iter');
        end
    case 'exGauss'
        for i = 1 : numel(sensorPix)
            % Update the waitbar every 100 pixels
            updateWaitbar(h, i, totalPix);

            % Get the current pixel index
            pixIn = sensorPix(i);
            % Get the pixel histogram
            pixHist = IRFtraces(:, i);
            % Find the maximum in each pixel
            [h01, param01(2)] = max(pixHist);
            % These are the initial guesses for the fit
            param01(1) = h01 * tau01 * exp(1);
            param01(3) = sigma01;
            param01(4) = tau01;
            % Run the fit based on the latest parameters
            [fh, fmu, fsigma, ftau, foffset, ~, fval01] = exgfit(param01);
            % Run the fit based on the average IRF parameter
            [h02, mu02, sigma02, tau02,offset02,~,fval02] = exgfit(param0);
            % Check which fit is better
            if fval01 < fval02
                % If the fit with the latest parameters is better, keep it
                Fit.h(pixIn) = fh;
                Fit.mu(pixIn) = fmu;
                Fit.sigma(pixIn) = fsigma;
                Fit.tau(pixIn) = ftau;
                Fit.offset(pixIn) = foffset;
            else
                % If the fit with the average IRF parameters is better, 
                % keep it
                Fit.h(pixIn) = h02;
                Fit.mu(pixIn) = mu02;
                Fit.sigma(pixIn) = sigma02;
                Fit.tau(pixIn) = tau02;
                Fit.offset(pixIn) = offset02;
            end
            % Update the sliding window of the last 100 fit parameters
            in = mod(i - 1, 100) + 1;
            sigmaAvg(in) = Fit.sigma(pixIn);
            sigma01 = mean(sigmaAvg);
            tauAvg(in) = Fit.tau(pixIn);
            tau01 = mean(tauAvg);
            % Model the curve and fin its peak
            IRFmodel = exGauss(pixIndex, ...
                               Fit.h(pixIn), ...
                               Fit.mu(pixIn), ...
                               Fit.sigma(pixIn), ...
                               Fit.tau(pixIn), ...
                               Fit.offset(pixIn));
            % Find the peak and its value
            [p, v] = dip_subpixelmaxima(...
                            dip_image(IRFmodel), ...
                            dip_image(IRFmodel > 0.1 * max(IRFmodel)), ...
                            'gaussian');
            % Store the peak position
            peakPos(pixIn) = p(v == max(v)) + 1;
        end
    case 'gauss'
        for i = 1 : numel(sensorPix)
            % Update the waitbar every 100 pixels
            updateWaitbar(h, i, totalPix);

            % Get the current pixel index
            pixIn = sensorPix(i);
            % Get the pixel histogram
            pixHist = IRFtraces(:, i);
            % Find the maximum in each pixel
            [param01(1), param01(2)] = max(pixHist);
            % These are the initial guesses for the fit
            param01(3) = sigma01;
            % Run the fit based on the latest parameters
            [fh, fmu, fsigma, ~, fval01] = gfit(param01);
            % Run the fit based on the average IRF parameter
            [h02, mu02, sigma02, ~, fval02] = gfit(param0);
            % Check which fit is better
            if fval01 < fval02
                % If the fit with the latest parameters is better, keep it
                Fit.h(pixIn) = fh;
                peakPos(pixIn) = fmu;
                Fit.sigma(pixIn) = fsigma;
            else
                % If the fit with the average IRF parameters is better, 
                % keep it
                Fit.h(pixIn) = h02;
                peakPos(pixIn) = mu02;
                Fit.sigma(pixIn) = sigma02;
            end
            % Update the sliding window of the last 100 fit parameters
            in = mod(i - 1, 100) + 1;
            sigmaAvg(in) = Fit.sigma(pixIn);
            sigma01 = mean(sigmaAvg);
        end
end
delete(h)

%% Create a zero offset to the peak position
peakPos = peakPos - mean(peakPos(:), 'omitnan');
end


function updateWaitbar(h, i, total)
    % Update the waitbar every 100 pixels
    if ~mod(i, 100)
        if ishandle(h)
            waitbar(i / total, ...
                    h, ...
                    {'Correcting IRF shifts...', ...
                     sprintf('Elapsed time %s', duration(0, 0, toc)), ...
                     sprintf('Remaining time %s', ...
                             duration(0, 0, ...
                                      toc * (total - i) / i))})
        end
    end
end