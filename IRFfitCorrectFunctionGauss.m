% Load the IRF data
%reconstructImage;
global pixHist
global pixIndex

%filename = '~/experiments/2020/201103_IRF/diode_mag_1670Amp_10MHz_2998V_29kHz_3.5mmIris_1800s_m1.spc';
filename = '~/experiments/2020/201103_IRF/NKT_mag_74perc_9-75MHz_2998V_33kHz_3.5mmIris_1800s_m1.spc';
input.frameTime = Inf;
input.Xstart = 1168;
input.Xend = 4008;
input.XYstep = 8;
input.Ystart = 846;
input.Yend = 3802;
input.Tstart0 = 1500;
input.Tend0 = 3500;
input.Tstart = 2020;
input.Tend = 2119;
input.Tstep = 1;
input.saveIRF = false;
%input.IRFfile = '~/experiments/2020/201103_IRF/diode_mag_167Amp_10MHz_2998V_3.5mmIris.gauss.mat';
input.IRFfile = '~/experiments/2020/201103_IRF/NKT_mag_74perc_9-75MHz_2998V_33kHz_3.5mmIris.gauss.mat';

input.IRFicsFile = [input.IRFfile(1 : end - 3), 'ics'];
input.fitFile = [input.IRFfile(1 : end - 3), 'png'];

% Read the IRF dataset
[XYZimage, param] = SPC2image(filename, '', input);

% Make a time-projection of the IRF dataset
image1 = sum(XYZimage, 3);

sensorIm = dip_image(image1);
% Find the imaging area
sensor = threshold(sensorIm, 'background');
% Sensor needs removing areas around the edges, where the signal is a bit
% weak
sensor = erosion(sensor);
%sensor = dip_growregions(erosion(sensor), [], [], 2, 3, 'low_first');
% Label all regions
%sensorLab = label(sensor);
% Find the center of the sensor
%sensorMeas = measure(sensorLab, sensorIm, {'center', 'Mass', 'radius'});
% Find the largest region. Ignore the rest.
%sensorMeas = sensorMeas(sensorMeas.Mass == max(sensorMeas.Mass));

% Matrix of coordinates
X = xx(size(sensor), 'corner');
Y = yy(size(sensor), 'corner');
% Matrix of coordinates in Matlab variables
Xmat = double(X) + 1;
Ymat = flipud(double(Y)) + 1;

% Distance from center
%Cdist = sqrt((X - sensorMeas.Center(1)) .^ 2 + (Y - sensorMeas.Center(2)) .^ 2);
% Sensor gate
%sensorGate = Cdist < sensorMeas.Radius(3);
sensorGate = flipud(logical(sensor));

%% Run the lifetime fit for pixels with more that 10 photons
%  reshape XYZ image into an array where the first dimension is the time
%  and the other dimensions are all the pixels
TPimage = reshape(XYZimage, ...
                  size(XYZimage, 1) * size(XYZimage, 2), ...
                  size(XYZimage, 3))';
clear('XYZimage')
% Extract the IRFs from the sensor only
IRFtraces = double(TPimage(:, sensorGate(:)));
clear('TPimage')
% Calculate the average IRF across the sensor
avgTrace = mean(IRFtraces, 2);

% Create matrix of shifts
% IRFshift = NaN(size(image1));

% Run through the pixels of the sensor
sensorPix = find(sensorGate)';

%% Find the position of the peak in each pixel based on the fit of the
%  average IRF trace
% Find the maximum in each pixel
pixHist = avgTrace;
[h0, mu0] = max(pixHist);
sigma0 = 6;
% These are the initial guesses for the fit
param0 = [h0, mu0, sigma0];

% Create an index of the bins
pixIndex = (1 : size(IRFtraces, 1))';

% Run the fit
[h0, mu0, sigma0] = gfit(param0);
sigmaAvg = ones(100, 1) * sigma0;
sigma01 = sigma0;
Fit.h = NaN(size(sensorGate));  Fit.h(sensorGate) = h0;
Fit.mu = Fit.h;                 Fit.mu(sensorGate) = mu0;
Fit.sigma = Fit.h;              Fit.sigma(sensorGate) = sigma0;
param0 = [h0, mu0, sigma0];
param01 = param0;

%% Go through all the pixels and find the shift between the average IRF and
%  the pixel IRF
% Create a waitbar that's larger than normal
h = waitbar(0, 'Correcting IRF shifts...');
h.Position(4) = h.Position(4) * 1.3;
% Record the time of the analysis start
tic;
for i = 1 : numel(sensorPix)
    % Update the waitbar every 100 pixels
    if ~mod(i, 100)
        if ishandle(h)
            waitbar(i / numel(sensorPix), ...
                    h, ...
                    {'Correcting IRF shifts...', ...
                     sprintf('Elapsed time %s', duration(0, 0, toc)), ...
                     sprintf('Remaining time %s', ...
                             duration(0, 0, ...
                                      toc * (numel(sensorPix) - i) / i))})
        end
    end
    % Get the current pixel index
    pixIn = sensorPix(i);
    % Get the pixel histogram
    pixHist = IRFtraces(:, i);
    % Find the maximum in each pixel
    [param01(1), param01(2)] = max(pixHist);
    % These are the initial guesses for the fit
    param01(3) = sigma01;
    % Run the fit
    [fh, fmu, fsigma, ~, fval01] = gfit(param01);
    [h02, mu02, sigma02, ~, fval02] = gfit(param0);
    % Check which fit is better
    if fval01 < fval02
        Fit.h(pixIn) = fh;
        Fit.mu(pixIn) = fmu;
        Fit.sigma(pixIn) = fsigma;
    else
        Fit.h(pixIn) = h02;
        Fit.mu(pixIn) = mu02;
        Fit.sigma(pixIn) = sigma02;
    end
    in = mod(i - 1, 100) + 1;
    sigmaAvg(in) = Fit.sigma(pixIn);
    sigma01 = mean(sigmaAvg);
    
    % Calculate the shift of the IRF in the pixel compared to the average
    % IRF
    %IRFshift(Ymat(pixIn), Xmat(pixIn)) = ...
    %    findshift(avgTrace, squeeze(IRFtraces(i - 1, :)), 'iter');
end
delete(h)

%% Calculate the shift differential
% Position of the average peak
%pos0 = dip_subpixelmaxima(dip_image(avgTrace), ...
%                          dip_image(avgTrace > 0.1 * max(avgTrace)), ...
%                          'gaussian') + 1;
%IRFshift = Fit.mu - pos0;
IRFshift = Fit.mu - mean(Fit.mu(:), 'omitnan');

%off = mean([max(IRFshiftCorr(:)), min(IRFshiftCorr(:))]);figure('Units', 'Pixels', 'Position', [10, 10, 800, 700], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 24, 21], 'PaperSize', [24, 21]); surf((IRFshiftCorr - off) * diff(param.tac([1, 2])) * 1000, 'EdgeColor', 'none'); view(2); box on; title('IRF Peak Shift Corrected'); set(gca, 'ZLim', [-70, 70], 'CLim', [-70, 70], 'FontSize', 16);axis square;cb = colorbar;set(cb.Label, 'String', 'IRF Peak Shift [ps]', 'FontSize', 16); cb.FontSize = 16;savefig('~/experiments/2020/201103_IRF/CORR_IRF_Gauss.fig');print('~/experiments/2020/201103_IRF/CORR_IRF_Gauss.png', '-dpng')
%offN = mean([max(IRFshiftNCorr(:)), min(IRFshiftNCorr(:))]);figure('Units', 'Pixels', 'Position', [10, 10, 800, 700], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 24, 21], 'PaperSize', [24,21]); surf((IRFshiftNCorr - off) * diff(param.tac([1, 2])) * 1000, 'EdgeColor', 'none'); view(2); box on; title('IRF Peak Shift Raw'); set(gca, 'ZLim', [-70, 70], 'CLim', [-70, 70], 'FontSize', 16);axis square;cb = colorbar;set(cb.Label, 'String', 'IRF Peak Shift [ps]', 'FontSize', 16); cb.FontSize = 16;savefig('~/experiments/2020/201103_IRF/RAW_IRF_Gauss.fig');print('~/experiments/2020/201103_IRF/RAW_IRF_Gauss.png', '-dpng')

%% Fit the IRFshift the with a 2D polynomials and get the residuals
% Create a figure
figure('Units', 'Normalized', 'Outerposition', [0 0 1 1])
% XY coordinates
XYcoord = [Xmat(~isnan(IRFshift)), Ymat(~isnan(IRFshift))];
% Fitted shifts
Zshifts = IRFshift(~isnan(IRFshift));
% Number of models
nrModels = 6;
% axes handles
ha = zeros(3, nrModels);
% histogram handles
hh = zeros(1, nrModels);
% IRF shift models
shiftModel = {'poly22', 'poly33', 'poly44', 'poly55', 'lowess'};
% IRF offset matrix
for i = 1 : (nrModels - 1)
    MCPshift.(shiftModel{i}).offset = zeros(size(IRFshift));
    MCPshift.(shiftModel{i}).model = sfit;
    MCPshift.(shiftModel{i}).error.mean = 0;
    MCPshift.(shiftModel{i}).error.std = 0;
    if i == 1
        MCPshift.(shiftModel{i}).name = '2nd Degree 2D Polynomial';
    elseif i == 2
        MCPshift.(shiftModel{i}).name = '3rd Degree 2D Polynomial';
    elseif i == nrModels - 1
        MCPshift.(shiftModel{i}).name = 'Linear smooting regression';
    else
        MCPshift.(shiftModel{i}).name = ...
            sprintf('%dth Degree 2D Polynomial', i + 1);
    end
end
for i = 1 : nrModels
    if i == 1
        ha(1, 1) = subplot(3, nrModels, 1);
        surf(IRFshift, 'EdgeColor', 'none')
        zl1 = get(ha(1, 1), 'ZLim');
        cl1 = get(ha(1, 1), 'CLim');
        title('IRF Shift')
        ha(3, 1) = subplot(3, nrModels, 2 * nrModels + 1);
        hh(1) = histogram(IRFshift(~isnan(IRFshift)));
    else
        % surface fit
        sf = fit(XYcoord, Zshifts, shiftModel{i - 1}, 'Normalize', 'on');
        % Store the surface fit
        MCPshift.(shiftModel{i - 1}).model = sf;
        % Evaluate the surface fit at all pixel positions
        MCPshift.(shiftModel{i - 1}).offset = sf(Xmat, Ymat);
        ha(1, i) = subplot(3, nrModels, i);
        plot(sf)%, [Xmat(~isnan(IRFshift)), Ymat(~isnan(IRFshift))], IRFshift(~isnan(IRFshift)))
        hold on
        % Plot the surface of the calculated IRF shift
        surf(IRFshift, 'EdgeColor', 'none')
        % Give a title to the plot
        title(MCPshift.(shiftModel{i - 1}).name)
        % Create a plot of the residuals
        ha(2, i) = subplot(3, nrModels, i + nrModels);
        %plot(sf)
        % fit and data difference
        dShifts = sf(Xmat(~isnan(IRFshift)), Ymat(~isnan(IRFshift))) - ...
                  Zshifts;
        dImage = NaN(size(IRFshift));
        dImage(~isnan(IRFshift)) = dShifts;
        surf(dImage, 'EdgeColor', 'none')
        title('Residual')
        % Calculate the errors of the fit
        MCPshift.(shiftModel{i - 1}).error.mean = mean(dShifts);
        MCPshift.(shiftModel{i - 1}).error.std = std(dShifts);
        % Plot the histogram of differences
        ha(3, i) = subplot(3, nrModels, 2 * nrModels + i);
        hh(i) = histogram(dShifts(:));
    end
end

% Adapt the scales
% Top row
xl = cell2mat(get(ha(1, :), 'XLim'));
xl = [max(xl(:, 1)), min(xl(:, 2))];
set(ha(1, :), 'XLim', xl)
yl = cell2mat(get(ha(1, :), 'YLim'));
yl = [max(yl(:, 1)), min(yl(:, 2))];
set(ha(1, :), 'YLim', yl)
set(ha(1, :), 'ZLim', zl1)
set(ha(1, :), 'CLim', cl1)
% Middle row
xl = cell2mat(get(ha(2, 2 : end), 'XLim'));
xl = [max(xl(:, 1)), min(xl(:, 2))];
set(ha(2, 2 : end), 'XLim', xl)
yl = cell2mat(get(ha(2, 2 : end), 'YLim'));
yl = [max(yl(:, 1)), min(yl(:, 2))];
set(ha(2, 2 : end), 'YLim', yl)
zl = cell2mat(get(ha(2, 2 : end), 'ZLim'));
zl = [min(zl(:, 1)), max(zl(:, 2))];
set(ha(2, 2 : end), 'ZLim', zl)
cl = cell2mat(get(ha(2, 2 : end), 'CLim'));
cl = [min(cl(:, 1)), max(cl(:, 2))];
set(ha(2, 2 : end), 'CLim', cl)
% Bottom row
% Resample the histogram
% Bin limits
bl = [min(cellfun(@min, get(hh, 'BinEdges'))), ...
      max(cellfun(@max, get(hh, 'BinEdges')))];
% Bin width
bw =  max(cell2mat(get(hh, 'BinWidth')));
% Correct the bin limits to the bin width
bl = bw * [floor(bl(1) / bw), ceil(bl(2) / bw)];
% Bin edges
be = bl(1) : bw : bl(2);
for i = 1 : numel(hh)
    axes(ha(3, i)); %#ok<LAXES>
    hh(i) = histogram(get(hh(i), 'Data'), be);
end
xl = cell2mat(get(ha(3, :), 'XLim'));
xl = [min(xl(:, 1)), max(xl(:, 2))];
set(ha(3, :), 'XLim', xl)
yl = cell2mat(get(ha(3, :), 'YLim'));
yl = [min(yl(:, 1)), max(yl(:, 2))];
set(ha(3, :), 'YLim', yl)
% Add text
for i = 1 : numel(hh)
    axes(ha(3, i)); %#ok<LAXES>
    title('Residual histogram')
    text(xl(2), yl(2), sprintf('Mean: %0.3f\nSTD: %0.3f', mean(get(hh(i), 'Data')), std(get(hh(i), 'Data'))), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
end

print(gcf, input.fitFile, '-dpng')
save(input.IRFfile, 'MCPshift', 'input', 'Fit')

%% Create shift masks
hf = zeros(1, nrModels);
hs = hf;
cl = [Inf, -Inf];
input.MCPfitFile = cell(1, nrModels - 1);
for i = 1 : (nrModels - 1)
    hf(i) = figure;
    mcpshift = MCPshift.(shiftModel{i}).offset;
    mcpshift(~sensorGate) = NaN;
    ha(4, i) = axes('FontSize', 16);
    hs(i) = surf(mcpshift, 'EdgeColor', 'none');
    view(2);
    ha(5, i) = colorbar;
    box on
    axis equal
    title(MCPshift.(shiftModel{i}).name)
    xlabel('X pixel');
    ylabel('Y pixel');
    set(ha(4, i), 'XLim', [0, size(image1, 2)]);
    set(ha(4, i), 'YLim', [0, size(image1, 1)]);
    set(get(ha(5, i), 'Label'), 'String', 'IRF Timing Skew [bins]', ...
                                'FontSize', get(ha(4, i), 'FontSize'));
    cl(1) = min([cl(1), get(ha(5, i), 'Limits')]);
    cl(2) = max([cl(2), get(ha(5, i), 'Limits')]);
    input.MCPfitFile{1, i} = [input.fitFile(1 : end - 3), ...
                              shiftModel{i}, '.bin.png'];
    input.MCPfitFile{2, i} = [input.fitFile(1 : end - 3), ...
                              shiftModel{i}, '.time.png'];
end
set(ha(4, 1 : (nrModels - 1)), 'ZLim', cl, 'CLim', cl)
bin2ps = diff(param.tac([1, 2])) * 1000;
cl = cl * bin2ps;
for i = 1 : (nrModels - 1)
    print(hf(i), input.MCPfitFile{1, i}, '-dpng')
    set(get(ha(5, i), 'Label'), 'String', 'IRF Timing Skew [ps]');
    set(hs(i), 'ZData', get(hs(i), 'Zdata') * bin2ps)
    set(ha(4, i), 'ZLim', cl, 'CLim', cl)
    set(ha(5, i), 'Limits', cl)
    print(hf(i), input.MCPfitFile{2, i}, '-dpng')
end

%% Run the MCP shift correction on the same IRF data to obtain a corrected
%  IRF
% Update the Tstart and Tend
input.Tstart = input.Tstart0;
input.Tend = input.Tend0;
% Read the IRF dataset, this time correcting the MCP shift
XYZimage = SPC2image(filename, input.IRFfile, input);
% Gate out the sensor part of the XYZ image stack
XYZimage(~repmat(sensorGate, [1 1 size(XYZimage, 3)])) = 0;
% Use a X and Y sum projections of the XYZ stack to extract the IRF
IRF = sum(sum(double(XYZimage), 1), 2);
% Scale IRF to UINT16
IRF = IRF / max(IRF) * double(intmax('uint16'));
% TAC bin width
binWidth = 1e3 * diff(param.tac([1, 2]));
% Save the IRF
saveIRF(IRF, [1 1], input.IRFicsFile, binWidth)