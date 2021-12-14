% Load the IRF data
%reconstructImage;

%filename = '~/experiments/2020/201103_IRF/diode_mag_1670Amp_10MHz_2998V_29kHz_3.5mmIris_1800s_m1.spc';
%filename = '~/experiments/2020/201103_IRF/NKT_mag_74perc_9-75MHz_2998V_33kHz_3.5mmIris_1800s_m1.spc';
%filename = '~/experiments/2020/201130_IRF_NKT/11.1MHz_71perc_NKT_f100mm_lens_1800s_55kHz_A_m1.spc';
%filename = '~/experiments/2020/201204_IRF/NKT_11-14MHz_69perc_f100mm_50kHz_lowMag_60x_1800s_A_m1.spc';
%filename = '/home/jakub/experiments/2020/201211_IRF/90xlowMag/11-1MHz_71.7perc_55kHz_1800s_A_m1.spc';
%filename = '/home/jakub/experiments/2020/201211_IRF/60xlowMag/11-1MHz_71.5perc_55kHz_1800s_A_m1.spc';
%filename = '/home/jakub/experiments/2020/201211_IRF/60xhighMag/11-1MHz_71.5perc_60kHz_1800s_A_m1.spc';
%filename = '/home/jakub/experiments/2020/201211_IRF/90xhighMag/11-1MHz_72.5perc_60kHz_1800s_A_m1.spc';
filename = '/home/jakub/experiments/2021/210903_IRF/MCP/loMag_FITC/11-14MHz_70-4perc_60x_loMag_f100mm_80kHz_FITC_filter_FITC_1800s_A_m1.spc';
input.frameTime = Inf;

% Model type: gauss, exGauss, shift
input.type = 'gauss';
input.Xstart = 1148;
input.Xend = 4008;
input.XYstep = 8;
input.Ystart = 656;
input.Yend = 3802;
% NKT 9.75MHz lowMag
input.Tstart0 = 1800;
input.Tend0 = 3799;
input.TstartIRF = 1980;
input.TendIRF = 2159;
% NKT 11.1MHz lowMag
input.Tstart0 = 2900;
input.Tend0 = 3799;
input.TstartIRF = 3000;
input.TendIRF = 3079;
% NKT 11.1MHz highMag
input.Tstart0 = 2900;
input.Tend0 = 3799;
input.TstartIRF = 3085;
input.TendIRF = 3164;

%input.Tstart = 1535;
%input.Tend = 1634;
input.Tstart = input.TstartIRF;
input.Tend = input.TendIRF;
input.Tstep = 1;
input.saveIRF = false;
input.IRFfile = [filename(1 : end - 7), '.', input.type, '.mat']; %'~/experiments/2020/201103_IRF/diode_mag_167Amp_10MHz_2998V_3.5mmIris.mat';
input.IRFicsFile = [input.IRFfile(1 : end - 3), 'ics'];
input.fitFile = [input.IRFfile(1 : end - 3), 'png'];
input.IRFrawPeakFile = [input.IRFfile(1 : end - 3), 'rawPeak.png'];
input.IRFcorrPeakFile = [input.IRFfile(1 : end - 3), 'corrPeak.png'];

% Read the IRF dataset
[XYZimage, param] = SPC2image(filename, '', input);

% Make a time-projection of the IRF dataset
image1 = sum(XYZimage, 3);

sensorIm = dip_image(image1);
% Find the imaging area
sensor = threshold(sensorIm, 'background');
%sensor = threshold(sensorIm, 'fixed', 200);
% Sensor needs removing areas around the edges, where the signal is a bit
% weak
sensor = erosion(sensor, 11);
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

%% Find the peak positions of the IRF
peakPos = IRFpeakFind(XYZimage, sensor, input);

%% Get rid of peak position outliers
peakPosNO = peakPos;
peakPosNO(isoutlier(peakPos, 'movmedian', 9, 'ThresholdFactor', 5)) = NaN;

%% Fit the IRFshift the with a 2D polynomials and get the residuals
% Create a figure
figure('WindowState', 'maximized')
% XY coordinates
XYcoord = [Xmat(~isnan(peakPosNO)), Ymat(~isnan(peakPosNO))];
% Fitted shifts
Zshifts = peakPos(~isnan(peakPosNO));
% Number of models
nrModels = 6;
% axes handles
ha = zeros(6, nrModels);
% histogram handles
hh = zeros(1, nrModels);
% IRF shift models
shiftModel = {'poly22', 'poly33', 'poly44', 'poly55', 'lowess'};
% IRF offset matrix
for i = 1 : (nrModels - 1)
    MCPshift.(shiftModel{i}).offset = zeros(size(peakPos));
    MCPshift.(shiftModel{i}).model = sfit;
    MCPshift.(shiftModel{i}).error.mean = 0;
    MCPshift.(shiftModel{i}).error.std = 0;
    if i == 1
        MCPshift.(shiftModel{i}).name = '2nd Degree 2D Polynomial';
    elseif i == 2
        MCPshift.(shiftModel{i}).name = '3rd Degree 2D Polynomial';
    elseif i == nrModels - 1
        MCPshift.(shiftModel{i}).name = 'Linear Smoothing regression';
    else
        MCPshift.(shiftModel{i}).name = ...
            sprintf('%dth Degree 2D Polynomial', i + 1);
    end
end
for i = 1 : nrModels
    if i == 1
        ha(1, 1) = subplot(3, nrModels, 1);
        surf(peakPos, 'EdgeColor', 'none')
        zl1 = get(ha(1, 1), 'ZLim');
        cl1 = get(ha(1, 1), 'CLim');
        title('IRF Shift')
        ha(3, 1) = subplot(3, nrModels, 2 * nrModels + 1);
        hh(1) = histogram(peakPosNO(~isnan(peakPosNO)));
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
        surf(peakPos, 'EdgeColor', 'none')
        % Give a title to the plot
        title(MCPshift.(shiftModel{i - 1}).name)
        % Create a plot of the residuals
        ha(2, i) = subplot(3, nrModels, i + nrModels);
        %plot(sf)
        % fit and data difference
        dShifts = sf(Xmat(~isnan(peakPosNO)), ...
                     Ymat(~isnan(peakPosNO))) - Zshifts;
        dImage = NaN(size(peakPos));
        dImage(~isnan(peakPosNO)) = dShifts;
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
    % Set the graph size right - not sure why it shrinks
    refPos = get(ha(3, 1), 'Position');
    curPos = get(ha(3, i), 'Position');
    curPos(1) = curPos(1) + curPos(3) / 2 - refPos(3) / 2;
    curPos(2 : 4) = refPos(2 : 4);
    set(ha(3, i), 'Position', curPos);
end

print(gcf, input.fitFile, '-dpng')

%% Create shift masks
hf = zeros(1, nrModels + 2);
hs = hf;
cl = [Inf, -Inf];
input.MCPfitFile = cell(1, nrModels - 1);
for i = 1 : (nrModels - 1)
    hf(i) = figure('Units', 'pixels', ...
                   'Position', [10, 10, 800, 700], ...
                   'PaperUnits', 'centimeter', ...
                   'PaperPosition', [0 0 24, 21], ...
                   'PaperSize', [24, 21]);
    mcpshift = MCPshift.(shiftModel{i}).offset;
    mcpshift(~sensor) = NaN;
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

%% Create a figure of the IRF peak position
hf(end - 1) = figure('Units', 'pixels', ...
                     'Position', [10, 10, 800, 700], ...
                     'PaperUnits', 'centimeter', ...
                     'PaperPosition', [0 0 24, 21], ...
                     'PaperSize', [24, 21]);
ha(6, 1) = axes;
% Calculate center value and median absolute deviation
peakPosMd = median(peakPos(:) * bin2ps, 'omitnan');
peakPosMad = mad(peakPos(:) * bin2ps);
% Create the plot
hs(end - 1) = surf(peakPos * bin2ps - peakPosMd, 'EdgeColor', 'none');
view(2)
ha(6, 2) = colorbar;
box on
axis equal
set(ha(6, 1), 'FontSize', 16)
title('IRF Peak Position (Raw)')
xlabel('X pixel');
ylabel('Y pixel');
set(ha(6, 1), 'XLim', [0, size(image1, 2)]);
set(ha(6, 1), 'YLim', [0, size(image1, 1)]);
set(ha(6, 1), 'CLim', 5 * peakPosMad * [-1, 1]);
set(get(ha(6, 2), 'Label'), 'String', 'IRF Timing Skew [ps]', ...
                            'FontSize', get(ha(6, 1), 'FontSize'));
print(hf(end - 1), input.IRFrawPeakFile, '-dpng')

% Save the IRF file
save(input.IRFfile, 'MCPshift', 'input')

%% Run the MCP shift correction on the same IRF data to obtain a corrected
%  IRF
% Update the Tstart and Tend
input.Tstart = input.Tstart0;
input.Tend = input.Tend0;
% Read the IRF dataset, this time correcting the MCP shift
XYZimage = SPC2image(filename, input.IRFfile, input);
% Find the peak positions of the IRF
peakPos = IRFpeakFind(XYZimage(:, :, (input.TstartIRF : input.TendIRF) - input.Tstart), ...
                      sensor, ...
                      input);
% Gate out the sensor part of the XYZ image stack
XYZimage(repmat(~sensor, [1 1 size(XYZimage, 3)])) = 0;
% Use a X and Y sum projections of the XYZ stack to extract the IRF
IRF = sum(sum(double(XYZimage), 1), 2);
% Scale IRF to UINT16
IRF = IRF / max(IRF) * double(intmax('uint16'));
% TAC bin width
binWidth = diff(param.tac([1, 2]));
% Save the IRF
saveIRF(IRF, [1 1], input.IRFicsFile, binWidth)

% Save the IRF file
save(input.IRFfile, '-append', 'IRF')



%% Create a figure of the IRF peak position
hf(end) = figure('Units', 'pixels', ...
                     'Position', [10, 10, 800, 700], ...
                     'PaperUnits', 'centimeter', ...
                     'PaperPosition', [0 0 24, 21], ...
                     'PaperSize', [24, 21]);
ha(6, 3) = axes;
% Calculate center value and median absolute deviation
peakPosMd = median(peakPos(:) * bin2ps, 'omitnan');
peakPosMad = mad(peakPos(:) * bin2ps);
% Create the plot
hs(end) = surf(peakPos * bin2ps - peakPosMd, 'EdgeColor', 'none');
view(2)
ha(6, 4) = colorbar;
box on
axis equal
set(ha(6, 3), 'FontSize', 16)
title('IRF Peak Position (Corrected)')
xlabel('X pixel');
ylabel('Y pixel');
set(ha(6, 3), 'XLim', [0, size(image1, 2)]);
set(ha(6, 3), 'YLim', [0, size(image1, 1)]);
set(ha(6, 3), 'CLim', get(ha(6, 1), 'CLim'))
set(get(ha(6, 4), 'Label'), 'String', 'IRF Timing Skew [ps]', ...
                            'FontSize', get(ha(6, 3), 'FontSize'));
print(hf(end), input.IRFcorrPeakFile, '-dpng')