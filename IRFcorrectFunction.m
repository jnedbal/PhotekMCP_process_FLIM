% Load the IRF data
reconstructImage;

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

% Extract the IRFs from the sensor only
IRFtraces = dip_image(TPimage(:, sensorGate(:)));
% Calculate the average IRF across the sensor
avgTrace = squeeze(mean(IRFtraces, [], 1));

% Create matrix of shifts
IRFshift = NaN(size(image1));

% Run through the pixels of the sensor
sensorPix = find(sensorGate)';
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
    % Calculate the shift of the IRF in the pixel compared to the average
    % IRF
    IRFshift(Ymat(pixIn), Xmat(pixIn)) = ...
        findshift(avgTrace, squeeze(IRFtraces(i - 1, :)), 'iter');
end
delete(h)

%% Fit the IRFshift the with a 2D polynomials and get the residuals
% Create a figure
figure('Units', 'Normalized', 'Outerposition', [0 0 1 1])
% XY coordinates
XYcoord = [Xmat(~isnan(IRFshift)), Ymat(~isnan(IRFshift))];
% Fitted shifts
Zshifts = IRFshift(~isnan(IRFshift));
% axes handles
ha = zeros(3, 5);
% histogram handles
hh = zeros(1, 5);
for i = 1 : 5
    if i == 1
        ha(1, 1) = subplot(3, 5, 1);
        surf(IRFshift, 'EdgeColor', 'none')
        title('IRF Shift')
        ha(3, 1) = subplot(3, 5, 10 + 1);
        hh(1) = histogram(IRFshift(~isnan(IRFshift)));
    else
        % surface fit
        sf = fit(XYcoord, Zshifts, sprintf('poly%d%d', i, i));
        ha(1, i) = subplot(3, 5, i);
        plot(sf)%, [Xmat(~isnan(IRFshift)), Ymat(~isnan(IRFshift))], IRFshift(~isnan(IRFshift)))
        hold on
        surf(IRFshift, 'EdgeColor', 'none')
        if i < 3
            title('2nd Degree 2D Polynomial')
        else
            title(sprintf('%dth Degree 2D Polynomial', i))
        end
        ha(2, i) = subplot(3, 5, i + 5);
        %plot(sf)
        % fit and data difference
        dShifts = sf(Xmat(~isnan(IRFshift)), Ymat(~isnan(IRFshift))) - ...
                  Zshifts;
        dImage = NaN(size(IRFshift));
        dImage(~isnan(IRFshift)) = dShifts;
        surf(dImage, 'EdgeColor', 'none')
        title('Residual')
        % Plot the histogram of differences
        ha(3, i) = subplot(3, 5, 10 + i);
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
zl = cell2mat(get(ha(1, :), 'ZLim'));
zl = [min(zl(:, 1)), max(zl(:, 2))];
set(ha(1, :), 'ZLim', zl)
cl = cell2mat(get(ha(1, :), 'CLim'));
cl = [min(cl(:, 1)), max(cl(:, 2))];
set(ha(1, :), 'CLim', cl)
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
    axes(ha(3, i));
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
    axes(ha(3, i));
    title('Residual histogram')
    text(xl(2), yl(2), sprintf('Mean: %0.3f\nSTD: %0.3f', mean(get(hh(i), 'Data')), std(get(hh(i), 'Data'))), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
end