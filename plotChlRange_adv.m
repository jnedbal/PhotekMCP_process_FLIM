close all


% Load the conversion from percentage to PPFD
load('~/experiments/2021/210122_powerCalibration/powerDensity.mat')
% Load the pixel size
load('~/experiments/2021/210204_AFpattern/pixelSize.mat')

%% Load FIFO data
% Create a wildcard for looking for the right type of files% Create a wildcard for looking for the right type of files
wildcard = sprintf('%s/experiments/%s/', ...
                   char(java.lang.System.getProperty('user.home')), ...
                   datetime('now', 'Format', 'yyyy'));
% Create a dialog box to load the file
%[filenamelist, pathname] = ...
%    uigetfile({'*m1.spc', 'SPC-files (*m1.spc)'; ...
%               '*.*',     'All Files (*.*)'}, ...
%               'Pick a file', ...
%               'MultiSelect', 'on', ...
%               wildcard);
[fnames, path] = ...
    uigetfile('*Tau.ics', 'Select the INPUT lifetime files', ...
              'MultiSelect', 'on', wildcard);
if isequal(fnames, 0)
    return
end

% Set the data to plot
dataTypes = {'Tau', 'Phasor', 'Photon Count', 'Amplitude'; ...
             'Tau', {'u', 'v'}, 'Intensity', 'A'; ...
             'tau', {'phasorU', 'phasorV'}, 'photonCount', 'amp'; ...
             'mean', '', 'sum', 'mean'; ...
             'Avg. Fluorescence Lifetime [ns]', '', 'Cumulative Photon Count [photons]', 'Average Amplitude [photons]'};
[indx, tf] = listdlg('ListString', dataTypes(1, :), 'InitialValue', 1);
if ~tf
    return
end
dataTypes = dataTypes(:, indx);

% Set the objective magnification
mags = {'loMag x60', 'loMag x90', 'hiMag x60', 'hiMag x90'};
[indx, tf] = listdlg('ListString', mags, 'InitialValue', 1, ...
                     'SelectionMode', 'single');
if ~tf
    return
end
mags = mags{indx};

% Set the tau limits
prompt = {'Lower Tau Limit [ns]:', 'Upper Tau Limit [ns]:'};
dlgtitle = 'Choose Lifetime Limits';
dims = [1 35];
definput = {'0.1', '2'};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer)
    return
end
tauLim = str2double(answer);

% path
tauBins = tauLim(1) : 0.01 : tauLim(2);
tauMap = jet;
tauMap(1, :) = 0;

% Load the correct scaling factor for the microscope setting
umPerPixel = pixelSize.(mags(1 : 5)).(mags(7 : 9)).avg;
% Length of scale bar in um
scaleBar = 20;


% Load the FLIM image
%flimIm = readim(fullfile(path, fnames{1}));
% Expand to the full stack
%flimIm = repmat(flimIm, [1, 1, numel(fnames)]);


% Get the laser power in percentage
laserPerc = getLaserPerc(fnames{1});

% Expand to the list of laser percentages
laserPerc = repmat(laserPerc, [numel(fnames), 1]);

% Record the increase or decrase of laser power
downs = any(regexp(fnames{1}, 'down'));
% Expand to the list of laser power ups and downs
downs = repmat(downs, [numel(fnames), 1]);


for type = dataTypes
    % phasors come in two batches
    if iscell(type{2})
        for i = 1 : numel(type{2})
            stack.(type{3}{i}) = getImageFrame(path, fnames{1}, type{2}{i});
            stack.(type{3}{i}) = repmat(stack.(type{3}{i}), [1, 1, numel(fnames)]);
        end
    else
        % The rest is just one image
        stack.(type{3}) = getImageFrame(path, fnames{1}, type{2});
        stack.(type{3}) = repmat(stack.(type{3}), [1, 1, numel(fnames)]);
    end
end


[decays, tauHist, input] = getDecayHist(path, fnames{1}, stack.tau(:, :, 0), tauBins);
% Expand the list of decays
decays = repmat(decays, [numel(fnames), 1]);
% Expand the list of histograms
tauHist = repmat(tauHist, [numel(fnames), 1]);

[repetition, frameNr, frameLims] = getFrameTime(fnames{1});
% Expand the list of repetitions
repetition = repmat(repetition, [numel(fnames), 1]);
% Expand the list of frame numbers
frameNr = repmat(frameNr, [numel(fnames), 1]);
% Expand the list of frame limits
frameLims = repmat(frameLims, [numel(fnames), 1]);


% Import the phasor plot
%    phasor{1, i} = readim(fullfile(fnames(1).folder, [fnames(1).name(1 : end - 7), 'u.ics']));
%    phasor{2, i} = readim(fullfile(fnames(1).folder, [fnames(1).name(1 : end - 7), 'v.ics']));
% Expand to the full stack
%    phasor{1, i} = repmat(phasor{1, i}, [1, 1, numel(fnames)]);
%    phasor{2, i} = repmat(phasor{2, i}, [1, 1, numel(fnames)]);

% Import the amplitude images
%amplitude = readim(fullfile(fnames(1).folder, [fnames(1).name(1 : end - 7), 'A.ics']));
% Expand to the full stack
%amplitude = repmat(amplitude, [1, 1, numel(fnames)]);

% Import the photon count images
%photonCount = readim(fullfile(fnames(1).folder, [fnames(1).name(1 : end - 7), 'Intensity.ics']));
% Expand to the full stack
%photonCount = repmat(photonCount, [1, 1, numel(fnames)]);

for j = 2 : numel(fnames)
    % Read reminder of the slides
    for type = dataTypes
        % phasors come in two batches
        if iscell(type{2})
            for i = 1 : numel(type{2})
                stack.(type{3}{i})(:, :, j - 1) = ...
                    getImageFrame(path, fnames{j}, type{2}{i});
            end
        else
            % The rest is just one image
            stack.(type{3})(:, :, j - 1) = ...
                getImageFrame(path, fnames{j}, type{2});
        end
    end
    % Get the laser power percentage
    laserPerc(j) = getLaserPerc(fnames{j});
    downs(j) = any(regexp(fnames{j}, 'down'));
    % Load the decay and histogram
    [decays(j, :), tauHist(j, :)] = getDecayHist(path, fnames{j}, stack.tau(:, :, j - 1), tauBins);
    % Get the frame times
    [repetition(j), frameNr(j), frameLims(j, :)] = getFrameTime(fnames{j});
    % Import the phasor plot
%        phasor{1, i}(:, :, j - 1) = readim(fullfile(fnames(j).folder, [fnames(j).name(1 : end - 7), 'u.ics']));
%        phasor{2, i}(:, :, j - 1) = readim(fullfile(fnames(j).folder, [fnames(j).name(1 : end - 7), 'v.ics']));
    % Import the amplitude
%    amplitude(:, :, j - 1) = readim(fullfile(fnames(j).folder, [fnames(j).name(1 : end - 7), 'A.ics']));
    % Import the photon count
%    photonCount(:, :, j - 1) = readim(fullfile(fnames(j).folder, [fnames(j).name(1 : end - 7), 'Intensity.ics']));
end

%% This section is about sorting the files
% Check if the files are seperated in experimental repetitions (They have
% _cXX in their filename
if all(~isnan(repetition))
    % Looks like the files should be sorted by the _cXX numbers and frame
    % times
    [~, frameIndex] = sort(repetition + 0.001 * frameNr);
    % Time offset for each _cXX number
    timeOffset = zeros(size(repetition));
    % Get the time out
    unFrameRep = sort(unique(repetition(frameIndex)))';
    for i = 2 : numel(unFrameRep)
        for j = 1 : (i - 1)
            timeOffset(repetition == unFrameRep(i)) = ...
                timeOffset(repetition == unFrameRep(i)) + ...
                max(frameLims(repetition == unFrameRep(j), 2));
        end
    end
    frameTime = timeOffset + frameLims(:, 1);
elseif sum(downs) > 0
    % Sort the laser by power going up and down
    ups = find(~downs);
    downs = find(downs);
    % Laser power increasing index
    [~, upIndex] = sort(laserPerc(ups));
    upIndex = ups(upIndex);
    % Laser power decreasing index
    [~, downIndex] = sort(laserPerc(downs), 'descend');
    downIndex = downs(downIndex);
    % Sort all the variables according to the up and down indices
    frameIndex = [upIndex, downIndex];
    % The frame time in this instance
    frameTime = zeros(size(repetition));
else
    % If no clear sorting mechanism is in place, just sort them 1 to end
    frameIndex = 1 : numel(fnames);
    % The frame time in this instance
    frameTime = zeros(size(repetition));
end
% sort all the inputs
fnames = fnames(frameIndex);
laserPerc = laserPerc(frameIndex);
%flimIm = flimIm(:, :, frameIndex - 1);
decays = decays(frameIndex, :);
tauHist = tauHist(frameIndex, :);
repetition = repetition(frameIndex);
frameNr = frameNr(frameIndex);
frameLims = frameLims(frameIndex, :);
frameTime = frameTime(frameIndex);

for type = dataTypes
    % phasors come in two batches
    if iscell(type{2})
        for i = 1 : numel(type{2})
            stack.(type{3}{i}) = stack.(type{3}{i})(:, :, frameIndex - 1);
        end
    else
        stack.(type{3}) = stack.(type{3})(:, :, frameIndex - 1);
    end
end
%    phasor{1, i} = phasor{1, i}(:, :, upDownIndex - 1);
%    phasor{2, i} = phasor{2, i}(:, :, upDownIndex - 1);
%amplitude = amplitude(:, :, frameIndex - 1);
%photonCount = photonCount(:, :, frameIndex - 1);
% Calculate the PPFD
%PPFD = polyval(model.PPFD, laserPerc);
PPFD = model.PPFD_interp(2, arrayfun(@(x)(find(model.PPFD_interp(1, :) == x)), laserPerc));

%% Run the animation
figSize = size(flimIm);
hf = figure('Units', 'pixels', 'Position', [10, 10, 2 * figSize(1) + 150, figSize(2) + 150]);
% Tau image figure
ha(1) = axes('Units', 'pixels', 'Position', [50, 100, figSize([1, 2])]);
him = imagesc(double(stack.tau(:, :, 0)));
colormap(tauMap)
set(ha(1), 'Box', 'on', 'XTick', [], 'YTick', [])
% Scale bar
hPPFD = text(figSize(1), figSize(2), ...
                sprintf('PPFD = %d %smol%sm%ss%s ', ...
                        round(PPFD(1)), char(956), char(183), ...
                        char([8315, 178, 183]), char([8315, 185])), ...
                'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w', ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'right');
rcet = rectangle('Position', ...
                    [20, 20, (scaleBar / umPerPixel), 5], ...
                    'FaceColor', 'w', 'EdgeColor', 'w');
scText = text(20 + scaleBar / 2 / umPerPixel, 28, ...
                sprintf('%d %sm', scaleBar, char(956)), ...
                'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w', ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center');
% Time
hTime = text(0, figSize(2), ...
                datestr(datenum(num2str(frameTime(1)), 'SS'), ' HH:MM'), ...
                'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w', ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left');
% Colorbar axes
ha(2) = axes('Units', 'pixels', 'Position', [50, 50, figSize(1), 50]);
imagesc(repmat(tauLim(1) : 0.001 : tauLim(2), 2, 1))
set(ha(2), 'XTick', [], 'YTick', [])
text(min(get(ha(2), 'Xlim')), mean(get(ha(2), 'Ylim')), ...
     sprintf(' %g ns', tauLim(1)), 'Color', 'w', 'FontSize', 14, ...
     'VerticalAlignment', 'middle', 'FontWeight', 'bold')
text(max(get(ha(2), 'Xlim')), mean(get(ha(2), 'Ylim')), ...
     sprintf('%g ns ', tauLim(2)), 'Color', 'w', 'FontSize', 14, ...
     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
     'FontWeight', 'bold')
% Decay axes
ha(3) = axes('Units', 'pixels', 'Position', [figSize(1) + 150, 100, figSize(1) - 50, figSize(2) / 2 - 30]);
hdecay = plot(input.tac - input.tac(1), decays(1, :), 'k.');
set(ha(3), 'XLim', [0, 20], 'YScale', 'log', 'FontSize', 14, ...
    'YLim', [1, 10 ^ ceil(log10(max(decays(:))))], ...
    'YLimMode', 'manual', 'YTick', 10 .^ (0 : 7))
text(max(get(ha(3), 'XLim')), max(get(ha(3), 'YLim')), ...
     'Cumulative Decay ', 'FontSize', 14, ...
     'FontWeight', 'bold', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top')
xlabel('Time [ns]')
ylabel('Photon Count')
% Histogram axes
ha(4) = axes('Units', 'pixels', 'Position', [figSize(1) + 150, figSize(2) / 2 + 130, figSize(1) - 50, figSize(2) / 2 - 30]);
hhist= plot(tauBins(1 : end - 1) + diff(tauBins([1, 2])) / 2, ...
            1.1 * tauHist(max(max(tauHist, [], 2)) == max(tauHist, [], 2), :), ...
            'k-', 'LineWidth', 2);
drawnow
set(ha(4), 'XLim', tauLim, 'FontSize', 14, 'YLimMode', 'manual')
text(max(get(ha(4), 'XLim')), max(get(ha(4), 'YLim')), ...
     'Fluorescence Lifetime Density ', 'FontSize', 14, ...
     'FontWeight', 'bold', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top')
xlabel('Fluorescene Lifetime [ns]')
ylabel('Pixel Count')
box on
tmppos = get(ha(4), 'Position');

videoFname = regexp(path, filesep, 'split');
videoFname = videoFname{find(~cellfun(@isempty, videoFname), 1, 'last')};
% If there are phasor data to plot
if any(strcmp(dataTypes(:), 'phasorU'))
    ha(5) = axes('Units', 'pixels', ...
                 'Position', tmppos, ...
                 'YLim', [0, 0.6], 'YLimMode', 'manual', ...
                 'XLim', [-1, 1] * 0.6 / tmppos(4) * tmppos(3) / 2, ...
                 'XLimMode', 'manual', 'FontSize', 14, 'Box', 'on', ...
                 'XTick', [], 'YTick', []);
    xlabel('u')
    ylabel('v')
    hold on
    % Plot the semi circle
    plot(0.5 * cos(0 : pi / 30 : pi), ...
         0.5 * sin(0 : pi / 30 : pi), 'LineWidth', 2)
    u = double(stack.phasorU(:, :, 0) - 0.5);
    u = u(u ~= -0.5);
    v = double(stack.phasorV(:, :, 0));
    v = v(v ~= 0);
    dscatter(u, v);
    text(max(get(ha(i, 5), 'XLim')), max(get(ha(i, 5), 'YLim')), ...
         'Phasor Plot ', 'FontSize', 14, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'top')
    vidObj2 = VideoWriter(fullfile(path, [videoFname, '.phasor.avi']));
    vidObj2.FrameRate = 1;
    open(vidObj2);
end

vidObj1 = VideoWriter(fullfile(path, [videoFname, '.tauHist.avi']));
vidObj1.FrameRate = 1;
open(vidObj1);
for j = 1 : numel(fnames)
    set(him, 'CData', double(stack.tau(:, :, j - 1)))
    set(ha(1), 'CLim', tauLim)
    hPPFD.String = sprintf('PPFD = %d %smol%sm%ss%s ', ...
                           round(PPFD(j)), char(956), ...
                           char(183), char([8315, 178, 183]), ...
                           char([8315, 185]));
    % Update time
    hTime.String = datestr(datenum(num2str(frameTime(j)), 'SS'), ' HH:MM');
    set(hdecay, 'YData', decays(j, :))
    set(hhist, 'YData', tauHist(j, :))
    if any(strcmp(dataTypes(:), 'phasorU'))
        % Hide phasor axes
        set(ha(5), 'Visible', 'off')
        set(get(ha(5), 'Children'), 'Visible', 'off')
        % Show histogram axes
        set(ha(4), 'Visible', 'on')
        set(get(ha(4), 'Children'), 'Visible', 'on')
    end
    drawnow
    % Save frame of histogram axes
    currFrame = getframe(hf);
    writeVideo(vidObj1, currFrame);
    if any(strcmp(dataTypes(:), 'phasorU'))
        % Hide histogram axes
        set(ha(4), 'Visible', 'off')
        set(get(ha(4), 'Children'), 'Visible', 'off')
        % Show phasor axes
        set(ha(5), 'Visible', 'on')
        set(get(ha(5), 'Children'), 'Visible', 'on')
        axes(ha(5))
        % Update the scatter
        delete(findobj(ha(5), 'type', 'Scatter'))
        u = double(stack.phasorU(:, :, 0) - 0.5);
        u = u(u ~= -0.5);
        v = double(stack.phasorV(:, :, 0));
        v = v(v ~= 0);
        dscatter(u, v);
        drawnow
        % Save frame of phasor axes
        currFrame = getframe(hf);
        writeVideo(vidObj2, currFrame);
    end
end
close(vidObj1);
if any(strcmp(dataTypes(:), 'phasorU'))
    close(vidObj2);
end

%% Make a graph of lifetime evolution
% Convert the image stack to doubles, ignoring the first two frames
for type = dataTypes
    % phasors come in two batches
    if ~iscell(type{2})
        stack.avg.(type{3}) = double(stack.(type{3}));
    end
    stack.avg.(type{3})(stack.avg.tau == 0) = NaN;
    if ~isempty(type{4})
        stack.avg.(type{3}) = eval([type{4}, '(reshape(stack.avg.(type{3}), size(stack.avg.(type{3}), 1) * size(stack.avg.(type{3}), 2), size(stack.avg.(type{3}), 3)), ''omitnan'')']);
    end

% pcAvg = double(photonCount(:, :, 2 : end));
% tauAvg = double(flimIm(:, :, 2 : end));
% pcAvg(tauAvg == 0) = NaN;
% tauAvg(tauAvg == 0) = NaN;
% pcAvg = sum(reshape(pcAvg, size(pcAvg, 1) * size(pcAvg, 2), size(pcAvg, 3)), 'omitnan');
% tauAvg = mean(reshape(tauAvg, size(tauAvg, 1) * size(tauAvg, 2), size(tauAvg, 3)), 'omitnan');
% timeX = (0 : (numel(tauAvg) - 1)) * 60;
    if strcmp(type{3}, 'tau')
        figure;
        yyaxis left
        plot(frameTime, stack.avg.tau, 'LineWidth', 2);
        xlabel('Time [s]');
        ylabel(type{5});
        set(gca, 'FontSize', 14)
        text(0.05 * max(get(gca, 'XLim')), ...
             diff(get(gca, 'YLim')) * 0.95 + min(get(gca, 'YLim')), ...
             hPPFD.String, 'FontSize', 14, 'FontWeight', 'bold', ...
             'VerticalAlignment', 'top')
    end
    if strcmp(type{3}, 'photonCount') || strcmp(type{3}, 'amp')
        yyaxis right
        plot(frameTime, stack.avg.(type{3}), 'LineWidth', 2);
        ylabel(type{5});
        print(gcf, fullfile(path, [videoFname, '.', type{3}, '.png']), '-dpng')
    end
end


function laserPerc = getLaserPerc(fname)
    percString = cell2mat(regexp(fname, '_[0-9][0-9]-[0-9]perc', 'match'));
    if isempty(percString)
        laserPerc = NaN;
    else
        laserPerc = ...
            sum(cell2mat(textscan(percString, '_%f-%f-perc')) .* [1, 0.1]);
    end
end

function imageOut = getImageFrame(path, fname, repString)
    imageOut = readim(fullfile(path, ...
        regexprep(fname, '_Tau.ics', ['_' repString '.ics'])));
end

function [decay, tauHist, input] = getDecayHist(path, fname, flimIm, tauBins)
    % Load the decays
    load(fullfile(path, [fname(1 : end - 18), '.mat']), 'XYZimage', 'input')
    % Calculate the decay
    decay = squeeze(sum(sum(XYZimage .* uint16((flimIm > 0))', 1), 2))';

    % Calculate Tau histograms
    tauHist = histcounts(double(flimIm), tauBins);
end

function [repetition, frameNr, frameLims] = getFrameTime(fname)
    repString = cell2mat(regexp(fname, '_c[0-9][0-9]\.', 'match'));
    if isempty(repString)
        repetition = NaN;
    else
        repetition = cell2mat(textscan(repString, '_c%f.'));
    end
    nrString = cell2mat(regexp(fname, '\.[0-9][0-9][0-9]_', 'match'));
    if isempty(nrString)
        frameNr = NaN;
    else
        frameNr = cell2mat(textscan(nrString, '.%f_'));
    end
    limString = cell2mat(regexp(fname, '_[0-9]+s-[0-9]+s_', 'match'));
    if isempty(nrString)
        frameLims = [NaN, NaN];
    else
        frameLims = cell2mat(textscan(limString, '_%fs-%fs_'));
    end
end