
global ha       % axes handles
global hl       % line handles
global hu       % uicontrol handles
global tac      
global fParam
global tauImage % Fluorescence lifetime image
global stopped  % Cursor stopped flag
global hSpin    % Spinner handle
global tauLim   % Limits of lifetime
global tauCLim  % Limits of lifetime colormap


%% Display the image
close all
% Get screen size
sSize = get(0, 'ScreenSize');
% Prepare figure size
fSize = [(sSize(3) - sSize(4)) / 2, sSize(4) * 0.1, sSize(4) * [1, 0.8]];
% Prepare axes size
aSize = [fSize(3) * 0.1, fSize(3) * 0.05, ...
         min([0.4 * fSize([4, 4]); 0.8 * fSize([3, 3])])];
% Create a figure
hf = figure('Units', 'pixels', 'Position', fSize);

% Create a pointer for the figure
set(hf, 'PointerShapeCData', createPointer, ...
        'PointerShapeHotSpot', [16 16]);

% Add some text elements
pos = [aSize(1), 1.5 * aSize(2) + aSize(4), 40, 20];
hu.text.XpixLab = uicontrol('Style', 'text', ...
                            'String', 'X:', ...
                            'Position', pos + [0, 30, 0, 0]);
hu.text.YpixLab = uicontrol('Style', 'text', ...
                            'String', 'Y:', ...
                            'Position', pos + [0, 0, 0, 0]);
hu.text.XpixVal = uicontrol('Style', 'text', ...
                            'String', '', ...
                            'Position', pos + [pos(3) + 10, 30, 0, 0], ...
                            'BackgroundColor', [1 1 1]);
hu.text.YpixVal = uicontrol('Style', 'text', ...
                            'String', '', ...
                            'Position', pos + [pos(3) + 10, 0, 0, 0], ...
                            'BackgroundColor', [1 1 1]);

% Create titles for lifetime fits
names = cell(1, size(fParam, 3));
order = zeros(size(names));
names{1} = 'Z:';
order(1) = numel(order) - 1;
names{end} = [char([967, 178]), ':'];
order(end) = numel(order);
if numel(names) == 4
    names{2} = 'A:';
    order(2) = 1;
    names{3} = [char(964), ':'];
    order(3) = 2;
else
    for i = 0 : 2 : numel(names) - 3
        names{i + 2} = ['A' char(8321 + i / 2) ':'];
        order(i + 2) = 1 + i;
        names{i + 3} = [char([964, 8321 + i / 2]) ':'];
        order(i + 3) = 2 + i;
    end
end

% Create starting position for lifetime fits
pos = [aSize(1) + 120, 1.5 * aSize(2) + aSize(4), 40, 20];

for i = 1 : numel(names)
    % i is the counting index
    tpos = pos + [0, 30 * (numel(names) - order(i)), 0, 0];
    hu.text.fitLab(i) = uicontrol('Style', 'text', ...
                                  'String', names(i), ...
                                  'Position', tpos);
    tpos = tpos + [tpos(3) + 10, 0, 0, 0];
    hu.text.fitVal(i) = uicontrol('Style', 'text', ...
                                  'String', num2str(i), ...
                                  'Position', tpos, ...
                                  'BackgroundColor', [1 1 1]);
end

%% Place four axes there
% Axes for intensity image
ha(1) = axes('Units', 'pixels', 'Position', aSize);
colormap(ha(1), 'gray')

% Axes for lifetime fit
aSize(1) = aSize(1) + fSize(3) / 2;
aSize(4) = aSize(4) * 0.6;
ha(2) = axes('Units', 'pixels', ...
             'Position', aSize, ...
             'XLim', tac([1, end]));
% Axes for residuals
aSize(2) = aSize(2) + aSize(4) * 7 / 6;
aSize(4) = aSize(4) * 0.5;
ha(3) = axes('Units', 'pixels', ...
             'Position', aSize, ...
             'XLim', tac([1, end]));

% Axes for lifetime image
ha(4) = axes('Units', 'pixels');
aSize(4) = aSize(4) / 0.5 / 0.6;
aSize(2) = aSize(2) + aSize(4) * 0.5 - 20;
set(ha(4), 'Position', aSize);

% Draw the maximum projection
axes(ha(1))
image1 = flipud(sum(XYZimage, 3));
hi(1) = imagesc(image1);
hold on
% Flip the Y axis, add the box
set(ha(1), 'YDir', 'reverse', 'Box', 'on', 'XAxisLocation', 'top')

% Draw two four dummy lines later used as a cross hair
for i = 1 : 4
    hl.ax1(i) = plot([0 0], [0 0], 'g-');
end

% Draw the maximum projection
axes(ha(2))
hold on
% draw two dummy lines on axes 2
hl.ax2(1) = plot([0 0], [0 0], 'b-');
hl.ax2(2) = plot([0 0], [0 0], 'r--');

% Add the X- and Y-axis labels
xlabel('TAC Microtime [ns]')
ylabel('Photon Count')

% Draw a dummy line of residuals
axes(ha(3))
hl.ax3(1) = plot([0 0], [0 0], 'b.');

% Draw the lifetime map
axes(ha(4))
hi(2) = imagesc(flipud(tauImage));
pos1 = get(ha(4), 'Position');
hold on
cmap = parula(256);
cmap(end, :) = [1 0 0]; % Make highest value red
cmap(1, :) = [0 0 0];   % Make lowest value black
colormap(ha(4), cmap)
ha(5) = colorbar('westoutside');
set(ha(4), 'Position', pos1)
set(ha(5), 'Units', 'Pixels')
pos2 = get(ha(5), 'Position');
pos2(1) = pos1(1) - pos1(3) * 0.2;
set(ha(5), 'Position', pos2);
% Label the colorbar
set(get(ha(5), 'Label'), 'String', 'Fluorescence Lifetime [ns]')

% Change the colormap to grayscale
% Flip the Y axis, add the box
set(ha(4), 'YDir', 'reverse', ...
           'Box', 'on', ...
           'XAxisLocation', 'top', ...
           'XLim', get(ha(1), 'XLim'), ...
           'YLim', get(ha(1), 'YLim'))



% Draw two four dummy lines later used as a cross hair
for i = 1 : 4
    hl.ax4(i) = plot([0 0], [0 0], 'g-');
end



%% Add lifetime scale spinners
pos3 = [pos2(1) - pos2(3) * 3, pos2(2), 60, 26];
hSpin(1) = spinner('Parent', gcf, ...
                   'Position', pos3, ...
                   'Callback', {@spinner_Callback, gcf}, ...
                   'Tag', 'lowTau');
pos3(2) = pos3(2) + pos2(4) - 26;
hSpin(2) = spinner('Parent', gcf, ...
                   'Position', pos3, ...
                   'Callback', {@spinner_Callback, gcf},...
                   'Tag', 'highTau');
% Set the spinners
tauLim = round(20 * [min(tauImage(:)), max(tauImage(:))]) / 20;
tauCLim = tauLim;

set(hSpin, 'SliderStep', [0.05, 0.05], ...
           'Min', tauLim(1), ...
           'Max', tauLim(2))
set(hSpin(1), 'Value', tauLim(1), ...
              'String', num2str(tauLim(1)))
set(hSpin(2), 'Value', tauLim(2), ...
              'String', num2str(tauLim(2)))
          
% Make the figure responsive to the mouse
% Call 'mousemove' function when mouse moves over the window
set(hf, 'WindowButtonMotionFcn', @mousemove);

% Call 'mouseclick' function when mouse is clicked over the window
set(hf, 'WindowButtonDownFcn', @mouseclick);
stopped = false;

%% Save the maximum projection image as a 16-bit tiff
fnames{5} = [fnames{4}(1 : end - 3), 'tiff'];
imwrite(uint16(image1), fnames{5})

function spinner_Callback(~)
global hSpin
global ha
global tauCLim

% get spinner data for the spinners
sdata = cell2mat(get(hSpin, 'Value'));

% Check which of the two spinners has changed
if sdata(1) == tauCLim(1) && sdata(2) == tauCLim(2)
    % If none has changed, just leave
    return
end

% Check if sdata1 is smaller than sdata2
if sdata(1) >= sdata(2)
    % Set them to the original values
    sdata(1) = tauCLim(1);
    sdata(2) = tauCLim(2);
    set(hSpin(1), 'Value', sdata(1), ...
                  'String', num2str(sdata(1)))
    set(hSpin(2), 'Value', sdata(2), ...
                  'String', num2str(sdata(2)))
end

tauCLim = sdata;
set(ha(4), 'CLim', tauCLim)
end


function mouseclick(~, ~)
global ha
global stopped

% Check if the mouse is in either the 1st or 4th axes
inAxes = false;
% Retrieve the coordinate of the mouse in either the 1st or 4th axes
for i = [1 4]
    coord = get(ha(i), 'Currentpoint');
    % Keep only the useful bit
    coord = coord(1, [1 2]);
    % Retrieve the limits of the axes
    xlim = get(ha(i), 'Xlim');
    ylim = get(ha(i), 'Ylim');
    % Check if pointer is outside of the limits
    outX = any(diff([xlim(1) coord(1) xlim(2)]) < 0);
    outY = any(diff([ylim(1) coord(2) ylim(2)]) < 0);
    if ~any([outX, outY])
        inAxes = true;
    end
end
% End the function if it is outside in any of the directions
if ~inAxes
    return
end

stopped = ~stopped;

mousemove(false, false);
end


function mousemove(~, ~)
% Use the global variable for the axes and line handles
global ha
global hl
global XYZimage
global hu
global tac
global ftac
global fImage
global fParam
global rImage
global stopped


% Check if the mouse is in either the 1st or 4th axes
inAxes = false;
% Retrieve the coordinate of the mouse in either the 1st or 4th axes
for i = [1 4]
    coord = get(ha(i), 'Currentpoint');
    % Keep only the useful bit
    coord = coord(1, [1 2]);
    % Retrieve the limits of the axes
    xlim = get(ha(i), 'Xlim');
    ylim = get(ha(i), 'Ylim');
    % Check if pointer is outside of the limits
    outX = any(diff([xlim(1) coord(1) xlim(2)]) < 0);
    outY = any(diff([ylim(1) coord(2) ylim(2)]) < 0);
    if ~any([outX, outY])
        inAxes = true;
        C = coord;
    end
end
% End the function if it is outside in any of the directions
if ~inAxes
    set(gcf, 'Pointer', 'arrow')
    return
end

% Set a camera focus pointer
set(gcf, 'Pointer', 'custom')


% If the motion is stopped, just exit the function
if stopped
    return
end

% Draw the cursor
set(hl.ax1(1), 'Xdata', C([1 1]));
set(hl.ax1(1), 'Ydata', [C(2) + 0.05 * diff(ylim), ylim(2)]);
set(hl.ax1(2), 'Xdata', C([1 1]));
set(hl.ax1(2), 'Ydata', [ylim(1), C(2) - 0.05 * diff(ylim)]);
set(hl.ax1(3), 'Xdata', [C(1) + 0.05 * diff(xlim), xlim(2)]);
set(hl.ax1(3), 'Ydata', C([2 2]));
set(hl.ax1(4), 'Xdata', [xlim(1), C(1) - 0.05 * diff(xlim)]);
set(hl.ax1(4), 'Ydata', C([2 2]));

for i = 1 : 4
    set(hl.ax4(i), 'XData', get(hl.ax1(i), 'XData'), ...
                   'YData', get(hl.ax1(i), 'YData'))
end
% Round the pixel number, subtract 1 to start from zero
C = round(C);

% Display the X and Y pixel positions
set(hu.text.XpixVal, 'String', num2str(C(1) - 1));
set(hu.text.YpixVal, 'String', num2str(C(2) - 1));

% Calculate the image row and column index based on the mouse position
rowIn = size(XYZimage, 1) + 1 - C(2);
colIn = C(1);
% Display the fit values
for i = 1 : size(fParam, 3)
    set(hu.text.fitVal(i), 'String', num2str(fParam(rowIn, colIn, i)))
end

% Draw the decay and fit
set(hl.ax2(1), 'XData', tac, ...
               'YData', squeeze(XYZimage(rowIn, colIn, 1 : end - 1)))
set(hl.ax2(2), 'XData', ftac, ...
               'YData', squeeze(fImage(rowIn, colIn, 1 : end - 1)))

set(ha(2), 'XLim', tac([1, end]))

% Draw the residuals
set(hl.ax3(1), 'XData', ftac, ...
               'YData', squeeze(rImage(rowIn, colIn, 1 : end - 1)))

set(ha(3), 'XLim', tac([1, end]))


end