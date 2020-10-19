
% Set some global variables
global XYZimage % image stack
global tac      % TAC bin range
global ftac     % Fitted TAC range
global fParam   % Fitted parameters
global fImage   % Fitted curve image
global rImage   % Residuals image
global tauImage % Fluorescence lifetime image
global data     % dummy variable for exportICS function
global fifo     % dummy variable for exportICS function

% Remove the old fifo variables
clear FIFO
close all

tic

%% Load FIFO data
% Create a wildcard for looking for the right type of files
wildcard = sprintf('/home/k1216365/experiments/%s/', ...
                   datetime('now', 'Format', 'yyyy'));
% Create a dialog box to load the file
[filename, pathname] = ...
    uigetfile({'*m1.spc', 'SPC-files (*.m1.spc)'; ...
               '*.*',     'All Files (*.*)'}, ...
               'Pick a file', ...
               'MultiSelect', 'off', ...
               wildcard);
% Check if the existing file is not provided
if isequal(filename, 0)
    return
end

% Build a cell of expected filenames
fnames = cell(1, 3);
for i = 1 : numel(fnames)
    filename(end - 4) = num2str(i);
    fnames{i} = fullfile(pathname, filename);
    assert(exist(fnames{i}, 'file') == 2, ...
           'File %s does not exist.', fnames{i})
end

for i = 1 : numel(fnames)
    FIFO(i) = SPC830read(fnames{i}); %#ok<SAGROW>
end

%% Expand the second channel of the FIFO data.
%  This is required because the FIFO macrotime from the different cards
%  shows relative offset to each other. This means that not simple
%  macrotime overlaps happen all the time, instead there are variable 0, 1,
%  2, ... counts off

%  Shift of the macrotime
shift = 2;

%  Create a matrix of fifo(2) macrotimes that is up to a few (shift) values
%  smaller and larger than the origininal one to look for instersections
%  with the other fifo macrotimes.
mT1 = repmat(FIFO(1).macroT, 2 * shift + 1, 1) + (-shift : shift)';

%  create an index lookup table for the expanded macrotime
%[ix, iy] = meshgrid(1 : size(mT2, 2), 1 : 2 * shift + 1);


%  Intersect the expanded fifo(1) macrotime matrix with the fifo(2)
%  macrotime matrix
[C, index12, index2] = compareEvents(mT1, FIFO(2).macroT);

%  Use values from the fifo(1) expanded matrix that have been found in the
%  intersection with the fifo(2) matrix and intersect them with the fifo(3)
%  macrotime matrix
[~, index13, I3] = compareEvents(mT1(:, index12), FIFO(3).macroT);

%  Get the indices of the coincident events for fifo(1) and fifo(2)
I1 = index12(index13);
I2 = index2(index13);


%% Create the image
% bins to combine the photons
bins = 0 : 8 : 4096;
% Create an XY image from the data in fifo(1) and fifo(2)
% Actually don't create the image, just get the indices of the bin elements
% for the input arrays
[~, ~, ~, binX, binY] = ...
    histcounts2(FIFO(2).microT(I2), FIFO(3).microT(I3), bins, bins);


% Create a 3D array to hold the time decay data. Use 16-bit unsigned
% integer
XYZimage = zeros(numel(bins) * [1 1 1], 'uint16');

% Convert the microtime into an addressable bin coordinate
binZ = bitshift(FIFO(1).microT(I1), -log2(bins(2))) + 1;
% TAC bin range in nanoseconds
tac = linspace(0, FIFO(1).set.SP.SP_TAC_R * 1e9, numel(bins) - 1);

% Run through every photon and throw it into the right 3D histogram bin
for i = 1 : numel(binX)
    XYZimage(binY(i), binX(i), binZ(i)) = ...
        XYZimage(binY(i), binX(i), binZ(i)) + 1;
end

%% Run the lifetime fit for pixels with more that 10 photons
%  reshape XYZ image into an array where the first dimension is the time
%  and the other dimensions are all the pixels
TPimage = reshape(XYZimage, ...
                  size(XYZimage, 1) * size(XYZimage, 2), ...
                  size(XYZimage, 3))';
TPimage = double(TPimage);
fit_index = find(sum(TPimage) > 10);

% Get the x increment
xincr = tac(end) / numel(bins);             % nanoseconds per bin
% Create a simulated prompt as we don't have one currently.
promptSigma = 0.05;
prompt = histc(tac(end) / 2 + promptSigma * randn(10000, 1), tac);
prompt = prompt / sum(prompt);
prompt = prompt(prompt > max(prompt) / 100);

% Get the start, where the transient starts rising
start = 195;
% Get the fit_start, where the fitting should start
fit_start = 6;
% Get the fit_end, where the fit should end
fit_end = size(TPimage, 1) - start - 2;

[paramsLMA, paramsRLD, fittedLMA] = ...
    mxSlimCurve(TPimage(start : end, fit_index), prompt, xincr, fit_start);

% Store the fitted results
fImage = zeros(size(fittedLMA, 1), size(TPimage, 2));
fImage(:, fit_index) = fittedLMA;
fImage = reshape(fImage', [size(XYZimage, 1), size(XYZimage, 2), size(fImage, 1)]);

% Create an image of residuals
rImage = abs((fImage - double(XYZimage(:, :, start : end))) ./ sqrt(fImage));

% Store the fitted parameters
fParam = zeros(size(paramsLMA, 1), size(TPimage, 2));
fParam(:, fit_index) = paramsLMA;
fParam = reshape(fParam', [size(XYZimage, 1), size(XYZimage, 2), size(fParam, 1)]);
% Add one to the chi^2 value, to make it realistic
fParam(:, :, end) = 1 + fParam(:, :, end);

% Store the fitted axis
ftac = tac(start : end);

% Generate lifetime image
tauImage = squeeze(fParam(:, :, 3));
tauImage(tauImage > 10) = NaN;

% Export the ICS file
ICSimage = permute(flipud(XYZimage(1 : end - 1, 1 : end - 1, 1 : end - 1)), [2, 3, 1]);
data.square = size(XYZimage, 1) - 1;
fifo = FIFO(1);
% Create a filename for the ICS file
fnames{4} = regexprep(fnames{1}, '_M1.SPC', '.ics', 'ignorecase');
exportICS(ICSimage, fnames{4})


%% Display the image
renderResults