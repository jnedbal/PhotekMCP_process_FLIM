function [XYZimage, input, intensityImage] = ...
    SPC2image(filename, MCPshiftFile, input)

% Create an empty filename if it does not exist
if nargin == 0
    filename = '';
end

% Create an empty IRF filename if it does not exist
if nargin < 2
    MCPshiftFile = '';
end

% Create an empty struct if no input is supplied
if nargin < 3
    input = struct('Xstart', 0);
end

% Check the input is a struct
assert(isstruct(input), 'Input needs to be a struct');

%  Shift of the macrotime
if ~isfield(input, 'shift')
    input.shift = 3;
end

% First coordinate for the X axis
if ~isfield(input, 'Xstart')
    input.Xstart = 0;
end

% Last coordinate for the X axis
if ~isfield(input, 'Xend')
    input.Xend = 4096;
end

% First coordinate for the Y axis
if ~isfield(input, 'Ystart')
    input.Ystart = 0;
end

% Last coordinate for the Y axis
if ~isfield(input, 'Yend')
    input.Yend = 4096;
end

% Spatial binning
if ~isfield(input, 'XYstep')
    input.XYstep = 8;
end

% First time bin on the microtime axis
if ~isfield(input, 'Tstart')
    input.Tstart = 0;
end

% Last time bin on the microtime axis
if ~isfield(input, 'Tend')
    input.Tend = 4096;
end

% Time axis binning
if ~isfield(input, 'Tstep')
    input.Tstep = 8;
end

% Frame integration time (Make it Inf for one frame only)
if ~isfield(input, 'frameTime')
    input.frameTime = Inf;
end

% File name with the source data
if ~isfield(input, 'filename') && ~exist(filename, 'file')
    % Create a wildcard for looking for the right type of files
    wildcard = sprintf('/home/jakub/experiments/%s/', ...
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
    input.filename = fullfile(pathname, filename);
else
    input.filename = filename;
end

% Get the source path name
input.pathname = fileparts(input.filename);

% Add the IRFfile if it exists
input.MCPshiftFile = MCPshiftFile;

% Add the randomization weighing factor for MCP shift correction.
if ~isfield(input, 'MCPrandomWeight')
    input.MCPrandomWeight = 0.5;
end

% IRF shift model
if ~isfield(input, 'MCPshiftModel')
    input.MCPshiftModel = 'poly55';
end

% Load the MCP shift dataset
if ~isempty(input.MCPshiftFile)
    load(input.MCPshiftFile, 'MCPshift')
end

% Should the stack be saved to a MAT file?
if ~isfield(input, 'saveMAT')
    input.saveMAT = true;
end

% Should the stack be saved to an ICS file?
if ~isfield(input, 'saveICS')
    input.saveICS = true;
end

%% Round the bin coordinates
input.Xstart = input.XYstep * round(input.Xstart / input.XYstep);
input.Xend = input.XYstep * round(input.Xend / input.XYstep);
input.Ystart = input.XYstep * round(input.Ystart / input.XYstep);
input.Yend = input.XYstep * round(input.Yend / input.XYstep);
input.Tstart = input.Tstep * round(input.Tstart / input.Tstep);
input.Tend = input.Tstep * round(input.Tend / input.Tstep);

% Spatial bins to combine the photons
input.Xbins = input.Xstart : input.XYstep : input.Xend;
input.Ybins = input.Ystart : input.XYstep : input.Yend;
% Spatial bins to combine the photons
input.Tbins = input.Tstart : input.Tstep : input.Tend;


%% Load the source files
% Build a cell of expected filenames
fnames = repmat({input.filename}, 1, 3);
for i = 1 : numel(fnames)
    fnames{i}(end - 4) = num2str(i);
    assert(exist(fnames{i}, 'file') == 2, ...
           'File %s does not exist.', fnames{i})
end


for i = 1 : numel(fnames)
    FIFO(i) = SPC830read(fnames{i}); %#ok<AGROW>
end

% Sometimes the B&H software makes empty one or more SPC files. This make
% further analysis impossible
% Check if the gaps are in the first bin
if any([FIFO(1).gap(1), FIFO(2).gap(1), FIFO(3).gap(1)])
    % It looks like one or more of the TCSPC cards were turned off
    % Finish without doing anything
    warndlg('Empty SPC file. Finishing early.', 'SPC2image', 'modal');
    % Create empty output
    XYZimage = NaN;
    intensityImage = NaN(numel(input.Ybins), numel(input.Xbins));
end

% Look for gaps in the FIFO data stream and try to reconsittute the
% macrotime data stream after the gaps
SPC2imageGapCorr;

% Get rid of data that is not needed to save memory
FIFO = rmfield(FIFO, {'rout', 'gap', 'mtov'});

% TAC bin size in seconds
binWidth = FIFO(1).set.SP.SP_TAC_TC * input.Tstep;
% TAC bin range in nanoseconds
input.tac = input.Tbins * binWidth * 1e9;


%% Go through the macrotimes to create individual files
frameBlocks = 0 : input.frameTime : FIFO(1).macroT(end) * FIFO(1).MTC;
% Assemble output filename list
input.outputStackName = repmat({input.filename}, numel(frameBlocks), 2);
for i = 1 : numel(frameBlocks)
    input.outputStackName{i, 1} = ...
        [input.outputStackName{i}(1 : end - 7), ...
         sprintf('.%03d_%ds-%ds.mat', ...
                 i - 1, ...
                 round(frameBlocks(i)), ...
                 round(min(frameBlocks(i) + input.frameTime, ...
                           FIFO(1).macroT(end) * FIFO(1).MTC)))];
	input.outputStackName{i, 2} = ...
        [input.outputStackName{i, 1}(1 : end - 3), 'ics'];
end
% Create a matrix of intensity images
if nargout > 2
    intensityImage = ...
        zeros(numel(input.Ybins), numel(input.Xbins), numel(frameBlocks));
end
% Assemble a 3D matrix of individual frames
for frame = 1 : numel(frameBlocks)
    %% Expand the second channel of the FIFO data.
    %  This is required because the FIFO macrotime from the different cards
    %  shows relative offset to each other. This means that not simple
    %  macrotime overlaps happen all the time, instead there are variable 
    %  0, 1, 2, ... counts off



    %  Index of macrotimes in the existing frame
    frameIndex = frameBlocks(frame) / FIFO(1).MTC <= FIFO(1).macroT & ...
                 FIFO(1).macroT < (frameBlocks(frame) + input.frameTime) / FIFO(1).MTC;

    %  Create a matrix of fifo(2) macrotimes that is up to a few (shift) 
    %  values smaller and larger than the origininal one to look for
    %  intersections with the other fifo macrotimes.
    mT1 = repmat(FIFO(1).macroT(frameIndex), 2 * input.shift + 1, 1) + ...
          (-input.shift : input.shift)';

    %  Intersect the expanded fifo(1) macrotime matrix with the fifo(2)
    %  macrotime matrix
    [~, index12, index2] = compareEvents(mT1, FIFO(2).macroT);

    %  Use values from the fifo(1) expanded matrix that have been found in
    %  the intersection with the fifo(2) matrix and intersect them with the
    %  fifo(3) macrotime matrix
    [~, index13, I3] = compareEvents(mT1(:, index12), FIFO(3).macroT);

    % Clear mT1 to save memory
    clear('mT1')

    %  Get the indices of the coincident events for fifo(1) and fifo(2)
    I1 = index12(index13);
    I2 = index2(index13);

    %% Create the image

    % Create an XY image from the data in fifo(1) and fifo(2)
    % Actually don't create the image, just get the indices of the bin
    % elements for the input arrays
    [~, ~, ~, binX, binY] = histcounts2(FIFO(2).microT(I2), ...
                                        FIFO(3).microT(I3), ...
                                        input.Xbins, ...
                                        input.Ybins);

    % Create a 3D array to hold the time decay data. Use 16-bit unsigned
    % integer
    XYZimage = zeros([numel(input.Ybins), ...
                      numel(input.Xbins), ...
                      numel(input.Tbins)], 'uint16');

    % Convert the microtime into an addressable bin coordinate
    binZ = bitshift(FIFO(1).microT(I1), -log2(input.Tstep)) + 1;

    % Make all values above the upper range of Tend zero (0)
    binZ(binZ > input.Tend) = 0;

    % Subtract the first time bin coordinate
    binZ = binZ - input.Tstart;

    % Get rid of photons outside the cropped area
    outPhots = binX == 0 | binY == 0 | binZ == 0;
    % Get rid of Z coordinate too close to the limits
    binX(outPhots) = [];
    binY(outPhots) = [];
    binZ(outPhots) = [];

    % Run through every photon and throw it into the right 3D histogram bin
    % Distinguish the situation when MCP shift is corrected and when it is
    % not
    if isempty(input.MCPshiftFile)
        % Run the photon assignment without MCP shift correction
        for i = 1 : numel(binX)
            Xcoord = binX(i);
            Ycoord = binY(i);
            Zcoord = binZ(i);
            XYZimage(Ycoord, Xcoord, Zcoord) = ...
                XYZimage(Ycoord, Xcoord, Zcoord) + 1;
        end
    else
        % Run the photon assignment with MCP shift correction
        % Simulate offset to the MCP shift
        binZoffset = (rand(size(binZ)) - 0.5) * input.MCPrandomWeight;
        % Number of Z bins
        Zlength = numel(input.Tbins);
        % Run the photon assignment without MCP shift correction
        for i = 1 : numel(binX)
            Xcoord = binX(i);
            Ycoord = binY(i);
            %
            Zcoord = max(min(binZ(i) - round(...
                MCPshift.(input.MCPshiftModel).offset(Ycoord, Xcoord) + ...
                binZoffset(i)), Zlength), 1);
            XYZimage(Ycoord, Xcoord, Zcoord) = ...
                XYZimage(Ycoord, Xcoord, Zcoord) + 1;
        end
    end
    % Save the frame
    if input.saveMAT
        save(input.outputStackName{frame, 1}, 'XYZimage', 'input')
    end
    % Save the frame as ICS
    if input.saveICS
        % Convert the stack to dip_image
        % Calculate the range of the histogram in seconds
        range = size(XYZimage, 3) * binWidth;
        % Save the ICS file
        exportICS2(XYZimage, ...
                   input.outputStackName{frame, 2}, ...
                   range, ...
                   binWidth * 1e12);
    end
    if nargout > 2
        intensityImage(:, :, frame) = sum(XYZimage, 3);
    end
end