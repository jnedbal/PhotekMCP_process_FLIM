function fifo = SPC830read(fname)
%SPC830read Import FIFO files from Becker&Hickl TCSPC cards
%   fifo = SPC830read(FNAME) returns a struct containing the setup
%   parameters, the macro times, micro times, gaps, macro time overflows
%   and gaps.
%
%   SPC830read supports FIFO data files from SPC-134, SPC-144, SPC-154 and
%   SPC-830 cards. It requires .spc and .set files in the same directory.
%   
%   Required input:
%       FNAME       : Link to the .spc filename.
%
%   Output:
%       fifo.set    : Setup header for the imported experiment.
%       fifo.MTC    : Macro time clock in seconds.
%       fifo.NRB    : Number of routing bits.
%       fifo.DINV   : Data invalid.
%       fifo.macroT : Macro time of photons.
%       fifo.microT : Micro time of photons.
%       fifo.rout   : Routing signals (inverted).
%       fifo.mtov   : Macro time overflow.
%       fifo.gap    : Possible recording gap due to 'FIFO Full'.
%   
%
%   See "The bh TCSPC Handbook"
%       (Page 485 in the Fourth Edition)

%   Copyright 2013 Jakub Nedbal
%   $Revision: 1.0.0.0 $  $Date: 2013/01/15 $


%% Import the setup file
% change the extention
setfile = horzcat(fname(1 : end - 4), '.set');

% List of information we want to get from the header:
hlist = {'Title', 'Version', 'Revision', 'Date', 'Time', 'Author', 'Company', 'Contents'};
vlist = cell(size(hlist));
rlist = true(size(hlist));

% Scan the set file for these occurences
fid = fopen(setfile);
tline = fgetl(fid);
while ischar(tline) && any(rlist)
    for i = find(rlist)
        if strfind(tline, hlist{i})
            vlist{i} = tline(strfind(tline, ' : ') + 3 : end);
            rlist(i) = false;
            break
        end
    end
    tline = fgetl(fid);
end

fifo.set = cell2struct(vlist', hlist');

% Read all the settings in the file
while ~feof(fid)
    tline = fgetl(fid);
    if regexp(tline, '  #')
        tline = textscan(tline, '%s', 'Delimiter', '#, []');
        tline = tline{:};
        % Don't save TR and WI paremeters
        if strcmp(tline{4}, 'TR') || strcmp(tline{4}, 'WI')
            continue
        end
        switch tline{7}
            case {'I', 'F'}
                tline{8} = str2double(tline{8});
            case 'B'
                tline{8} = logical(str2double(tline{8}));
        end
        fifo.set.(tline{4}).(tline{6}) = tline{8};
    end
end
fclose(fid);



%% Import the file content
fid = fopen(fname);

header = fread(fid, 4, '*uint8');

% Process the header.
% Bytes 0, 1, 2: Macro time clock in 0.1ns units ('500' for 50 ns macrotime
%                clock)
fifo.MTC = double(typecast([header(1 : 3); uint8(0)], 'uint32')) * 1e-10;

% Byte 4: bits 3 to 6 = number of routing bits
fifo.NRB = bitshift(bitand(header(4), 120), -3);

% Byte 4: bit 7 = 1 ('Data invalid')
fifo.DINV = logical(bitand(header(4), 128));

%Import the rest of the file
fifo.macroT = fread(fid, Inf, '*uint8');
fclose(fid);

%% Four bytes associate with each photon. Therefore we reshape the matrix
fifo.macroT = reshape(fifo.macroT, 4, numel(fifo.macroT) / 4);

%% The microtimes are in the third and fourth bytes
fifo.microT = fifo.macroT([3, 4], :);
fifo.macroT([3, 4], :) = [];

%% Get the ROUT as the 4 upper bites of the second byte of the macrotime
fifo.rout = bitshift(fifo.macroT(2, :), -4);

%% Make the upper four bits of the second byte of macroT zero
fifo.macroT(2, :) = bitand(fifo.macroT(2, :), 15);
fifo.macroT = double(typecast(fifo.macroT(:)', 'uint16'));

%% Invalid photons
inv = logical(bitand(fifo.microT(2, :), 128));

%% MTOV macro timer overflow
fifo.mtov = logical(bitand(fifo.microT(2, :), 64));

%% GAP "FIFO full error": Some earlier photons were presumably lost.
fifo.gap = logical(bitand(fifo.microT(2, :), 32));

%% Make the upper four bits of the second byte of microT zero
fifo.microT(2, :) = bitand(fifo.microT(2, :), 15);
fifo.microT = 4095 - typecast(fifo.microT(:)', 'uint16');

%% Use information about macrotime overflow to create a steadily increasing
%  macrotime. Each time 'inv' and 'mtov' are true, the value in macrotime
%  tells us how many times the macrotime timer has overflown.
%  We search for occurences where mtov and inv are true:
overflowMany = fifo.mtov & inv;

%  Where the macrotime overflows with a photon arrival, inv is false but
%  mtoc is true. We therefore search for these single overflows
overflowOne = fifo.mtov & ~inv;


%  We produce the cumulative sum of the number of overflows written in the
%  macrotime and the single overflows:
macroTadds = cumsum(overflowMany .* fifo.macroT + overflowOne) * 4096;

%  Now we add the macro time additions to the macro times:
fifo.macroT = fifo.macroT + macroTadds;

%% Dump the photons that are invalid
fifo.macroT(inv) = [];
fifo.microT(inv) = [];
fifo.rout(inv) = [];
fifo.mtov(inv) = [];
fifo.gap(inv) = [];

