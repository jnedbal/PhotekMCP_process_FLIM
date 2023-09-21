function saveIRF(IRFmap, pixel, fname, binWidth)
% Function to save IRF into an ICS file for use in TRI2.
%   IRFmap      is a x, y, t TDC histogram from the SPAD camera
%   pixel       is the pixel coordinate to save
%   fname       is the filename of the ICS file
%   binWidth    is the width of each time bin in nanoseconds

% Keep only the needed data
IRFmap = IRFmap(pixel(1), pixel(2), :);
% Scale the IRFmap to the range of uint16
% if max(IRFmap, [], 'all') > intmax('uint16')
%     IRFmap = IRFmap / max(IRFmap, 'all') * intmax('uint16');
% end
% Create a matrix for saving the IRF
IRF = newim(size(IRFmap), 'uint16');
% The data needs to be restructured
IRF(:, :, :) = dip_image(IRFmap, 'uint16');
% Get the range of the histogram in seconds
range = 1e-9 * size(IRF, 3) * binWidth;
% IRF needs to have the correct order of the coordinates. This means time
% along the first axis, and 1x1 pixel along the two other axes.
IRF = permute(IRF, [3, 2, 1]);
% Save the IRF
exportICS3(IRF, fname, range, binWidth);