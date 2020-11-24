function [value, index1, index2] = compareEvents(expandedMatrix, compareTo)

% Find the intersection between the two matrices
[value, index1, index2] = intersect(expandedMatrix, compareTo);

% The first matrix has been expanded to allow for some deviation of the
% values. Compress the index to only represent the original value.
index1 = ceil(index1 / size(expandedMatrix, 1));

% Convert the data to integers to save memory
index1 = uint32(index1);
index2 = uint32(index2);
