% This is a new version of the gap correcting software that will look for
% similar patterns of threes
% The FIFO data can be corrupt and if so, the code tries to jump over the
% gaps and reconstitute the rest of the data stream
if any(FIFO(1).gap) || any(FIFO(2).gap) || any(FIFO(3).gap)
    warndlg(['There are gaps in the data. ', ...
                     'It is a corrupt file'], ...
                    'SPC2image', 'modal');
    % Find a pattern in the time differentials following the gaps
    if sum(FIFO(1).gap) ~= sum(FIFO(2).gap) || ...
            sum(FIFO(1).gap) ~= sum(FIFO(3).gap)
        warndlg(['The gaps are not same on all three cards. ', ...
                         'Don''t know what to do!'], ...
                        'SPC2image', 'modal');
    else
        % Indices of the gaps in the streams
        %gi = [find(FIFO(1).gap); find(FIFO(2).gap); find(FIFO(3).gap)];
        gi{1} = find(FIFO(1).gap);
        gi{2} = find(FIFO(2).gap);
        gi{3} = find(FIFO(3).gap);
        % Maximum length of the stream 
        maxLen = ceil(max([numel(FIFO(2).macroT), numel(FIFO(3).macroT)] - 1) / 4);
        % Pattern LUT
        patLUT = struct;
        % There seem to the a problem, where the gap appears right at the 
        for i = 1 : numel(gi{1})
            % Look for the pattern in the delays between the photons in the
            % 1000 photons past the gap
            steps = 10 .^ (2 : 5);
            for j = 1 : numel(steps)
                % Find a pattern of four
                for k = 1 : 10
                    patLUT(i, j, k, 1).index = ...
                        gi{1}(i) + steps(j) + (k - 1) * 5 + (0 : 5);
                    % Create a pattern of three consecutive differences
                    patLUT(i, j, k, 1).pat = diff(FIFO(1).macroT(patLUT(i, j, k, 1).index));
                    % Create a permutation of the orders
                    pattern = [patLUT(i, j, k, 1).pat; ...
                               patLUT(i, j, k, 1).pat([2, 3, 4, 5, 1]); ...
                               patLUT(i, j, k, 1).pat([3, 4, 5, 1, 2]); ...
                               patLUT(i, j, k, 1).pat([4, 5, 1, 2, 3]); ...
                               patLUT(i, j, k, 1).pat([5, 1, 2, 3, 4])];
                    % Repeat the pattern so many times to fill the longest
                    % photons stream;
                    pattern = repmat(pattern, [1, maxLen]);
                    for m = [2, 3]
                        mDiff = diff(FIFO(m).macroT);
                        mDiff = abs(mDiff - pattern(:, 1 : numel(mDiff)));
                        mDiff = mDiff(:, 1 : end - 4) + ...
                                mDiff(:, 2 : end - 3) + ...
                                mDiff(:, 3 : end - 2) + ...
                                mDiff(:, 4 : end - 1) + ...
                                mDiff(:, 5 : end);
                        [R, C, V] = find(mDiff < 10);
                        
                        if ~isempty(R)
                            % fprintf('j: %d, k: %d m: %d; R: %d, C: %d\n', j, k, m, R, C)
                            patLUT(i, j, k, m).pat = diff(FIFO(m).macroT(C + (0 : 5)));
                            patLUT(i, j, k, m).index = C + (0 : 5);
                            
                        end
                    end
                end
                % Check if the search produced anything
                % Find the nonempty inputs in each of the search areas
                % Identify those that have nonempty one
                success = find(all(squeeze(any(cellfun(@numel, squeeze(struct2cell(patLUT(i, j, :, [2, 3])))), 1)), 2));
                for k = success'
                    % Find the index value
                    inVal = [FIFO(1).macroT(patLUT(i, j, k, 1).index(1)); ...
                             FIFO(2).macroT(patLUT(i, j, k, 2).index(1)); ...
                             FIFO(3).macroT(patLUT(i, j, k, 3).index(1))];
                    % Find the maximum index value
                    maxVal = max(inVal);
                    % Get the addition value
                    addVal = maxVal - inVal;
                    % Try several little offsets to find whats best value
                    % to add
                    %  Create a matrix of fifo(2) macrotimes that is up to a few (shift) 
                    %  values smaller and larger than the origininal one to look for
                    %  intersections with the other fifo macrotimes.
                    mT1 = repmat(FIFO(1).macroT(patLUT(i, j, k, 1).index(1) : end) + addVal(1), ...
                                 2 * input.shift + 1, 1) + ...
                                 (-input.shift : input.shift)';
                    % Add value indices
                    valAddIn = -input.shift : input.shift;
                    % Create an empty matrix of matches
                    matchCount = zeros(2, numel(valAddIn));
                    % flag for an OK match
                    OKmatch = false(2, 1);
                    for fifoIn = [2, 3]
                        % Try to find the best value to add
                        for valAdd = 1 : numel(valAddIn)
                            compMatX = FIFO(fifoIn).macroT(patLUT(i, j, k, fifoIn).index(1) : end) + addVal(fifoIn) + valAddIn(valAdd);
                            matchCount(fifoIn - 1, valAdd) = numel(compareEvents(mT1, compMatX));
                        end
                        % If the match count is more than half the
                        % remaining photons, assume the search was
                        % successful
                        OKmatch(fifoIn - 1) = (size(compMatX, 2) / max(matchCount(fifoIn - 1, :))) < 2;
                        % Find the best comparison
                        %[~, mIn] = max(matchCount);
                        %addVal(fifoIn) = addVal(fifoIn) + valAddIn(mIn);
                        % Correct the rest of the stream
                        %FIFO(fifoIn).macroT(sh(fifoIn) : end) = FIFO(fifoIn).macroT(sh(fifoIn) : end) + addVal(fifoIn);
                    end
                    % Check if both matches are OK
                    if all(OKmatch)
                        % Change all the indices of the macrotimes
                        % downstream of the gap
                        for fifoIn = 1 : 3
                            if fifoIn > 1
                                % Find the best comparison
                                [~, mIn] = max(matchCount(fifoIn - 1, :));
                                addVal(fifoIn) = addVal(fifoIn) + valAddIn(mIn);
                            end
                            FIFO(fifoIn).macroT(patLUT(i, j, k, fifoIn).index(1) : end) = FIFO(fifoIn).macroT(patLUT(i, j, k, fifoIn).index(1) : end) + addVal(fifoIn);
                        end
                        clear('compMatX', 'pattern', 'mDiff', 'mT1')
                        break
                    end
                end
                if all(OKmatch)
                    break
                end
            end
        end
    end
end