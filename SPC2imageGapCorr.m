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
        gi = [find(FIFO(1).gap); find(FIFO(2).gap); find(FIFO(3).gap)];
        % There seem to the a problem, where the gap appears right at the 
        for i = 1 : size(gi, 2)
            % Look for the pattern in the delays between the photons in the
            % 1000 photons past the gap
            % MacroTime gaps
            mtg1 = diff(FIFO(1).macroT(gi(1, i) + (1 : 1000)));
            mtg2 = diff(FIFO(2).macroT(gi(2, i) + (1 : 1000)));
            mtg3 = diff(FIFO(3).macroT(gi(3, i) + (1 : 1000)));
            % Calculate the shifts between the patterns
            sh12 = round(findshift(dip_image(mtg1), dip_image(mtg2), 'iter'));
            sh13 = round(findshift(dip_image(mtg1), dip_image(mtg3), 'iter'));
            % Combined shift position
            sh(1) = max([0, sh12, sh13]);
            sh(2) = sh(1) - sh12;
            sh(3) = sh(1) - sh13;
            % Start from the shifted positions
            mtg1 = diff(FIFO(1).macroT(gi(1, i) + (0 : 99) + sh(1)));
            mtg2 = diff(FIFO(2).macroT(gi(2, i) + (0 : 99) + sh(2)));
            mtg3 = diff(FIFO(3).macroT(gi(3, i) + (0 : 99) + sh(3)));
            % Calculate the shifts between the patterns
            sh12A = round(findshift(dip_image(mtg1), dip_image(mtg2), 'iter'));
            sh13A = round(findshift(dip_image(mtg1), dip_image(mtg3), 'iter'));
            % Combined shift position
            sh1A = max([0, sh12A, sh13A]);
            sh2A = sh1A - sh12A;
            sh3A = sh1A - sh13A;
            % Store the shifts
            sh(1) = sh(1) + sh1A;
            sh(2) = sh(2) + sh2A;
            sh(3) = sh(3) + sh3A;
            % Start from the shifted positions
            mtg1 = diff(FIFO(1).macroT(gi(1, i) + (0 : 99) + sh(1)));
            mtg2 = diff(FIFO(2).macroT(gi(2, i) + (0 : 99) + sh(2)));
            mtg3 = diff(FIFO(3).macroT(gi(3, i) + (0 : 99) + sh(3)));
            % Hopefully, by now, there is a reasonable match between the
            % three channel macrotime differentials
            % With this, try to brute crack where exactly the three streams
            % meet by finding two consecutive differentials
            for j = 9 : (numel(mtg1) - 10)
                % try to find the exact offset to the point by searching
                % the 17 values around the central position. Keep going
                % until you find the match for all three channels.
                j2 = find(sum(abs(mtg1([j, j + 1])' - [mtg2(j - 8 : j + 8); mtg2((j - 8 : j + 8) + 1)])) < 2 * input.shift);
                j3 = find(sum(abs(mtg1([j, j + 1])' - [mtg3(j - 8 : j + 8); mtg3((j - 8 : j + 8) + 1)])) < 2 * input.shift);
                % Check if a match has been found
                if ~isempty(j2) && ~isempty(j3)
                    % Add the successfully found indices
                    sh(1) = gi(1, i) + sh(1) + j;
                    sh(2) = gi(2, i) + sh(2) + j2;
                    sh(3) = gi(3, i) + sh(3) + j3;
                    % Find the index value
                    inVal = [FIFO(1).macroT(sh(1)); ...
                             FIFO(2).macroT(sh(2)); ...
                             FIFO(3).macroT(sh(3))];
                    % Find the maximum index value
                    maxVal = max(inVal);
                    % Get the addition value
                    addVal = maxVal - inVal;
                    % Change all the indices of the macrotime 1 downstream
                    % of the gap
                    FIFO(1).macroT(sh(1) : end) = FIFO(1).macroT(sh(1) : end) + addVal(1);
                    % Try several little offsets to find whats best value
                    % to add
                    %  Create a matrix of fifo(2) macrotimes that is up to a few (shift) 
                    %  values smaller and larger than the origininal one to look for
                    %  intersections with the other fifo macrotimes.
                    mT1 = repmat(FIFO(1).macroT(sh(1) : end), ....
                                 2 * input.shift + 1, 1) + ...
                            (-input.shift : input.shift)';
                    % Add value indices
                    valAddIn = -input.shift : input.shift;
                    % Create an empty matrix of matches
                    matchCount = zeros(1, numel(valAddIn));
                    for fifoIn = [2, 3]
                        % Try to find the best value to add
                        for valAdd = 1 : numel(valAddIn)
                            compMatX = FIFO(fifoIn).macroT(sh(fifoIn) : end) + addVal(fifoIn) + valAddIn(valAdd);
                            matchCount(valAdd) = numel(compareEvents(mT1, compMatX));
                        end
                        % Find the best comparison
                        [~, mIn] = max(matchCount);
                        addVal(fifoIn) = addVal(fifoIn) + valAddIn(mIn);
                        % Correct the rest of the stream
                        FIFO(fifoIn).macroT(sh(fifoIn) : end) = FIFO(fifoIn).macroT(sh(fifoIn) : end) + addVal(fifoIn);
                    end
                    clear('compMatX')
                    break
                end
            end
        end
    end
end