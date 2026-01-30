function [num_drainage_timepoints, drainage_indices, contributing_rows, cluster_num_lakes, cluster_num_timepoints] = ...
    count_simultaneous_drainage(days_drainages)
% function for counting up clusters of drainage events: loops through each
% individual drainage event to find the days of that drainage (first function), 
% and then asks how many other drainage events' timings fall entirely within that 
% same drainage event's time window (second function)
% --> yields unique clusters of lake drainages

% outputs
num_drainage_timepoints = zeros(size(days_drainages, 1), 1);
drainage_indices = cell(size(days_drainages, 1), 1);
contributing_rows = cell(size(days_drainages, 1), 1);    
cluster_num_lakes = zeros(size(days_drainages, 1), 1);
cluster_num_timepoints = zeros(size(days_drainages, 1), 1);
    
% loop over each lake, looking for the drainage event for that lake
    for row = 1:size(days_drainages, 1)      
        % indices of that lake's drainage (==1)
        drainage_indices{row} = find(days_drainages(row, :) == 1);
        % number of 24-hr periods of that lake's drainage
        num_drainage_timepoints(row) = length(drainage_indices{row});
        % find the contributing rows for each drainage event
        contributing_rows{row} = find_contributing_rows(days_drainages, drainage_indices{row}, row);
        % update cluster summary variables
        cluster_num_lakes(row) = numel(contributing_rows{row});
        cluster_num_timepoints(row) = numel(drainage_indices{row});
    end
end

% see if any others lakes also drain within that same time window:
function contributing_rows = find_contributing_rows(matrix, drainage_indices, current_row)
    % Find rows with no 1's outside the event columns of the current row
    contributing_rows = [];
    for row = 1:size(matrix,1) % for each lake drainage
        if all(all(isnan(matrix(row, [1:drainage_indices(1)-1, drainage_indices(end)+1:end])))) && row ~= current_row
            contributing_rows = [contributing_rows, row];
            % add another row (lake) if drainage happens at same time
            % AND if that row is ~= the current lake
        end
    end
    if all(all(isnan(matrix(current_row, [1:drainage_indices(1)-1, drainage_indices(end)+1:end]))))
        contributing_rows = [contributing_rows, current_row];
        % also append the current_row (lake) to set of lakes
    end
end
