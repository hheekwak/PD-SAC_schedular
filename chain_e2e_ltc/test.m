% Example cell array of different-length arrays
i = 1;
i_chains = st_ts(i, 1:ts_ch_size(i));

% Step 1: Sort tasks in individual arrays within the cell array
t_sorted_chains = cellfun(@sort, i_chains, 'UniformOutput', false);

% Step 2: Save the original indices
org_idx = 1:length(i_chains);

% Step 3: Sort by values first and then by length
% Pad arrays with inf to make them the same length for comparison
[~, sorted_idx] = sortrows(cell2mat(cellfun(@(x) [x, inf(1, max(cellfun(@length, t_sorted_chains))-length(x))], t_sorted_chains, 'UniformOutput', false)'));

% Step 4: Sort the original cell array using the sorted indices
pr_chain = i_chains(sorted_idx);

% Step 5: Store the original indices in the order of the sorted cell
sorted_original_indices = org_idx(sorted_idx);

% Display the sorted cell array and the original indices
disp('Sorted cell array with sorted items:');
disp(sorted_cell);
disp('Original indices in sorted order:');
disp(sorted_original_indices);