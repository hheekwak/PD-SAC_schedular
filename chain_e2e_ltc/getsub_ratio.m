function [p_ratio, n_ratio] = getsub_ratio(row_col, mtx)
    % GETSUB_RATIO Get schdulability ratio of each specific-sized row for each column in a matrix
    %   INPUTS: row_col : size of count window and number of chains
    %           mtx : original matrix

    p_ratio = NaN(length(mtx)/row_col(1), row_col(2));       % schedulability ratio in each case
    n_ratio = NaN(length(mtx)/row_col(1), row_col(2));       % non-convergence rate in each case
    for i = 1: length(mtx)/row_col(1)
        s = (i - 1) * row_col(1) + 1;
        e = i * row_col(1);
        
        for j = 1: row_col(2)
            p_ratio(i, j) = sum(mtx(s:e, j) == 1)/row_col(1);
            n_ratio(i, j) = sum(mtx(s:e, j) == -1)/row_col(1);
        end
    end

