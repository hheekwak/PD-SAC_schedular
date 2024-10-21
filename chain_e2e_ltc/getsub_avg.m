function [out_avg] = getsub_avg(row_col, mtx)
    % GETSUB_AVG Get average of each specific-sized submatrix in a matrix
    %   INPUTS: row_col : size of submatrix that is average window
    %           mtx : original matrix

    % block-sized filter
    filter = ones(row_col(1), row_col(2)); 
    
    % Scale
    n = row_col(1) * row_col(2);

    % result array 
    out_avg = zeros(size(mtx)./ size(filter));

    % averaging
    for i = 1: row_col(1) : size(mtx, 1) - row_col(1) + 1 % row
        for j = 1: row_col(2): size(mtx, 2) - row_col(2) + 1 % column
            % extract the region 
            rgn = mtx(i:i + row_col(1) - 1, j: j+ row_col(2) - 1);
            % summation
            summation = sum(sum(rgn .* filter));
            % average and fill the result array
            out_avg((i+row_col(1)-1)/row_col(1), (j+row_col(2)-1)/row_col(2)) = summation / n;  
        end
    end
end
