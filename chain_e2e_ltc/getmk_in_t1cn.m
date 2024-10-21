function [mk_set] = getmk_in_t1cn(taskset, cn_st, sumC)
    % GETMK_IN_T1CN calculate lanient range of (m, k) for each task 
    %               in a chain with deadline miss termination as a single chain
    %   INPUTS: a taskset: taskset
    %           chain structure that shows task order in the chain: cn_st
    %           sum of execution time of the chain: sumC_cn

    mk_set = repmat({NaN}, height(taskset), 1);
    subset_cn = taskset(ismember(taskset.ID, cn_st), ["ID" "C" "I" "T" "ED"]);

    for i = 1: length(cn_st)
        nd = i - 1;
        df = 1;
        if nd == 0
            nd = length(cn_st);
            df = 0;
        end
        ii_cum_int = getcummaxint(i, nd, subset_cn, cn_st, df);
        ii_ltc = subset_cn{(subset_cn{:,"ID"} == cn_st(i)), "ED"} + sumC + ii_cum_int;
        
        % length of task "C" is added to cover the cases of deadline miss at the next job execution
        K = ceil(ii_ltc/subset_cn{(subset_cn{:,"ID"} == cn_st(i)), "T"});
        mk_set{taskset{:,"ID"} == cn_st(i)} =[(K-1), K];
    end

        