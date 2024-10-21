function [mk_set] = getmk_in_t2cn(taskset, cn_st, latency)
    % GETMK_IN_T2CN calculate lanient range of (m, k) for each task 
    %               in a chain without deadline miss termination as a single chain
    %   INPUTS: a taskset: taskset
    %           chain structure that shows task order in the chain: cn_st
    %           The longest latency from the execution starting points of τi in chain's jth instance
    %           to the completion point of τi−1 in chain's (j+1)th instance: latency

    %mk_set = cell(height(taskset),1);
    mk_set = repmat({NaN}, height(taskset), 1);
    subset_cn = taskset(ismember(taskset.ID, cn_st), ["ID" "C" "I" "T" "ED"]);

    for i = 1: length(cn_st)
        ii_ltc = subset_cn{(subset_cn{:,"ID"} == cn_st(i)), "ED"} + latency;
        K = ceil(ii_ltc/subset_cn{(subset_cn{:,"ID"} == cn_st(i)), "T"});
        mk_set{(taskset{:,"ID"} == cn_st(i))} =[(K-1), K];
    end