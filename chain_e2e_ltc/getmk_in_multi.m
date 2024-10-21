function [mk_set] = getmk_in_multi(taskset, cn_st, sumC, WCRT, cnT)
    % GETMK_IN_MULTI calculate lanient range of (m, k) for each task 
    %               in a multi-chain model
    %   INPUTS: One taskset: taskset
    %           chain structure that shows task order in the chain: cn_st
    %           sum of execution time of the chain: sumC
    %           worst-case response time of the chain in multi-chain model: WCRT
    %           period of the chain: cnT

    %mk_set = cell(height(taskset),1);
    mk_set = repmat({NaN}, height(taskset), 1);
    if WCRT == -1
        for i = 1: length(cn_st)
            mk_set{taskset{:,"ID"} == cn_st(i)} =[-1, -1];  % if WCRT has not converged, no m,k information is available
        end
    else
        subset_cn = taskset(ismember(taskset.ID, cn_st), ["ID" "C" "T"]);

        if WCRT >= cnT
            for i = 1: length(cn_st)
                W = 2 * WCRT - sumC - subset_cn{(subset_cn{:,"ID"} == cn_st(i)), "C"};
                m = ceil(W / subset_cn{(subset_cn{:,"ID"} == cn_st(i)), "T"});
                mk_set{taskset{:,"ID"} == cn_st(i)} =[m, (m+1)];  
            end
        else
            for i = 1: length(cn_st)
                W = cnT + WCRT - sumC - subset_cn{(subset_cn{:,"ID"} == cn_st(i)), "C"};  
                m = ceil(W / subset_cn{(subset_cn{:,"ID"} == cn_st(i)), "T"});
                mk_set{taskset{:,"ID"} == cn_st(i)} =[m, (m+1)];  
            end
        end      
    end

        