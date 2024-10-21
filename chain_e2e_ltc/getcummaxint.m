function cummaxint = getcummaxint(st, nd, cn_taskset, chn_struct, diff_cn_ins)
    % GETCUMMAXINT calculate cumulatively calculated maximum interval between two tasks 
    %   INPUTS: The index of calculation starting task in the chain: st
    %           The index of calculation ending task in the chain: nd
    %           A taskset with tasks consist the chain: cn_taskset
    %               column 1: ID task id to build chain structure
    %               column 2: C execution time which is upper bound of deadline miss 
    %               column 3: I maximum possible idle time of a task
    %           Corresponding chain structure: chn_struct
    %           If start and end tasks crosses over chain instances (and how many): diff_cn_ins 
    %               0: in the same instance, 
    %               1: start from one instance end in the next instance

    if (st > nd & diff_cn_ins == 0) | diff_cn_ins < 0
        ME = MException('MyComponent:paramenters set error', ...
            'Starting index %d is later than ending index %d in the same chain instance or diff_cn_ins %d set error', ...
            st, nd, diff_cn_ins);
        throw(ME)
    end

    %[7 9 6 3]

    if st == nd & diff_cn_ins == 0
        cummaxint = 0;
    else
        if nd == 1
            prec_cumint = getcummaxint(st, height(cn_taskset), cn_taskset, chn_struct, diff_cn_ins - 1); 
        else
            prec_cumint = getcummaxint(st, nd - 1, cn_taskset, chn_struct, diff_cn_ins); 
        end

        cummaxint = max(prec_cumint + cn_taskset{cn_taskset{:,"ID"} == chn_struct(nd), "C"}, cn_taskset{cn_taskset{:,"ID"} == chn_struct(nd), "I"});

    end



    
                 