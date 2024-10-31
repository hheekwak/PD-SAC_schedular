%% Experiments
% Parameters
clc; close all; clear;

save_mk = 0;    % save (m,k) of all tasks as a file

ts_size = 15; % size of a taskset (number of tasks in a taskset)
n_ts = 1000;  % number of tasksets for each case (e.g. utilization, max T position...)

% manually set chain configuration (chain structure)
cn_st1 = [1 2 3 4 5]; 
cn_st2 = [6 7 8 9 10];
cn_st3 = [11 12 13 14 15];
%cn_st4 = [10 11 12];
%cn_st5 = [13 14 15];

st_ts = {cn_st1 cn_st2 cn_st3}; % cn_st4 cn_st5}; % structure of taskset e.g.{cn_st1 cn_st2}

%% data read
% also consider re-nameing names of result and figure file 

% read data
data1 = load('../data/chain_case17_1000.txt');
data1 = array2table(data1, 'VariableNames', {'T' 'C' 'D' 'ID' 'Prior'});

% add ED, upperbound of potential delay from its release point to execution
% add empty column for I, upperbound of idle time of each task

data1.ED = data1.T - data1.C;
data1.I(:) = NaN;

%% variables

% split tasksets and set result arrays
taskset = cell(height(data1)/ts_size, 1);
sumC_cn = zeros(length(taskset), length(st_ts));

                              
e2e_ltc_cn = NaN(length(taskset), length(st_ts));      % Type 1: end-to-end latency in a single chain, including ED,DM 
e2e_ltc_cn_wo_ED = NaN(length(taskset), length(st_ts)); % Type 1: end-to-end latency in a single chain, including DM, not ED 
maxDM_cn = zeros(length(taskset), length(st_ts));      % Type 1: max length of deadline miss for type1 WCRT 
T_cn_t1 = zeros(length(taskset), length(st_ts));       % Type 1: set T of each chain with its max latency and T_1 with Deadline miss termination 
WCRT_chain_1 = NaN(length(taskset), length(st_ts));    % Type 1: WCRT with deadline miss termination 
schd_able_1 = NaN(length(taskset), length(st_ts));     % Type 1: schedulability of each chain 

e2e_ltc_wo_DM = NaN(length(taskset), length(st_ts));      % Type 2: end-to-end latency in a single chain, including ED, not DM 
e2e_ltc_wo_DM_ED = NaN(length(taskset), length(st_ts));   % Type 2: end-to-end latency in a single chain, no ED, no DM considered 
T_cn_t2 = zeros(length(taskset), length(st_ts));          % Type 2: set T of each chain with its max latency and T_1 without Deadline miss termination 
WCRT_chain_2 = NaN(length(taskset), length(st_ts));       % Type 2: WCRT without deadline miss termination 
schd_able_2 = NaN(length(taskset), length(st_ts));        % Type 2: schedulability of each chain 

utilization = NaN(length(taskset), 4);
mean_util = NaN(length(taskset)/n_ts, 4);
median_util = NaN(length(taskset)/n_ts, 4);
count_unsch_ts = NaN(length(taskset)/n_ts, 4);


%maxT_cn = zeros(length(taskset), length(st_ts));      % get max length T of each chain 
%maxT_pos = -1 * ones(length(taskset), length(st_ts)); % to save the position of max T in chain
hyperperiods = zeros(length(taskset), length(st_ts));  % to save the hyperperiods of each chain
T_cn_h = zeros(length(taskset), length(st_ts));          % set T of each chain with its max latency and hyperperiod
WCRT_chain_h1 = NaN(length(taskset), length(st_ts));     % Type1 and chainT = hyperperiod
schd_able_h1 = NaN(length(taskset), length(st_ts));       % schedulability of each chain
WCRT_chain_h2 = NaN(length(taskset), length(st_ts));     % Type2 and chainT = hyperperiod
schd_able_h2 = NaN(length(taskset), length(st_ts));       % schedulability of each chain
%TF = false(length(taskset), length(st_ts));            % check if T of chain is longer than any T in chain

%% calculate
for i = 1: height(taskset)
    taskset(i) = {data1((i-1)*ts_size+1:i*ts_size,["ID" "T" "C" "ED" "I"])};
    
    % set up arrays to save cumulative interval for each chain
    cum_interval = NaN(1, length(st_ts));     % Type 1: arr for cumulative interval with deadline miss
    
    cum_int_wo_dm = NaN(1, length(st_ts));    % Type 2: arr for cumulative interval without deadline miss for latency (exclude idle time of task1)
    max_idle_in_chn = NaN(1, length(st_ts));  % Type 2: arr for longest max idle in a chain for (m,k) (include idle time of task1)
    
    % calculate sum of execution time for each chain & max idle time for each task
    for j = 1: length(st_ts)
        for k = 1: length(st_ts{j})     % st_ts{j} jth chain configuration e.g.[3 2 5 6 7 8 1]
            sumC_cn(i, j) = sumC_cn(i, j) + taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "C"};
            maxDM_cn(i, j) = max(maxDM_cn(i, j), taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "C"});
            %if maxT_cn (i,j) < taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"} 
            %    maxT_cn (i,j) = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"};
            %    maxT_pos(i,j) = k;
            %end

            %maxT_cn(i, j) = max(maxT_cn(i, j), taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"});
            
            % Method 2: hyper periods for each chain
            if k == 1
                 hyperperiods(i, j) = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"};
            end
            hyperperiods(i, j) = lcm(hyperperiods(i, j), taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"});
        end    

        % max idle time for each task
        for k = 1: length(st_ts{j}) 
            if taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"} > sumC_cn(i, j)
                taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "I"} = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"} - sumC_cn(i, j);
            else
                taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "I"} = 0;
            end
        end

        % Cumulatively calculated maximal interval of each chain
        % range of cumulative interval, from nth and to nth task in a chain
        % the order in a chain = index of chain structure
        st_order = 1;
        nd_order = length(st_ts{j});
        subset_cn = taskset{i}(ismember(taskset{i}.ID, st_ts{j}), ["ID" "C" "I" "ED"]);
        
        cum_interval(j) = getcummaxint(st_order, nd_order, subset_cn, st_ts{j}, 0);     % Type 1: cumulative interval with deadline miss
        cum_int_wo_dm(j) = max(taskset{i}{ismember(taskset{i}.ID, st_ts{j}(2:end)), "I"});  % Type 2: idle time of task1 doesn't affect to end-to-end latency
        max_idle_in_chn(j) = max(taskset{i}{ismember(taskset{i}.ID, st_ts{j}), "I"});       % Type 2: longest max idle time in a chain for m,k calculation

        % Latency of a chain 
        %       e2e_ltc_cn:       Type 1 end-to-end latency as a single chain (from the release time of start-task)
        %                               used for setting chian T for WCRT, considering execution delay and deadline miss
        %       e2e_ltc_cn_wo_ED: Type 1 latency between execution initiation of start-task and competion of end-task 
        %                               used for setting chain C for WCRT, deadline miss termination considered, excluded execution delay
        %       e2e_ltc_wo_DM:    Type 2 end-to-end latency as a single chain 
        %                               used for setting chain T for WCRT, execution delay included, deadline miss termination not considered
        %       e2e_ltc_wo_ED:    Type 2 latency between execution initiation of start-task and compeletion of end-task
        %                               used for setting chain C for WCRT (execution delay, deadline miss) not considered 
        
        e2e_ltc_cn(i, j) = subset_cn{(subset_cn{:,"ID"} == st_ts{j}(st_order)), "ED"} + sumC_cn(i, j) + cum_interval(j);
        e2e_ltc_cn_wo_ED(i, j) = sumC_cn(i, j) + cum_interval(j);

        e2e_ltc_wo_DM(i, j) = subset_cn{(subset_cn{:,"ID"} == st_ts{j}(st_order)), "ED"} + sumC_cn(i, j) + cum_int_wo_dm(j);   
        e2e_ltc_wo_DM_ED(i, j) = sumC_cn(i, j) + cum_int_wo_dm(j);         
        
        % Find lenient (m,k) of each task within each chain in single-chain model
        % Type 1: deadline miss termination considered
        taskset{i}.("mk_s_ty1_"+string(j)) = getmk_in_t1cn(taskset{i}, st_ts{j}, sumC_cn(i,j));
        % Type 2: deadline miss termination NOT considered
        taskset{i}.("mk_s_ty2_"+string(j)) = getmk_in_t2cn(taskset{i}, st_ts{j}, (sumC_cn(i,j) + max_idle_in_chn(j)));
    

        % Set T(period) of each chain from its max e2e latency with and without DM and T of task1
        cn_T1 = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(1)), "T"};           % T_1 of each chain 
        
        T_cn_t1(i, j) = ceil(e2e_ltc_cn(i, j)/cn_T1) * cn_T1;       % Type 1: with deadline miss termination 
        T_cn_t2(i, j) = ceil(e2e_ltc_wo_DM(i, j)/cn_T1) * cn_T1;    % Type 2: without deadline miss termination
        
        % Method2 Set T(period) of each chain from tasks' hyper period
        T_cn_h(i, j) = ceil(e2e_ltc_cn(i, j)/hyperperiods(i, j)) * hyperperiods(i, j);
       
    end

    % WCRT with T of chains and C (latency or sumC) of chains
    % Type1 (with deadline miss termination) chainT from max e2e latency and T of task1
    R_1 = WCRT_t1chain(maxDM_cn(i,:), e2e_ltc_cn_wo_ED(i,:), T_cn_t1(i,:));
    WCRT_chain_1(i,:) = R_1(:,1)';
    schd_able_1(i,:) = R_1(:,2)';

    % Type2 (No deadline miss termination) chainT from max e2e latency and T of task1
    R_2 = WCRT_t2chain(sumC_cn(i,:), e2e_ltc_wo_DM_ED(i,:), T_cn_t2(i,:));
    WCRT_chain_2(i,:) = R_2(:,1)';
    schd_able_2(i,:) = R_2(:,2)';

    % Method2 hyper period Type 1
    R_h1 = WCRT_c_h1(maxDM_cn(i,:), sumC_cn(i,:), e2e_ltc_cn_wo_ED(i,:), T_cn_h(i,:));
    WCRT_chain_h1(i,:) = R_h1(:,1)';
    schd_able_h1(i,:) = R_h1(:,2)';

    % Method2 hyper period Type 2
    R_h2 = WCRT_c_h2(sumC_cn(i,:), e2e_ltc_wo_DM_ED(i,:), T_cn_h(i,:));
    WCRT_chain_h2(i,:) = R_h2(:,1)';
    schd_able_h2(i,:) = R_h2(:,2)';

    % Find lenient (m,k) of each task within each chain in multi-chain model
    for j = 1: length(st_ts) 
        % Find max(m,k) of each task with WCRT in a multi-chain model
        % Type 1: deadline miss termination considered
        taskset{i}.("mk_m_ty1_"+string(j)) = getmk_in_multi(taskset{i}, st_ts{j}, sumC_cn(i,j), WCRT_chain_1(i,j), T_cn_t1(i,j));
        % Type 2: deadline miss termination NOT considered
        taskset{i}.("mk_m_ty2_"+string(j)) = getmk_in_multi(taskset{i}, st_ts{j}, sumC_cn(i,j), WCRT_chain_2(i,j), T_cn_t2(i,j));
        % Method2 Type 1: deadline miss termination considered
        taskset{i}.("mk_m_ty1_h"+string(j)) = getmk_in_multi(taskset{i}, st_ts{j}, sumC_cn(i,j), WCRT_chain_h1(i,j), T_cn_h(i,j));
        % Method2 Type 2: deadline miss termination NOT considered
        taskset{i}.("mk_m_ty2_h"+string(j)) = getmk_in_multi(taskset{i}, st_ts{j}, sumC_cn(i,j), WCRT_chain_h2(i,j), T_cn_h(i,j));

    end

    % Find the minimum utilization of taskset
    sumUtil_1 = 0;
    sumUtil_2 = 0;
    sumUtil_h1 = 0;
    sumUtil_h2 = 0;
    % Type 1 with lanient m,k applied
    for j = 1: length(st_ts)
        for k = 1: length(st_ts{j})
            T = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"};
            C = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "C"};
            K_1 = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "mk_m_ty1_"+string(j)}{1}(2);
            % min utilization is not available if it is not schedulable
            if K_1 == -1
                sumUtil_1 = -1;
                break;
            else
                sumUtil_1 = sumUtil_1 + ((C/T) * 1/K_1);
            end
        end
    end
    
    % Type 2 with lanient m,k applied
    for j = 1: length(st_ts)
        for k = 1: length(st_ts{j})
            T = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"};
            C = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "C"};
            K_2 = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "mk_m_ty2_"+string(j)}{1}(2);
            % min utilization is not available if it is not schedulable
            if K_2 == -1
                sumUtil_2 = -1;
                break;
            else
                sumUtil_2 = sumUtil_2 + ((C/T) * 1/K_2);
            end
        end
    end
    % Type 1 hyperperiod with lanient m,k applied
    for j = 1: length(st_ts)
        for k = 1: length(st_ts{j})
            T = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"};
            C = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "C"};
            K_h1 = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "mk_m_ty1_h"+string(j)}{1}(2);
            % min utilization is not available if it is not schedulable
            if K_h1 == -1
                sumUtil_h1 = -1;
                break;
            else
                sumUtil_h1 = sumUtil_h1 + ((C/T) * 1/K_h1);
            end
        end
    end
    % Type 2 hyperperiod with lanient m,k applied
    for j = 1: length(st_ts)
        for k = 1: length(st_ts{j})
            T = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "T"};
            C = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "C"};
            K_h2 = taskset{i}{(taskset{i}{:,"ID"} == st_ts{j}(k)), "mk_m_ty2_h"+string(j)}{1}(2);
            % min utilization is not available if it is not schedulable
            if K_h2 == -1
                sumUtil_h2 = -1;
                break;
            else
                sumUtil_h2 = sumUtil_h2 + ((C/T) * 1/K_h2);
            end
        end
    end

    utilization(i, 1) = sumUtil_1;
    utilization(i, 2) = sumUtil_2;
    utilization(i, 3) = sumUtil_h1;
    utilization(i, 4) = sumUtil_h2;
    
end

%% averages and ratios by utilization

if n_ts ~= 1
    avg_unit = [n_ts 1];

    for u = 1: (length(taskset)/n_ts)
        row_s = (u-1) * n_ts + 1;
        row_e = u * n_ts;
        for m = 1: 4
            % sub-section for each utilization and each type or method
            sub_util = utilization(row_s:row_e, m);
            % Count the number of -1s in the array
            count_unsch_ts(u, m) = sum(sub_util == -1);
            % Filter out unschedulable taskset utilization value (-1)
            filtered_util = sub_util(sub_util ~= -1);
            % Calculate the average of schedulable taskset's utilization
            mean_util(u, m) = mean(filtered_util);
            % Find the median value of schedulable taskset's utilization
            median_util(u, m) = median(filtered_util);
        end
    end 

    % out_avg_t1 = getsub_avg(avg_unit, e2e_ltc_cn); % latency with deadline miss termination
    % out_avg_t2 = getsub_avg(avg_unit, e2e_ltc_wo_DM); % latency without deadline miss termination
    %exe_avg = getsub_avg(avg_unit, sumC_cn);
    %prd_avg = getsub_avg(avg_unit, avgT_cn);
    
    
    rt_unit = [n_ts length(st_ts)];
    % p_ratio: schedulable cases, n_ratio: not converged cases 
    [p_ratio_1, n_ratio_1] = getsub_ratio(rt_unit, schd_able_1); % Type 1: with DM termination
    d_ratio_1 = 1 - p_ratio_1 - n_ratio_1;                       % Type 1: WCRT > T_cn but converged
    
    [p_ratio_2, n_ratio_2] = getsub_ratio(rt_unit, schd_able_2); % Type 2: without DM termination 
    d_ratio_2 = 1 - p_ratio_2 - n_ratio_2;                       % Type 2: WCRT > T_cn but converged
    
    [p_ratio_h1, n_ratio_h1] = getsub_ratio(rt_unit, schd_able_h1);
    d_ratio_h1 = 1 - p_ratio_h1 - n_ratio_h1;     % WCRT > T_cn but converged
    
    [p_ratio_h2, n_ratio_h2] = getsub_ratio(rt_unit, schd_able_h2);
    d_ratio_h2 = 1 - p_ratio_h2 - n_ratio_h2;     % WCRT > T_cn but converged
end

%% save result as txt files

% Type 1: a multi-chain model with deadline miss termination

fileID_1_m = fopen('../output/cn_e2eltc_case17_1000_type1_multi.txt', 'w');

fprintf(fileID_1_m,'count of unschedulable, Mean and Median of min. util (m,k considered) per utilization \n');
fprintf(fileID_1_m,'%8s %8s %8s \r\n', 'count', 'Mean', 'Median');
for util = 1: length(count_unsch_ts) 
        fprintf(fileID_1_m,'%8d %8.4f %8.4f \r\n', count_unsch_ts(util, 1), mean_util(util, 1), median_util(util, 1));
end

fprintf(fileID_1_m,'\n Ratio of chain completed per utilization \n');
fprintf(fileID_1_m,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(p_ratio_1)
    for j = 1:length(st_ts)
        fprintf(fileID_1_m,'%8.4f', p_ratio_1(i,j));
    end
    fprintf(fileID_1_m,'\r\n');
end

fprintf(fileID_1_m,'\n Ratio of chain delayed per utilization \n');
fprintf(fileID_1_m,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(d_ratio_1)
    for j = 1:length(st_ts)
        fprintf(fileID_1_m,'%8.4f', d_ratio_1(i,j));
    end
    fprintf(fileID_1_m,'\r\n');
end

fprintf(fileID_1_m,'\n Ratio of chain not-converged per utilization \n');
fprintf(fileID_1_m,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(n_ratio_1)
    for j = 1:length(st_ts)
        fprintf(fileID_1_m,'%8.4f', n_ratio_1(i,j));
    end
    fprintf(fileID_1_m,'\r\n');
end

fprintf(fileID_1_m,'\n Indiv. Chain WCRT in multi-chain model with deadline miss termination \n ');
fprintf(fileID_1_m,'%8s %8s %8s \r\n','chain 1', 'chain 2', 'chain 3');
for i = 1:length(taskset)
    for j = 1:length(st_ts)
        fprintf(fileID_1_m,'%8.1f', WCRT_chain_1(i,j));
    end
    fprintf(fileID_1_m,'\r\n');
end

% (m,k) for each task (multi, type 1)
if save_mk == 1
    fprintf(fileID_1_m,'\n Indiv task (m,k) in a multi-chain model with deadline miss termination \n %6s %4s %4s %4s %8s \r\n','#set', '#chn', '#tsk', '#tid', '(m,k)');
    for i = 1:length(taskset)
        for j = 1:length(st_ts)
            mk_data = taskset{i}.("mk_m_ty1_"+string(j));
            for k = 1: length(st_ts{j})
                fprintf(fileID_1_m,'%6d %4d %4d %4d %4d %4d \r\n', i, j, k, st_ts{j}(k), mk_data{(taskset{i}{:,"ID"} == st_ts{j}(k))});
            end
        end
    end
end

fclose(fileID_1_m);

% Type 2: a multi-chain model without deadline miss termination 
fileID_2_m = fopen('../output/cn_e2eltc_case17_1000_type2_multi.txt', 'w');

fprintf(fileID_2_m,'count of unschedulable, Mean and Median of min. util (m,k considered) per utilization \n');
fprintf(fileID_2_m,'%8s %8s %8s \r\n', 'count', 'Mean', 'Median');
for util = 1: length(count_unsch_ts) 
        fprintf(fileID_2_m,'%8d %8.4f %8.4f \r\n', count_unsch_ts(util, 2), mean_util(util, 2), median_util(util, 2));
end

fprintf(fileID_2_m,'\n Ratio of chain completed per utilization \n');
fprintf(fileID_2_m,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(p_ratio_2)
    for j = 1: length(st_ts)
    fprintf(fileID_2_m,'%8.4f', p_ratio_2(i,j));
    end
    fprintf(fileID_2_m,'\r\n');
end

fprintf(fileID_2_m,'\n Ratio of chain delayed per utilization \n');
fprintf(fileID_2_m,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(d_ratio_2)
    for j = 1: length(st_ts)
        fprintf(fileID_2_m,'%8.4f',d_ratio_2(i,j));
    end
    fprintf(fileID_2_m,'\r\n');
end

fprintf(fileID_2_m,'\n Ratio of chain not-converged per utilization \n');
fprintf(fileID_2_m,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(n_ratio_2)
    for j = 1: length(st_ts)
        fprintf(fileID_2_m,'%8.4f',n_ratio_2(i,j));
    end
    fprintf(fileID_2_m,'\r\n');
end

fprintf(fileID_2_m,'\n Indiv. Chain WCRT in multi-chain model without deadline miss termination \n ');
fprintf(fileID_2_m,'%8s %8s %8s \r\n','chain 1', 'chain 2', 'chain 3');
for i = 1:length(taskset)
    for j = 1: length(st_ts)
        fprintf(fileID_2_m,'%8.1f',WCRT_chain_2(i,j));
    end
    fprintf(fileID_2_m,'\r\n');
end

% (m,k) for each task (multi, type 2)
if save_mk == 1
    fprintf(fileID_2_m, '\n Indiv task (m,k) in a multi-chain model without deadline miss termination \n %6s %4s %4s %4s %8s \r\n','#set', '#chn', '#tsk', '#tid', '(m,k)');
    for i = 1:length(taskset)
        for j = 1:length(st_ts)
            mk_data = taskset{i}.("mk_m_ty2_"+string(j));
            for k = 1: length(st_ts{j})
                fprintf(fileID_2_m,'%6d %4d %4d %4d %4d %4d \r\n', i, j, k, st_ts{j}(k), mk_data{(taskset{i}{:,"ID"} == st_ts{j}(k))});
            end
        end
    end
end

fclose(fileID_2_m);

% a single-chain model with deadline miss termination (single, type 1)
fileID_1_s = fopen('../output/cn_e2eltc_case17_1000_type1_single.txt', 'w');

% Average latency per utilization (single, type 1)
% fprintf(fileID_1_s,'\n Avg. Latency per case with deadline miss termination \n');
% fprintf(fileID_1_s,'%8s %8s %8s %8s %8s \r\n','avg cn1', 'avg cn2', 'avg cn3', 'avg cn4', 'avg cn5');
% for i = 1:length(out_avg_t1)
%     for j = 1: length(st_ts)
%         fprintf(fileID_1_s,'%8.1f', out_avg_t1(i,j));
%     end
%     fprintf(fileID_1_s,'\r\n');
% end

% Indivicual latency (single, type 1)
fprintf(fileID_1_s,'\n Indiv. Chain Latency as single chain with deadline miss termination \n');
fprintf(fileID_1_s,'%8s %8s %8s %8s %8s \r\n','chain 1', 'chain 2', 'chain 3', 'chain 4', 'chain 5');
for i = 1:length(taskset)
    for j = 1: length(st_ts)
        fprintf(fileID_1_s,'%8.1f', e2e_ltc_cn(i,j));
    end
    fprintf(fileID_1_s,'\r\n');
end

% (m,k) for each task (single, type 1)
if save_mk == 1
    fprintf(fileID_1_s,'\n Indiv task (m,k) in each chain with deadline miss termination \n %6s %4s %4s %4s %8s \r\n','#set', '#chn', '#tsk', '#tid', '(m,k)');
    for i = 1:length(taskset)
        for j = 1:length(st_ts)
            mk_data = taskset{i}.("mk_s_ty1_"+string(j));
            for k = 1: length(st_ts{j})
                fprintf(fileID_1_s,'%6d %4d %4d %4d %4d %4d \r\n', i, j, k, st_ts{j}(k), mk_data{(taskset{i}{:,"ID"} == st_ts{j}(k))});
            end
        end
    end
end

fclose(fileID_1_s);

% a single-chain model without deadline miss termination (single, type 2)
fileID_2_s = fopen('../output/cn_e2eltc_case17_1000_type2_single.txt', 'w');

% Average latency per utilization (single, type 2)
% fprintf(fileID_2_s,'\n Avg. Latency per case with deadline miss termination \n');
% fprintf(fileID_2_s,'%8s %8s %8s %8s %8s \r\n','avg cn1', 'avg cn2', 'avg cn3', 'avg cn4', 'avg cn5');
% for i = 1:length(out_avg_t2)
%     for j = 1: length(st_ts)
%         fprintf(fileID_2_s,'%8.1f',out_avg_t2(i,j));
%     end
%     fprintf(fileID_2_s,'\r\n');
% end

% Individual latency (single, type 2)
fprintf(fileID_2_s,'\n Indiv. Chain Latency as single chain without deadline miss termination \n %8s %8s %8s %8s %8s \r\n','chain 1', 'chain 2', 'chain 3', 'chain 4', 'chain 5');
for i = 1:length(taskset)
    for j = 1: length(st_ts)
        fprintf(fileID_2_s,'%8.1f', e2e_ltc_wo_DM(i,j));
    end
    fprintf(fileID_2_s,'\r\n');
end

% (m,k) for each task (single, type 2)
if save_mk == 1
    fprintf(fileID_2_s,'\n Indiv task (m,k) in each chain without deadline miss termination \n %6s %4s %4s %4s %8s \r\n','#set', '#chn', '#tsk', '#tid', '(m,k)');
    for i = 1:length(taskset)
        for j = 1:length(st_ts)
            mk_data = taskset{i}.("mk_s_ty2_"+string(j));
            for k = 1: length(st_ts{j})
                fprintf(fileID_2_s,'%6d %4d %4d %4d %4d %4d \r\n', i, j, k, st_ts{j}(k), mk_data{(taskset{i}{:,"ID"} == st_ts{j}(k))});
            end
        end
    end
end
fclose(fileID_2_s);

% Mothod 2 hyperperiod type 1
fileID_h1 = fopen('../output/cn_e2eltc_case17_1000_type1_hyp.txt', 'w');

fprintf(fileID_h1,'count of unschedulable, Mean and Median of min. util (m,k considered) per utilization \n');
fprintf(fileID_h1,'%8s %8s %8s \r\n', 'count', 'Mean', 'Median');
for util = 1: length(count_unsch_ts) 
        fprintf(fileID_h1,'%8d %8.4f %8.4f \r\n', count_unsch_ts(util, 3), mean_util(util, 3), median_util(util, 3));
end

fprintf(fileID_h1,'\n Ratio of chain completed per utilization \n');
fprintf(fileID_h1,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(p_ratio_h1)
    fprintf(fileID_h1,'%8.4f %8.4f %8.4f \r\n',p_ratio_h1(i,:));
end

fprintf(fileID_h1,'\n Ratio of chain delayed per utilization \n');
fprintf(fileID_h1,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(d_ratio_h1)
    fprintf(fileID_h1,'%8.4f %8.4f %8.4f \r\n',d_ratio_h1(i,:));
end

fprintf(fileID_h1,'\n Ratio of chain not-converged per utilization \n');
fprintf(fileID_h1,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(n_ratio_h1)
    fprintf(fileID_h1,'%8.4f %8.4f %8.4f \r\n',n_ratio_h1(i,:));
end
fclose(fileID_h1);

% Mothod 2 hyperperiod type 2
fileID_h2 = fopen('../output/cn_e2eltc_case17_1000_type2_hyp.txt', 'w');

fprintf(fileID_h2,'count of unschedulable, Mean and Median of min. util (m,k considered) per utilization \n');
fprintf(fileID_h2,'%8s %8s %8s \r\n', 'count', 'Mean', 'Median');
for util = 1: length(count_unsch_ts) 
        fprintf(fileID_h2,'%8d %8.4f %8.4f \r\n', count_unsch_ts(util, 4), mean_util(util, 4), median_util(util, 4));
end

fprintf(fileID_h2,'\n Ratio of chain completed per utilization \n');
fprintf(fileID_h2,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(p_ratio_h2)
    fprintf(fileID_h2,'%8.4f %8.4f %8.4f \r\n',p_ratio_h2(i,:));
end

fprintf(fileID_h2,'\n Ratio of chain delayed per utilization \n');
fprintf(fileID_h2,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(d_ratio_h2)
    fprintf(fileID_h2,'%8.4f %8.4f %8.4f \r\n',d_ratio_h2(i,:));
end

fprintf(fileID_h2,'\n Ratio of chain not-converged per utilization \n');
fprintf(fileID_h2,'%8s %8s %8s \r\n',' cn1', 'cn2', 'cn3');
for i = 1:length(n_ratio_h2)
    fprintf(fileID_h2,'%8.4f %8.4f %8.4f \r\n',n_ratio_h2(i,:));
end
fclose(fileID_h2);




% sum of execution time of each chain
% fileID = fopen('../output/cn_e2eltc_case13_1000.txt', 'w');
% fprintf(fileID,'\n Indiv Chain Sum of Exe. \n %8s %8s %8s %8s %8s \r\n','chain 1', 'chain 2', 'chain 3', 'chain 4', 'chain 5');
% for i = 1:length(taskset)
%     %fprintf(fileID,'%8.1f \r\n',sumC_cn(i,:));     % 1 chain
%     %fprintf(fileID,'%8.1f %8.1f \r\n',sumC_cn(i,1), sumC_cn(i,2));   % 2 chains
%     fprintf(fileID,'%8.1f %8.1f %8.1f %8.1f %8.1f \r\n',sumC_cn(i,1), sumC_cn(i,2), sumC_cn(i,3), sumC_cn(i,4), sumC_cn(i,5));   % 5 chains
% end
% 
% fclose(fileID);

%% WCRT : L-chain Latency along H-chain's Latency in dual-chain model
% % WARNING: CHANGE the FILE NAMES that save FIGURES from plotWCRTvsL_h.m
% for i = 1:length(taskset)
%     WCRT_values = plotWCRTvsL_h(cn1_T1(i), e2e_ltc_cn(i,1), sumC_cn(i,1), e2e_ltc_cn(i,2));
% end

%% graph of end-to-end latency by utilization, chain 
%newcolors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};

% figure(1); grid on; hold on;
% Tpos = linspace(1, 5, 5);   % data1 linspace
% xticks([1 2 3 4 5]);
% ylim([0 600])
% plot(Tpos, out_avg(:,1), '-o', 'Color', newcolors{1}, 'LineWidth', 2, 'MarkerSize', 10);
% saveas(gcf,'../figures/cn_e2eltc_case7_Tpos_avg.png');
% 
% 
% figure(2); grid on; hold on;
% s= swarmchart(maxT_pos, e2e_ltc_cn, 10, sumC_cn, "filled");
% colormap(jet);
% colorbar;
% xticks([1 2 3 4 5]);
% s.XJitter = 'rand';
% s.XJitterWidth = 0.5;
% ylabel('End-to-end latency'); xlabel('Position of longest T in chain');
% saveas(gcf,'../figures/cn_e2eltc_case7_Tpos_all.png');


%figure(1); grid on; hold on;
%U = linspace(0.60, 0.9, 7);   % data1 linspace
%plot(U, out_avg(:,1), '-o', 'Color', newcolors{1}, 'LineWidth', 2, 'MarkerSize', 10);

%lgn=legend('Utilization');
%lgn.FontSize = 12;
%ylabel('End-to-end latency'); xlabel('Utilization');

% save figure as a file
%saveas(gcf,'../figures/cn_e2eltc_case3_U.png');

% end-to-end latency by sum of execution time for each chain 
%figure(2); grid on; hold on;
%scatter(sumC_cn, e2e_ltc_cn, 'Color', newcolors{2});
%lgn=legend('Sum C per chain');
%lgn.FontSize = 12;
%ylabel('End-to-end latency'); xlabel('Sum of execution time per chain');
%saveas(gcf,'../figures/cn_e2eltc_case3_C.png');

% end-to-end latency by max of period for each chain 
%figure(3); grid on; hold on;
%scatter(maxT_cn, e2e_ltc_cn, 'Color', newcolors{3});
%lgn=legend('Max T per chain');
%lgn.FontSize = 12;
%ylabel('End-to-end latency'); xlabel('Max T per chain');
%saveas(gcf,'../figures/cn_e2eltc_case3_maxT.png');

% end-to-end latency by max T and sum C for each chain
%figure(4); grid on; hold on;
%scatter(maxT_cn, e2e_ltc_cn, 10, sumC_cn, "filled");
%colormap(jet);

%colorbar; % colorbar for the color scale
%ylabel('End-to-end latency'); xlabel('Max T per chain');

%saveas(gcf,'../figures/cn_e2eltc_case3_TandC.png');

% end-to-end latency by max T and its position in each chain
%figure(5); grid on; hold on;
%scatter(maxT_cn, e2e_ltc_cn, 10, maxT_pos, "filled");
%colormap(jet);

%colorbar; % colorbar for the color scale
%ylabel('End-to-end latency'); xlabel('Max T per chain');

%saveas(gcf,'../figures/cn_e2eltc_case3_Tandpos.png');




