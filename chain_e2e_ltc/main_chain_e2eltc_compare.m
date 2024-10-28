%% Experiments
% Parameters
clc; close all; clear;

is_timing = 0;     % Process time measurement: when 1, comment out Line 147- 192)
is_compare = 0;    % Compare reaction time difference between PD-SAC ann TA-CEC 

% Bring task-counts and chain-counts of each taskset as arrays 
fileID_sz = fopen('../data/rtas2021data/102124_t_ch_each_ts_180_1000.txt', 'r');
line = fgets(fileID_sz);                    % Read a task count line from the file as a string
ts_t_size = eval(line);                     % Convert the string to a MATLAB array
line = fgets(fileID_sz);                    % Read a chain count line from the file as a string
ts_ch_size = eval(line);                    % Convert the string to a MATLAB array
fclose(fileID_sz);                          % Close the file

% number of tasksets for each case (e.g. utilization, max T position...)
n_ts = length(ts_t_size);  

% Bring chain configurations in 
fileID = fopen('../data/rtas2021data/102124_chains_180_1000.txt', 'r');
% Initialize a cell array to store the arrays
st_ts = {};

% Read each line one by one

for set = 1: n_ts
    for chn = 1: ts_ch_size(set)
    % Read a line from the file
    line = fgetl(fileID);
    
    % Split the line by white space 
    splitLine = strsplit(line);
    
    % Convert the second (third) to last items to double, build a chain structure
    cn_st = str2double(splitLine(3:end));  % It counts starting white space as an item 
    
    % Store the chain configuration in the structure cell array
    st_ts{set, chn} = cn_st;
    end
end


%% data read
% also consider re-nameing names of result and figure file 

% read data
data1 = load('../data/rtas2021data/102124_taskset_180_1000.txt');
data1 = array2table(data1, 'VariableNames', {'T' 'C' 'D' 'ID' 'Prior'});

% add ED, upperbound of potential delay from its release point to execution
% add empty column for I, upperbound of idle time of each task

data1.ED = data1.T - data1.C;
data1.I(:) = NaN;

%% variables

% set split taskset cell and result arrays
taskset = cell(length(ts_t_size), 1);
prio_chain = cell(size(st_ts));
sorted_org_idx = cell(length(ts_t_size), 1);
sumC_cn = zeros(length(taskset), max(ts_ch_size));          % use longest row of st_ts

e2e_ltc_wo_DM = NaN(length(taskset), max(ts_ch_size));      % Type 2: end-to-end latency in a single chain, including ED, not DM 
e2e_ltc_wo_DM_ED = NaN(length(taskset), max(ts_ch_size));   % Type 2: end-to-end latency in a single chain, no ED, no DM considered 
T_cn_t2 = zeros(length(taskset), max(ts_ch_size));          % Type 2: set T of each chain with its max latency and T_1 without Deadline miss termination 
WCRT_chain_2 = NaN(length(taskset), max(ts_ch_size));       % Type 2: WCRT without deadline miss termination 
%schd_able_2 = NaN(length(taskset), max(ts_ch_size));        % Type 2: schedulability of each chain 

react_t = NaN(length(taskset), max(ts_ch_size));            % Type 2: without deadline miss termination 
utilization = NaN(1, length(ts_t_size));                    % Type 2: minimum utilization for each taskset

%% calculate
% Start timing
tic;
prev_ts_ln = 0;
for i = 1: height(taskset)
    taskset(i) = {data1(prev_ts_ln + 1: prev_ts_ln + ts_t_size(i),["ID" "T" "C" "ED" "I"])};
    prev_ts_ln = prev_ts_ln + ts_t_size(i);
    % set up arrays to save cumulative interval for each chain
    cum_int_wo_dm = NaN(1, ts_ch_size(i));    % Type 2: arr for cumulative interval without deadline miss for latency (exclude idle time of task1)
    max_idle_in_chn = NaN(1, ts_ch_size(i));  % Type 2: arr for longest max idle in a chain for (m,k) (include idle time of task1)

    
    % Sort the chain and make a priority of chains, also save the original chain index order for comparison to rtas2021
    % Step 0: Copy current task's chain configurations 
    i_chains = st_ts(i, 1:ts_ch_size(i));
    % Step 1: Sort tasks in individual chains within the cell array
    t_sorted_chains = cellfun(@sort, i_chains, 'UniformOutput', false);
    % Step 2: Save the original indices
    org_idx = 1:length(i_chains);
    % Step 3: Sort by values first and then by length
    % Pad arrays with -inf to make them the same length for comparison but shorter has higher priority  
    [~, sorted_idx] = sortrows(cell2mat(cellfun(@(x) [x, -inf(1, max(cellfun(@length, t_sorted_chains))-length(x))], t_sorted_chains, 'UniformOutput', false)'));
    % Step 4: Prioritize the original chain composition using the sorted indices
    prio_chain(i,1:ts_ch_size(i)) = i_chains(sorted_idx);
    % Step 5: Store the original indices in the order of the priority
    sorted_org_idx{i} =  org_idx(sorted_idx);

    % calculate sum of execution time for each chain & max idle time for each task
    for j = 1: ts_ch_size(i)
        for k = 1: length(prio_chain{i, j})     % ith taskset, jth chain configuration e.g.[3 2 5 6 7 8 1]
            sumC_cn(i, j) = sumC_cn(i, j) + taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(k)), "C"};
        end    

        % max idle time for each task
        for k = 1: length(prio_chain{i, j}) 
            if taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(k)), "T"} > sumC_cn(i, j)
                taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(k)), "I"} = taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(k)), "T"} - sumC_cn(i, j);
            else
                taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(k)), "I"} = 0;
            end
        end

        % Cumulatively calculated maximal interval of each chain
        % range of cumulative interval, from nth and to nth task in a chain
        % the order in a chain = index of chain structure
        st_order = 1;
        nd_order = length(prio_chain{i, j});
        subset_cn = taskset{i}(ismember(taskset{i}.ID, prio_chain{i, j}), ["ID" "C" "I" "ED"]);
        
        cum_int_wo_dm(j) = max(taskset{i}{ismember(taskset{i}.ID, prio_chain{i, j}(2:end)), "I"});  % Type 2: idle time of task1 doesn't affect to end-to-end latency
        max_idle_in_chn(j) = max(taskset{i}{ismember(taskset{i}.ID, prio_chain{i, j}), "I"});       % Type 2: longest max idle time in a chain for m,k calculation

        % Latency of a chain 
        %       e2e_ltc_wo_DM:    Type 2 end-to-end latency as a single chain 
        %                               used for setting chain T for WCRT, execution delay included, deadline miss termination not considered
        %       e2e_ltc_wo_ED:    Type 2 latency between execution initiation of start-task and competion of end-task
        %                               used for setting chain C for WCRT (execution delay, deadline miss) not considered 
        
        e2e_ltc_wo_DM(i, j) = subset_cn{(subset_cn{:,"ID"} == prio_chain{i, j}(st_order)), "ED"} + sumC_cn(i, j) + cum_int_wo_dm(j);   
        e2e_ltc_wo_DM_ED(i, j) = sumC_cn(i, j) + cum_int_wo_dm(j);         

        % Set T(period) of each chain from its max e2e latency with and without DM and T of task1
        cn_T1 = taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(1)), "T"};           % T_1 of each chain 
        T_cn_t2(i, j) = ceil(e2e_ltc_wo_DM(i, j)/cn_T1) * cn_T1;    % Type 2: without deadline miss termination
        
    end


    % WCRT with T of chains and C (latency or sumC) of chains
    % Type2 (No deadline miss termination) chainT from max e2e latency and T of task1
    R_2 = WCRT_t2chain(sumC_cn(i,:), e2e_ltc_wo_DM_ED(i,:), T_cn_t2(i,:));
    WCRT_chain_2(i,:) = R_2(:,1)';
    %schd_able_2(i,:) = R_2(:,2)';
 
% WHEN TIMIMG, COMMENT OUT FROM Line 147- 192 (current line is 146)

    % Find lenient (m,k) of each task within each chain in multi-chain model
    for j = 1: ts_ch_size(i)
        % Find max(m,k) of each task with WCRT in a multi-chain model
        % Type 2: deadline miss termination NOT considered
        taskset{i}.("mk_m_ty2_"+string(j)) = getmk_in_multi(taskset{i}, prio_chain{i, j}, sumC_cn(i,j), WCRT_chain_2(i,j), T_cn_t2(i,j));
    end

    % Collect tasks that does not consist any chain
    solo_t_id = [taskset{i}.ID(isnan(taskset{i}.I))];

    % Find the minimum utilization of taskset
    sumUtil = 0;
    % Singly executing tasks (no m,k considered)
    for s = 1:length(solo_t_id)
        T = taskset{i}{(taskset{i}{:,"ID"} == solo_t_id(s)), "T"};
        C = taskset{i}{(taskset{i}{:,"ID"} == solo_t_id(s)), "C"};
        sumUtil = sumUtil + (C/T);
    end
    % Chainly executing tasks (lanient m,k applied)
    for j = 1: ts_ch_size(i)
        for k = 1: length(prio_chain{i, j}) 
            T = taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(k)), "T"};
            C = taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(k)), "C"};
            K = taskset{i}{(taskset{i}{:,"ID"} == prio_chain{i, j}(k)), "mk_m_ty2_"+string(j)}{1}(2);
            % new utilization is not available if it is not schedulable
            if K == -1
                sumUtil = -1;
                break;
            else
                sumUtil = sumUtil + ((C/T) * 1/K);
            end
        end
    end
    utilization(i) = sumUtil;

    % Reaction Time in this model 
    % This is for comparison between Type 2 multi-chain model and rtas2021
  
    for j = 1: ts_ch_size(i)
        if WCRT_chain_2(i,j) >= T_cn_t2(i, j)
            react_t(i, j) = 2 * WCRT_chain_2(i,j);
        else
            react_t(i, j) = T_cn_t2(i, j) + WCRT_chain_2(i,j);
        end
    end 

end
% End timing
elapsedTime = toc;
if is_timing == 1
    fileID_time = fopen('../output/102124_elapsedTime_180_1000.txt', 'w');
    fprintf(fileID_time,'Elapsed time: %8.4f', elapsedTime);
end

% Count the number of -1s in the array
count_unsch_ts = sum(utilization == -1);
% Filter out unschedulable taskset utilization value (-1)
filtered_util = utilization(utilization ~= -1);
% Calculate the average of schedulable taskset's utilization
average_util = mean(filtered_util);
% Fine the median value of schedulable taskset's utilization
median_util = median(filtered_util);

%% Compare reaction time of PDSAC scheduler and TA_CEC
if (is_timing == 0 & is_compare == 1)
    % Bring chain configurations in 
    fileID_cp = fopen('../data/rtas2021data/102124_reaction_t_180_1000.txt', 'r');
    % Initialize an array to store reaction time
    ta_cec = NaN(size(react_t));
    
    % Read each line one by one
    for set = 1: n_ts
        array = [];
        for chn = 1: ts_ch_size(set)
            % Read a line from the file
            line = fgetl(fileID_cp);
            
            % Split the line by white space 
            splitLine = strsplit(line);
            
            % Convert the second (third) to double, store the value
            array = [array, str2double(splitLine(3))];  % It counts starting white space as an item 
        end
        ta_cec(set, 1:ts_ch_size(set)) = array(sorted_org_idx{set});
    end
    
    % Calculate the difference between corresponding reaction times
    difference = react_t - ta_cec;
    % Convert the differences to a 1D vector
    diff_vector = difference(:);
    % Remove NaN values
    diff_vector = diff_vector(~isnan(diff_vector));
    % Create the box plot
    boxplot(diff_vector);
    title('Box Plot of Differences Between Values (NaNs Removed)');
    ylabel('Difference Value');
    saveas(gcf, '102124_boxplot_difference_180_1000.png'); % Save as PNG file
end
%% save result as txt files
if is_timing == 0
    % For Data Analysis Comparison To RTAS 2021 Reaction time and WCRT
    fileID = fopen('../output/102124_cn_e2eltc_compare_180_1000.txt', 'w');
    fprintf(fileID,'\n Count of not schedulable tasksets: %d ', count_unsch_ts);
    fprintf(fileID,'\n Mean of utilization: %8.3f ', average_util);
    fprintf(fileID,'\n Median of utilization: %8.3f ', median_util);
    fprintf(fileID,'\n Indiv. Chain Reaction time and WCRT in multi-chain model without deadline miss termination \n ');
    fprintf(fileID,'%8s %8s %8s %8s \r\n', 'chnPrio', 'ReactT', 'WCRT', 'origID');
    for i = 1:length(taskset)
        fprintf(fileID, 'Min. utilization: %.3f \r\n', utilization(i));
        for j = 1:ts_ch_size(i)
            fprintf(fileID,'%8d %8.3f %8.3f %8d \r\n', j, react_t(i,j), WCRT_chain_2(i,j), sorted_org_idx{i}(j));
        end
        fprintf(fileID,'\r\n');
    end
    fclose(fileID);
    
    % (m,k) for each task (multi, type 2)
    fileID_2_m = fopen('../output/102124_cn_e2eltc_compare_mk_180_1000.txt', 'w');
    
    fprintf(fileID_2_m, '\n Indiv task (m,k) in a multi-chain model without deadline miss termination \n %6s %4s %4s %4s %8s \r\n','#set', '#chn', '#tsk', '#tid', '(m,k)');
    for i = 1:length(taskset)
        for j = 1:ts_ch_size(i)
            mk_data = taskset{i}.("mk_m_ty2_"+string(j));
            for k = 1: length(prio_chain{i, j})
                fprintf(fileID_2_m,'%6d %4d %4d %4d %4d %4d \r\n', i, j, k, prio_chain{i, j}(k), mk_data{(taskset{i}{:,"ID"} == prio_chain{i, j}(k))});
            end
        end
    end
    
    fclose(fileID_2_m);
end