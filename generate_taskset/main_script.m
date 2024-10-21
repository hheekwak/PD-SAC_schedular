%% Experiments
% Generate (Hard real time) task sets
clc; close all; clear all;

sets = {};
N = 1000; % Number of tasksets of specific utilization (usually 1000)
N_range = [15 15]; % The range of tasks in a taskset
u_start = 0.3; u_finish = 0.9; unit = 0.05;
M = ceil((u_finish - u_start)/unit);
if M == 0 M = 1; end
U = linspace(u_start, u_finish, M);
iter = 0;
for u = 1 : M
    tmp_util_M = 0; tmp_util_m = 0; valid_task = 0; flag = 1; i = 1;
    cnt = 0;
    k = 0;
    while flag
        rm_flag = 1;
        iter = iter + 1;
        valid_task = 1;
        
        % Generate task set
        % num of tasks, T range, K range, Maximum utilization
        N_task = randi(N_range, 1);
   %     [Data_raw, util_M, util_m, valid] = Generate_Taskset(N_task, [10 1000], [1,1], U(u), -1); % T with range
        [Data_raw, util_M, util_m, valid] = Generate_Taskset(N_task, [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000], [1,1], U(u), -1); % T with options
        Data = Data_raw; % not sorted, if sort needed, use below line
        %Data = sortrows(Data_raw, [1, 2], {'ascend' 'descend'});
        ss = size(Data);
        
        % Define tasks
        for j = 1 : ss(1, 1)
            task(j) = MKTaskModel(j, Data(j,1), Data(j,2), Data(j,3), Data(j,4), Data(j,5), Data(j,6), Data(j,7), Data(j,8));
        end     

        % below commented out for chain-aware scheduler
        %{
        wcrt = WCRT_rm(task);
        
        for j = 1 : ss(1, 1)
            if strcmp(wcrt(j).schedulable, 'un-schedulable')
                rm_flag = 0;
                valid_task = 0;
            end
        end
        %} 
        
        if rm_flag == 1
            cnt = cnt + 1;
        end
        
        if valid == true %&& valid_task == 1
            sets(u, i).util_M = util_M;
            sets(u, i).util_m = util_m;
            tmp_util_M = tmp_util_M + util_M;
            tmp_util_m = tmp_util_m + util_m;
            sets(u, i).Data = Data(:,1:4);
            sets(u, i).num_tasks = size(Data, 1);
            i = i + 1;
        end
        
        if i == N + 1
            flag = 0;
        end
        k = k + 1;
    end
    rm_sched(u) = cnt/N;
    Max_util(u) = tmp_util_M/N;
    Min_util(u) = tmp_util_m/N;
    
    % disp(['Max Utilization of ',num2str(U(u)),': ',num2str(tmp_util_M/N)]);
    % disp(['Min Utilization of ',num2str(U(u)),': ',num2str(tmp_util_m/N)]);
end

% save sets as a txt file
ss = size(sets);
fileID = fopen('../data/chain_case13_1.txt', 'w');
fileNum = fopen('../data/chain_case13_1_num.txt', 'w');

for u = 1 : ss(1,1)
    for i = 1 : ss(1,2)
        data = sets(u, i).Data;
        index_column = (1:size(data, 1))';
        data = [data(:, 1:3), index_column, data(:, 4:end)];
        num = sets(u, i).num_tasks;
        max_util = sets(u, i).util_M;
        min_util = sets(u, i).util_m;
        fprintf(fileID, '%6.1f %7.4f %6.1f %4d %d \r\n', data');
        fprintf(fileNum, '%d %6.6f %6.6f\n', num, max_util, min_util);
    end
end
