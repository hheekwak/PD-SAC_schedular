%% Experiments
% Generate tasksets with one task TC-fixed, in different position in taskset
% Designed exact N tasksets with sz_taskset different positions of longest 
% Hard real time tasksets
clc; close all; clear all;

sets = {};
N = 5; % Number of tasksets of specific position of fixed task
sz_taskset =  5; % Size of a taskset
fix_util = 0.8;
fix_task_T = 500;
oth_T_range = [10 100];

epsilon = 0.001; % A small positive value to ensure C is not zero
upperC = fix_task_T * 0.2; % Limit upper bound of C as 20 percent of T 

for i = 1 : N
    % Generate N tasksets
        % fixed task, longest T 
        fix_task_C = epsilon + (upperC - epsilon) * rand();
        res_util = fix_util - (fix_task_C/fix_task_T);
        fix_task = [fix_task_T fix_task_C fix_task_T 0 1 0 0 0];

        % num of tasks, T range, K range, Maximum utilization
        N_task = sz_taskset - 1;
        [Data_raw, util_M, util_m, valid] = Generate_Taskset(N_task, oth_T_range, [1,1], res_util, -1);
    num = 0;
    % Distribute longest T in different positions
    for j = 1 : sz_taskset
        if j == 1
            Data = [fix_task; Data_raw];
        elseif j == sz_taskset
            Data = [Data_raw; fix_task];
        else
            Data = [Data_raw(1:j-1, :); fix_task; Data_raw(j:end, :)];
        end

        ss = size(Data);
        % Define tasks
        for k = 1 : ss(1, 1)
            task(k) = MKTaskModel(k, Data(k,1), Data(k,2), Data(k,3), Data(k,4), Data(k,5), Data(k,6), Data(k,7), Data(k,8));
        end     

    % sets(i, j)
        if valid == true 
            sets(i, j).Data = Data;
            sets(i, j).util = util_M;
            %sets(i, j).util_M = util_M;
            %sets(i, j).util_m= util_m;
            num = num + 1; 
            sets(i, j).num = num;
        end 
    end
 end

% save sets as a txt file
ss = size(sets);
fileID = fopen('../data/chain_case6_5samples.txt', 'w');
fileNum = fopen('../data/chain_case6_5samples_util.txt', 'w');

for j = 1 : ss(1,2)
    max_util = 0;     % max and min utilizations of 1000 samples
    min_util = 1;
    for i = 1 : ss(1,1)
        data = sets(i, j).Data;
        num = sets(i, j).num;
        util = sets(i, j).util;
        max_util = max(max_util, util);
        min_util = min(min_util, util);
        %max_util = sets(i, j).util_M;
        %min_util = sets(i, j).util_m;
        fprintf(fileID, '%6.1f %7.3f %6.1f %d %d %6.6f %6.6f %6.6f\r\n', data');
        if num == 1
            fprintf(fileNum, '%4d %6.6f\n', i, util);
        end
    end
    fprintf(fileNum, '\n=============================\n');
    fprintf(fileNum, '%d %6.6f %6.6f\n', num, max_util, min_util);
end
