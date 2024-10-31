function [ sets, util_M, util_m, valid ] = Generate_Taskset(num, t_range, k_range, util, mode)
    % mode (1) : apply same (m, k) constraints for all tasks
    % mode (-1): apply hard real-time (m,k)=(0,1) constraints for all tasks
    % check if t_range is range to open pick(1) or array of specific numbers(2)

    
    % T / C / D / m / k / offset / sporadic / jitter
    sets = zeros(num, 5);
    
%     T = (t_range(2)-t_range(1)).*round(rand([1 num]), 1)+t_range(1);
%    T = randi(t_range, num, 1); % t_range pick is open from the range (1)
     T = randsample(t_range, num, true); % t_range pick from given numners(2) 

    sets(:, 1) = T';
    sets(:, 3) = T';
    
    % Limit task utilization
%     flag = 1;
%     while flag
%         util_per = UUniFast(num, util);
%         if util_per < 0.5
%             break;
%         end
%     end
    
    util_per = UUniFast(num, util); % randomly generated utilization of each task
    valid = true;
    if ~isempty(find(util_per > 1.0))
%         disp('Not correct generation');
        valid = false;
    end
    util_M = 0; util_m = 0;
    
    K = randi(k_range, 1, 1);
    M = randi([1 9], 1, 1);
    for i = 1 : num
        sets(i, 2) = round(util_per(i)*sets(i, 1),3); 
        if sets(i, 2) < 0.001
            valid = false;      % if execution time is 0, it is not valid
            break;
%           disp('Not correct generation:');
        end
        if mode == 0
            sets(i, 5) = randi(k_range, 1, 1);
            sets(i, 4) = randi([1 k_range(1)-1], 1, 1); % setting random m 
%             sets(i, 4) = randi([1 2], 1, 1);
%             sets(i, 4) = randi([1 ceil(sets(i,5)*0.5)], 1, 1);
%             sets(i, 5) = 2;
%             sets(i, 4) = 1;
        elseif mode == -1       % hard real-time
            sets(i, 5) = 1;
            sets(i, 4) = 0;
        else
            sets(i, 5) = K;
            sets(i, 4) = M;
        end
        sets(i, 6) = 0;
        sets(i, 7) = 0;
        sets(i, 8) = 0;
        
        util_M = util_M + sets(i, 2)/sets(i, 1);
        if mode ~= -1
            util_m = util_m + sets(i, 2)/sets(i, 1)*sets(i, 4)/sets(i, 5);
        else
            util_m = util_M;
        end
    end
    
end

