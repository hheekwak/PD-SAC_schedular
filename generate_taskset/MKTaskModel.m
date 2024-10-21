classdef MKTaskModel
    % (m, k) task model
    % - priority: fixed / dynamic priorities, use a missed model in order
    % to exploit response time analysis
    
    properties
        id, T, C, D, m, k, l, s, w, meet_idx, miss_cnt, cur_C, timer, ready, priorities, ...
          p_time, p_stack, exe_slot, miss_slot, rt_task, util_idx, util, ...
          mk_ratio, mk_idx, failure, consecutive_miss, act_util, mode, WCRT, offset, ...
          sporadic, spo_timer, spo_n_init, concat_release, jitter
    end
    
    methods
        % Initialize parameters
        function obj = MKTaskModel(id, T, C, D, m, k, offset, sporadic, jitter)
            if nargin > 0
                obj.id = id;
                obj.m = m;
                obj.k = k;
                s = '';
                for i = 1 : k
                    s = strcat(s, 'N');
                end
                obj.s = s;
                obj.T = T;
                obj.C = C;
                obj.D = D;
                obj.l = m + 1;    % Number of job-classes
                obj.w = max(floor(k/(k-(k-m))) - 1, 1);
                obj.cur_C = 0;
                obj.meet_idx = 1;
                obj.miss_cnt = 0;
                obj.timer = zeros(1, obj.l);
                obj.ready = 0;
                obj.priorities = zeros(1, m+1);
                obj.p_time = [];
                obj.p_stack = [];
                obj.exe_slot = [];
                obj.miss_slot = [];
                obj.rt_task = {};
                obj.util_idx = [];
                obj.util = [];
                obj.mk_ratio = 0;
                obj.mk_idx = 1;
                obj.failure = 0;
                obj.consecutive_miss = 0;
                obj.WCRT = zeros(1, m+1);
                obj.act_util = 0;
                obj.mode = 1;
                obj.offset = offset;
                obj.sporadic = sporadic;
                obj.spo_timer = 'N';
                obj.spo_n_init = [];
                obj.concat_release = [];
                obj.jitter = jitter;
            end
        end
        
        function obj = setTimer(obj, idx)
            obj.timer(idx) = obj.T;
        end
        
        function obj = exeTimer(obj, idx, scale)
            obj.timer(idx) = obj.timer(idx) - scale;
        end
        
        function obj = updateFailure(obj)
            m = 0;
            for i = 1 : obj.k
                if obj.s(i) == '0'
                    m = m + 1;
                end
            end
            
            if m > obj.k - obj.m
                obj.failure = obj.failure + 1;
%                disp(['Failure task: ', num2str(obj.id)]);
            end
        end
        
        function obj = updateMKratio(obj)
            m = 0;
            for i = 1 : obj.k
                if obj.s(i) == '1'
                    m = m + 1;
                end
            end
            
            if m ~= 0
                tmp = obj.m/m;

                obj.mk_ratio = (obj.mk_idx*obj.mk_ratio + tmp)/(obj.mk_idx + 1);
                obj.mk_idx = obj.mk_idx + 1;
            else
                
            end
        end
        
        function obj = updateConsecutiveMiss(obj)
            m = 0; flag = 0; max_m = 0;
            for i = 1 : obj.k
                if obj.s(i) == '0' && flag == 1
                    m = m + 1;
                elseif obj.s(i) == '0' && flag == 0
                    m = m + 1;
                    flag = 1;
                else
                    if m > max_m
                        max_m = m;
                    end
                    flag = 0;
                    m = 0;
                end
            end
            
            if isempty(find(obj.s == '1'))
                obj.consecutive_miss = obj.consecutive_miss + 1;
            else
                if max_m > (obj.k - obj.m)
                    obj.consecutive_miss = obj.consecutive_miss + 1;
                end
            end
        end
        
        function obj = updateTrace(obj, complete)
            obj.s(1:end-1) = obj.s(2:end);
            if complete == 1
                obj.s(end) = '1';
            elseif complete == 0
                obj.s(end) = '0';
            end
            obj = updateMKratio(obj);
            obj = updateFailure(obj);
            obj = updateConsecutiveMiss(obj);
        end
        
        function obj = updateMeetIdx(obj)            
            % Modified Mode 1 - check from the lastest sequence
            tmp_s = obj.s(2:end);
            p_cnt = 0;
            for i = 1 : length(tmp_s)
                if tmp_s(length(tmp_s) - i + 1) == '1'
                    p_cnt = p_cnt + 1;
                else
                    break;
                end
            end
            
            % Count the latest consecutive deadline misses
            n_cnt = 0;
            for i = 1 : length(tmp_s)
                if tmp_s(length(tmp_s) - i + 1) == '0'
                    n_cnt = n_cnt + 1;
                else
                    break;
                end
            end
             
            if n_cnt == 0
                obj.meet_idx = p_cnt + 1;
            else
                if n_cnt == obj.w % Threshold to change priority to the first job-class
                    obj.meet_idx = 1;
                end 
            end
            
            if obj.meet_idx >= length(obj.priorities)
                obj.meet_idx = length(obj.priorities);
            elseif obj.meet_idx <= 1
                obj.meet_idx = 1;
            end
        end
        
        function obj = completeProcess(obj, t, rtask)
            if obj.ready < 1
                disp('debug');
            end
            obj.timer(obj.ready) = 0; % unset timer      
            obj.cur_C = 0;   % clear current exe
            obj.ready = 0;   % clear ready param
            obj = updateTrace(obj, 1);
            
            obj = updateMeetIdx(obj);
            
            obj.rt_task(rtask.child_id) = rtask;
        end
        
        function obj = terminateProcess(obj, t, rtask)
            obj = updateTrace(obj, 0);
            obj = updateMeetIdx(obj);
            obj.cur_C = 0;
            obj.ready = 0;
            
            obj.miss_cnt = obj.miss_cnt + 1;            
            obj.rt_task(rtask.child_id) = rtask;
        end

        function obj = concateReleaseTime(obj)
            for i = 1 : length(obj.rt_task)
                obj.concat_release = [obj.concat_release obj.rt_task(i).init_act obj.rt_task(i).no_sporadic];
%                 disp(obj.concat_release);
            end
        end
    end
    
end
