function R = WCRT_c_h2_common_task(sumC, latency, chainT, hp_block_time, lp_block_time)
    % WCRT_c calculates worst case respose time in chain-scale in a taskset
    %   INPUTS: vector of sum of all task execution time of each chain: sumC
    %           vector of latency of each chain (without phsing delay): latency 
    %           vector of period of each chain set with max latencies and hyper period: chainT
    %           matrix of sum of blocking time of common tasks with specific higher periority chains: hp_block_time
    %                           row- higher(compared) priority chain idx, col- current chain idx
    %           vector of sum of blocking time of common tasks with lower priority chains: lp_block_time
    %   OUTPUT: vector of worst case response times R(:,1) and its 
    %           schedulability: R(:,2) 
    %           1: schedulable, 0: not schedulable, -1:RT is not converged  

    n = length(sumC);      % number of chains
    R = [zeros(n, 1), ones(n, 1)];
    schedulability = 1;
    for i = 1:n                           % i = current chain index
        if i > 1 && R(i-1, 2) == -1       % previous chain WCRT may not have converged
            R(i, 1) = -1;                 % rest may not be possible to be completed
            R(i, 2) = -1;

        else
            % Initial WCRT for a chain is its Latency
            Ri_old = 0;
            Ri_new= latency(i);
        
            % Iteratively calculate the WCRT of the chain
            max_iterations = 1000;  % Maximum number of iterations to avoid infinite loop
            tolerance = 1e-6;       % Convergence tolerance
            iter = 0;
        
            while abs(Ri_new - Ri_old) > tolerance && iter < max_iterations
                Ri_old = Ri_new;
                interference = 0;

                % Interference from higher-priority chains
                for j = 1:i-1             % j = higher priority chain index
                    interference = interference + ceil(Ri_old / chainT(j)) * sumC(j);
                end
                % Update WCRT
                Ri_new = latency(i) + interference + sum(hp_block_time(:,i));
                
                % Check if the response time exceeds the task's deadline
                if Ri_new > chainT(i)
                    schedulability = 0;     % cannot meet its deadline
                end
                if Ri_new == inf
                    Ri_new = -1;
                    schedulability = -1;
                    break;
                end
                iter = iter + 1;
            end

            % If WCRT < chain T, find WCRT with refined wcet (sumC) of the chain 
            if schedulability == 1
                % Initial WCRT for a chain is its sumC
                Ri_old = 0;
                Ri_new = sumC(i) + lp_block_time(i);
        
                % Iteratively calculate the WCRT of the chain
                iter = 0;
        
                while abs(Ri_new - Ri_old) > tolerance && iter < max_iterations
                    Ri_old = Ri_new;
                    interference = 0;
                    % Interference from higher-priority chains
                    for j = 1:i-1
                        interference = interference + ceil(Ri_old / chainT(j)) * sumC(j);
                    end
                    % Update WCRT
                    Ri_new = sumC(i) + lp_block_time(i) + interference + sum(hp_block_time(:,i));
                    
                    % Check if the response time exceeds the task's deadline
                    if Ri_new > chainT(i)
                        schedulability = 0;     % cannot meet its deadline
                    end
                    iter = iter + 1;
                end
            end

            R(i, 1) = Ri_new;               
            R(i, 2) = schedulability;
    
            if iter == max_iterations
                R(i, 1) = -1; 
                R(i, 2) = -1; % Maximum number of iterations reached. WCRT may not have converged.
            end
        end
    end    
end