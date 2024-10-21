function R = WCRT_t2chain(sumC, latency, chainT)
    % WCRT_t2chain calculates worst case respose time 
    %       in chain-scale in a taskset, without deadline miss termination
    %   INPUTS: vector of sum of tasks' execution time of each chain: sumC
    %               - for higher-priority chain's execution time
    %           vector of maximum end-to-end latencis of each chain: latency
    %               - for corresponding chain's execution time considering
    %               idle times
    %           vector of period of each chain set with max latencies: chainT
    %   OUTPUT: vector of worst case response times R(:,1) and its 
    %           schedulability: R(:,2) 
    %           1: schedulable, 0: not schedulable, -1:RT is not converged  

    n = length(sumC);      % number of chains
    R = [zeros(n, 1), ones(n, 1)];
    schedulability = 1;
    for i = 1:n
        if i > 1 && R(i-1, 2) == -1;      % previous chain WCRT may not have converged
            R(i, 1) = -1;                 % rest may not be possible to be completed
            R(i, 2) = -1;

        else
            % Initial WCRT for a chain is its Latency
            if i == 1               % execution time of highest chain is sum of tasks execution time
                Ri_old = sumC(i);   % idle time is not considered for highest-priority chain
                Ri_new = sumC(i);   
            else 
                Ri_old = 0;
                Ri_new= latency(i);
            end
        
            % Iteratively calculate the WCRT of the chain
            max_iterations = 1000;  % Maximum number of iterations to avoid infinite loop
            tolerance = 1e-6;       % Convergence tolerance
            iter = 0;
        
            while abs(Ri_new - Ri_old) > tolerance && iter < max_iterations
                Ri_old = Ri_new;
                interference = 0;
                % Interference from higher-priority chains
                for j = 1:i-1
                    interference = interference + ceil(Ri_old / chainT(j)) * sumC(j);
                end
                % Update WCRT
                Ri_new = latency(i) + interference;
                
                % Check if the response time exceeds the task's deadline
                if Ri_new > chainT(i)
                    schedulability = 0;     % cannot meet its deadline
                end
                iter = iter + 1;
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