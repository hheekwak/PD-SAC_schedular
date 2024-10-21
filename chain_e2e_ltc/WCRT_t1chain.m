function R = WCRT_t1chain(exeC, maxDM, latency, chainT)
    % WCRT_t1chain calculates worst case respose time 
    %       in chain-scale in a taskset, with deadline miss termination
    %   INPUTS: vector of sum of tasks' execution time of each chain at the beginning
    %               - for every WCRT calculation, it is updated as a new execution time : exeC
    %           vector of maximum deadline miss (=c) of tasks in each chain: maxDM
    %           vector of maximum end-to-end latencis of each chain
    %               - for corresponding chain's execution time considering: latency
    %               - idle times and deadline misses are considered
    %           vector of period of each chain set with max latencies: chainT
    %   OUTPUT: vector of worst case response times R(:,1) and its 
    %           schedulability: R(:,2) 
    %           1: schedulable, 0: not schedulable, -1:RT is not converged  

    n = length(exeC);      % number of chains
    R = [zeros(n, 1), ones(n, 1)];
    schedulability = 1;
    for i = 1:n
        if i > 1 && (R(i-1, 2) == -1 || R(i-1, 1) == inf) % previous chain WCRT may not have converged
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
                added_DM = 0;
                % Interference from higher-priority chains
                for j = 1:i-1
                    interference = interference + ceil(Ri_old / chainT(j)) * exeC(j); 
                    added_DM = added_DM +  ceil(Ri_old / chainT(j)) * 2 * maxDM(i); % assume every interference brings 2 deadline misses 
                end
                % Update WCRT and update exeC
                Ri_new = latency(i) + interference + added_DM;
                exeC(i) = latency(i) + added_DM;
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

            R(i, 1) = Ri_new;               
            R(i, 2) = schedulability;
    
            if iter == max_iterations
                R(i, 1) = -1; 
                R(i, 2) = -1; % Maximum number of iterations reached. WCRT may not have converged.
            end
        end
    end    
end