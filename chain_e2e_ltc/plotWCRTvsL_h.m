function WCRT_values = plotWCRTvsL_h(cnh_T1, L_h_max, L_h_opt, L_l)
    % plotWCRTvsL_h draws plot of L-chain Latency along H-chain's Latency
    %   INPUTS: T(period) of task 1 of H-chain: cnh_T1
    %           max latency of the higher priority chain: L_h_max 
    %           opt latency of the higher priority chain: L_h_opt 
    %           latency of the lower priority chain: L_l

    % Define a range of L_h values
    L_h_values = linspace(L_h_opt, L_h_max, 100);  % Adjust range and number of points as needed
    WCRT_values = zeros(size(L_h_values));  % Initialize WCRT values
    diff = zeros(size(L_h_values));        % for analysis

    % Calculate WCRT for each L_h value
    for i = 1:length(L_h_values)
        L_h = L_h_values(i);
        
        %Calculate T_h from L_h
        cnh_T = ceil(L_h/cnh_T1) * cnh_T1;       % H-chain's period, depending on its latency 
        
        diff(i) = cnh_T - L_h;            % how different Latency and Period of H-chain

        % Call the WCRT function with fixed L_l and current L_h
        WCRT_values(i) = WCRT_c(cnh_T, L_h, L_l);  
    end
 
     % Plot the results
     figure(1);
     plot(L_h_values, WCRT_values, '-o');
     ylim([0 700000])
     xlabel('Latency of Higher Priority Chain (L_h)');
     ylabel('Worst Case Response Time (WCRT)');
     title('WCRT vs. Latency of Higher Priority Chain (L_h)');
     grid on;

     saveas(gcf,'../figures/WCRTvsL_h_case7_shortest_1_100.png');
     
     figure(2);
     plot(L_h_values, diff, '-*');
     ylim([0 1000])
     xlabel('Ltc of H-chain');
     ylabel('how close to L-T');
     title('Difference between L-T of Higher Priority Chain');
     grid on;

     saveas(gcf,'../figures/WCRTvsL_h_case7_shortest_1_100_diff.png');
 end
 