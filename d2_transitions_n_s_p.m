function [naccess, parking_probability, guessed_price, tau, VOT_k1, VOT_k2, VOT_k3, VOT_k4, E_p_vot, travel_distance_cost, penalty_distance, ...
    total_costs, final_parking_fee, final_delta_searching_veh_available_parking_spots, ...
    delta_searching_veh_available_parking_spots_i_plus_one] = ...
    d2_transitions_n_s_p(A, N, L, v, t, matrix_ns_s_VOT, available_parking_spots, searching_veh, v_all_t, parking_price_oscillation, A_total, ...
    old_parking_pricing, old_delta_searching_veh_available_parking_spots, Ns_3_k1, Ns_3_k2, Ns_3_k3, Ns_3_k4, decideforparking_3, starttosearch_3)
% N=4;    %searchers
% A=3;    %available parking
% v=30;   %km/h
% L=1;    %network size
% t=1/60; %time slice

% if N<=A
%     if v*t/L<=1
%         naccess=N*(1-(1-v*t/L)^A);
%     else
%         naccess=N;
%     end
% end
% 
% if N>A
%     if v*t/L<=1/N
%         naccess=N*(1-(1-v*t/L)^A);
%     elseif v*t/L<=A/N
%         naccess=(N-A-N*(1-1/N)^A)/(1-A)*(N*v*t/L-A)+A;
%     else
%         naccess=A;
%     end
% end

% -------------------------------------------------------------------------

a = size(searching_veh,1);

parking_pricing_switched_on = c13_input_switch_on_parking_pricing;

if parking_pricing_switched_on == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Change, Paper Review, February 2017  
    
    only_occupancy = 0;
    % call function to compute demand-responsive parking pricing
    [~, delta_searching_veh_available_parking_spots, final_delta_searching_veh_available_parking_spots] = ...
        d7_parking_pricing_demand_responsive(searching_veh, available_parking_spots, parking_price_oscillation, ...
        old_parking_pricing, old_delta_searching_veh_available_parking_spots, only_occupancy);
  
    % call function to compute parking probability
    [final_parking_fee, parking_probability, guessed_price, tau, VOT_k1, VOT_k2, VOT_k3, VOT_k4, E_p_vot, travel_distance_cost, ...
        penalty_distance, total_costs, delta_searching_veh_available_parking_spots_i_plus_one] = ...
                d8_parking_pricing_probability(old_parking_pricing, searching_veh, available_parking_spots, delta_searching_veh_available_parking_spots, ...
                matrix_ns_s_VOT, v_all_t, t, A, L, N, Ns_3_k1, Ns_3_k2, Ns_3_k3, Ns_3_k4, decideforparking_3, starttosearch_3, parking_price_oscillation, A_total, 0);           
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif parking_pricing_switched_on == 0 
    parking_probability = zeros(1,4);
    guessed_price = 0;
    tau = 0; 
    VOT_k1 = 0;
    VOT_k2 = 0;
    VOT_k3 = 0;
    VOT_k4 = 0;
    E_p_vot = 0;
    penalty_distance = 0; 
    total_costs = 0;
    final_parking_fee = zeros(size(searching_veh,1),1);
    travel_distance_cost = 0;
    delta_searching_veh_available_parking_spots_i_plus_one = 0;
    final_delta_searching_veh_available_parking_spots = 0;
    
elseif parking_pricing_switched_on == 2
    
    only_occupancy = 0;
    [~, delta_searching_veh_available_parking_spots, final_delta_searching_veh_available_parking_spots] = ...
        d7_parking_pricing_demand_responsive(searching_veh, available_parking_spots, parking_price_oscillation, ...
        old_parking_pricing, old_delta_searching_veh_available_parking_spots, only_occupancy);
    
    % Constant on-street parking pricing:
    parking_pricing = c8_input_parking_price * ones(size(searching_veh,1),1);
    
    % call function to compute parking probability
    [final_parking_fee, parking_probability, guessed_price, tau, VOT_k1, VOT_k2, VOT_k3, VOT_k4, E_p_vot, travel_distance_cost, ...
        penalty_distance, total_costs, delta_searching_veh_available_parking_spots_i_plus_one] = ...
                d8_parking_pricing_probability(parking_pricing, searching_veh, available_parking_spots, delta_searching_veh_available_parking_spots, ...
                matrix_ns_s_VOT, v_all_t, t, A, L, N, Ns_3_k1, Ns_3_k2, Ns_3_k3, Ns_3_k4, decideforparking_3, starttosearch_3, parking_price_oscillation, A_total, 1);  

elseif parking_pricing_switched_on == 3
    
    only_occupancy = 1;
    [~, delta_searching_veh_available_parking_spots, final_delta_searching_veh_available_parking_spots] = ...
        d7_parking_pricing_demand_responsive(searching_veh, available_parking_spots, parking_price_oscillation, ...
        old_parking_pricing, old_delta_searching_veh_available_parking_spots, only_occupancy);   
    
    % call function to compute parking probability
    [final_parking_fee, parking_probability, guessed_price, tau, VOT_k1, VOT_k2, VOT_k3, VOT_k4, E_p_vot, travel_distance_cost, ...
        penalty_distance, total_costs, delta_searching_veh_available_parking_spots_i_plus_one] = ...
                d8_parking_pricing_probability(old_parking_pricing, searching_veh, available_parking_spots, delta_searching_veh_available_parking_spots, ...
                matrix_ns_s_VOT, v_all_t, t, A, L, N, Ns_3_k1, Ns_3_k2, Ns_3_k3, Ns_3_k4, decideforparking_3, starttosearch_3, parking_price_oscillation, A_total, 2);            
            
end    
    
if N <= A
    if v*t/L <= 1
        naccess = (N*(1 - (1 - v*t/L)^A));

    else
        naccess = N;
    end
end

if N > A
    if v*t/L <= 1/N
        naccess = (N*(1 - (1 - v*t/L)^A));

    elseif v*t/L <= A/N
        naccess = ((N - A - N*(1 - 1/N)^A) / (1 - A)*(N*v*t/L - A) + A);

    else
        naccess = A;
    end
end

end