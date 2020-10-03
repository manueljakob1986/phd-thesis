function [n_vehicles_dgp] = d10_transitions_n_s_dgp(Ns_3_k1, Ns_3_k2, Ns_3_k3, Ns_3_k4, share_off_street_parking_decision, decideforparking_3_VOT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This describes the new transition event:
% Number of vehicles that transition from “searching” to “deciding for off-street parking” 
% during time slice i (Switch to off-street parking).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parking garage information (Change in August 2017, MJA):
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parking_garage_information = c24_input_parking_garage_information;
% switch_on_parking_garage_information = c25_input_switch_on_parking_garage_information;
% initial_garage_capacity = c19_input_number_of_parking_garages * c20_input_capacity_garage;
% total_garage_capacity = h2_getGlobal_parking_garage_capacity;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global_s_dgp_penalty = h2_getGlobal_s_dgp_penalty;
penalty_term = min(1, 1/(Ns_3_k1 + Ns_3_k2 + Ns_3_k3 + Ns_3_k4)^global_s_dgp_penalty);

n_vehicles_dgp(1,1) = share_off_street_parking_decision(1,1) * (Ns_3_k1 - decideforparking_3_VOT(1,1)) * penalty_term;
n_vehicles_dgp(1,2) = share_off_street_parking_decision(1,2) * (Ns_3_k2 - decideforparking_3_VOT(1,2)) * penalty_term;
n_vehicles_dgp(1,3) = share_off_street_parking_decision(1,3) * (Ns_3_k3 - decideforparking_3_VOT(1,3)) * penalty_term;
n_vehicles_dgp(1,4) = share_off_street_parking_decision(1,4) * (Ns_3_k4 - decideforparking_3_VOT(1,4)) * penalty_term;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parking garage information (Change in August 2017, MJA):
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if switch_on_parking_garage_information == 1
%     
%     n_vehicles_dgp(1,1) = n_vehicles_dgp(1,1) * ...
%         (parking_garage_information * total_garage_capacity / initial_garage_capacity ...
%         + (1 - parking_garage_information));
%     
%     n_vehicles_dgp(1,2) = n_vehicles_dgp(1,2) * ...
%         (parking_garage_information * total_garage_capacity / initial_garage_capacity ...
%         + (1 - parking_garage_information));
%     
%     n_vehicles_dgp(1,3) = n_vehicles_dgp(1,3) * ...
%         (parking_garage_information * total_garage_capacity / initial_garage_capacity ...
%         + (1 - parking_garage_information));
%     
%     n_vehicles_dgp(1,4) = n_vehicles_dgp(1,4) * ...
%         (parking_garage_information * total_garage_capacity / initial_garage_capacity ...
%         + (1 - parking_garage_information));
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end