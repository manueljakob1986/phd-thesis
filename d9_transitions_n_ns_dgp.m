function [ngarage, ngarage_VOT] = d9_transitions_n_ns_dgp(matrix_ns_s_VOT, share_off_street_parking_decision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This describes the new transition event:
% Number of vehicles that transition from “non-searching” to “off-street parking” 
% during time slice i (Decide for off-street parking from “non-searching”).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a2, b2, ~, ~, ~, ~] = c18_input_parking_duration_distribution_parameter;
% a2 = 1.6;
% b2 = 142;
%the average parking duration=a2*b2
timecontrol_op = c14_input_on_street_parking_duration;
% timecontrol_gp = c17_input_off_street_parking_duration;

ngarage_VOT = zeros(1,4);
gamma_gp = zeros(1,4);
ngarage = 0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parking garage information (Change in August 2017, MJA):
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parking_garage_information = c24_input_parking_garage_information;
% switch_on_parking_garage_information = c25_input_switch_on_parking_garage_information;
% initial_garage_capacity = c19_input_number_of_parking_garages * c20_input_capacity_garage;
% total_garage_capacity = h2_getGlobal_parking_garage_capacity;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:4
    
    gamma_gp(1,k) = (gamcdf(timecontrol_op,a2,b2) - gamcdf(0,a2,b2)) * ...
        share_off_street_parking_decision(1,k) ...
        + (1-(gamcdf(timecontrol_op,a2,b2) - gamcdf(0,a2,b2)));   
    
    ngarage_VOT(1,k) = gamma_gp(1,k) * matrix_ns_s_VOT(1,k);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parking garage information (Change in August 2017, MJA):
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if switch_on_parking_garage_information == 1
%         ngarage_VOT(1,k) = ngarage_VOT(1,k) * ...
%             (parking_garage_information * total_garage_capacity / initial_garage_capacity ...
%             + (1 - parking_garage_information));
%     end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ngarage = ngarage + ngarage_VOT(1,k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end