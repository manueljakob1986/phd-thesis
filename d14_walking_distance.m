function [walking_distance_gp, walking_distance_op] = d14_walking_distance(L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation of average walking distance from parking garage to final
% destination

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
walking_distance_gp = 2/3 * (c21_input_network_block_length / sqrt(pi * c19_input_number_of_parking_garages)) *...
   (-1/2 + sqrt(1/4 + L/(2*c21_input_network_block_length)));

walking_distance_op = 2/3 * (c21_input_network_block_length *(-1/2 + sqrt(1/4 + L/(2*c21_input_network_block_length))));

end