function [walking_distance] = d10_walking_distance(L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation of average walking distance from parking garage to final
% destination

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
walking_distance = 2/3 * (c18_input_network_block_length *(-1/2 + sqrt(1/4 + L/(2*c18_input_network_block_length))));

end