function network_block_length = c21_input_network_block_length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This input parameter describes the average block length in the network.
% Assumption: Network is a square.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% 0.0833 km = 83.3333 m (block length)
% Based on formula: L = 2b*(b+1)*Bl
% 1 = 2*2*(2+1)*Bl
% network_block_length = 0.0833;

% 0.076 km = 76 m (block length)
network_block_length = 0.076;

end