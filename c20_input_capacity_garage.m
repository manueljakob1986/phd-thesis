function capacity_garage = c20_input_capacity_garage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This input parameter describes the capacity of a parking garage.
% It is the initial capacity for time step i = 0.
% Assumption: Same capacity is assumed for all garages in the network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% capacity_garage = 166;

capacity_garage = h2_getGlobal_capacity_garage;

end