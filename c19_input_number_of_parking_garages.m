function number_of_parking_garages = c19_input_number_of_parking_garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This input parameter describes the number of a parking garages in the 
% network.
% Assumption: All parking garages are homogeneously distributed over the
% ring network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number_of_parking_garages = 2;

number_of_parking_garages = h2_getGlobal_number_of_parking_garages;

end