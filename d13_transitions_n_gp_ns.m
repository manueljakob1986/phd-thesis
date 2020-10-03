function ndepart_garage = d13_transitions_n_gp_ns(n_enter_garage, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This describes the new transition event:
% Number of vehicles that transition from “garage parking” to “non-searching” 
% during time slice i (Depart garage parking).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = size(n_enter_garage,1);
ndepart_garage = 0;
total_garage_capacity = h2_getGlobal_parking_garage_capacity;

for k = 1:4
    for i = 1:a
        td=(a+1-i)*60*t;% unit: minutes;
        n(i,k) = n_enter_garage(i,k) * c3_input_parkingduration(td,t,2);
        
        % n(i,1) is the number of vehicles who arrived in time slice i and leaving
        % in time slice a, in other words, the cars who arrived in time slice i and
        % stayed for a period of td.
               
        ndepart_garage = ndepart_garage + n(i,k);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Update total capacity:

if (total_garage_capacity + ndepart_garage) <= c19_input_number_of_parking_garages * c20_input_capacity_garage
    h1_setGlobal_parking_garage_capacity(total_garage_capacity + ndepart_garage); 
else
    h1_setGlobal_parking_garage_capacity(c19_input_number_of_parking_garages * c20_input_capacity_garage); 
end

end