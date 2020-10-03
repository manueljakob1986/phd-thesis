function [n_searching] = d12_transitions_n_dgp_s(n_enter_garage_total, n_want_to_garage_initially, Ndgp_k1, Ndgp_k2, Ndgp_k3, Ndgp_k4, A_3, initial_garage_capacity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This describes the new transition event:
% Number of vehicles that transition from “deciding to off-street parking” 
% to “searching for on-street parking” 
% during time slice i (Not access off-street parking).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_searching = zeros(1,4);
Ndgp_k = zeros(1,4);
Ndgp_k(1,1) = Ndgp_k1;
Ndgp_k(1,2) = Ndgp_k2;
Ndgp_k(1,3) = Ndgp_k3;
Ndgp_k(1,4) = Ndgp_k4;

total_garage_capacity = h2_getGlobal_parking_garage_capacity;

for k = 1:4
    if sum(Ndgp_k,2) == 0
        n_searching(1,k) = 0;
    else
        n_searching(1,k) = Ndgp_k(1,k)/sum(Ndgp_k,2) * (max(sum(n_want_to_garage_initially,2) - total_garage_capacity,0)) * A_3 / (initial_garage_capacity + A_3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Update total capacity:

if (total_garage_capacity - n_enter_garage_total) > 0
    % All vehicles "n_vehicles_dgp" enter garage and total capacity is reduced
    h1_setGlobal_parking_garage_capacity(total_garage_capacity - n_enter_garage_total);
    
else
    % The capacity is full.
    h1_setGlobal_parking_garage_capacity(0);
    
end

end