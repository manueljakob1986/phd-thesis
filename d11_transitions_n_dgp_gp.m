function [n_enter_garage, n_enter_garage_total, n_want_to_garage_initially] = d11_transitions_n_dgp_gp(speed, t, ADD, n_vehicles_ns_dgp, n_vehicles_s_dgp, C_op, C_gp, C_drive, Ndgp_k1, Ndgp_k2, Ndgp_k3, Ndgp_k4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This describes the new transition event:
% Number of vehicles that transition from “deciding to off-street parking” 
% to “off-street parking” 
% during time slice i (Access off-street parking).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_garage_capacity = h2_getGlobal_parking_garage_capacity;

i = size(n_vehicles_ns_dgp,1);
distancetraveled = zeros(i,2);
n_ns_s_dgp = zeros(1,4);
n_want_to_garage_initially = zeros(1,4);
% % gamma_dgp_gp = zeros(1,4);
n_enter_garage = zeros(1,4);

for k = 1:4
       
    for j=1:i
        distancetraveled(j,1)=sum(speed(j:i,1))*t;
        distancetraveled(j,2)=sum(speed(j:i-1,1))*t;
        
        if and(distancetraveled(j,1) >= (ADD-0.7), distancetraveled(j,2) <= (ADD-0.7))
            n_ns_s_dgp(1,k) = n_ns_s_dgp(1,k) + n_vehicles_ns_dgp(j,k) + n_vehicles_s_dgp(j,k);
        end
    end
    
% %     % "There may be additional search costs (or cruising costs) involved if the
% %     % driving arrives at a garage and finds out that the pricing is higher than
% %     % what they expected, for example.”
% %     
% %     [a2, b2] = c18_input_parking_duration_distribution_parameter;
% %     % a2 = 2;
% %     % b2 = 5;
% %     %the average parking duration=a2*b2
% %     timecontrol_op = c14_input_on_street_parking_duration;
% %     % timecontrol_gp = c17_input_off_street_parking_duration;
% %     

    % vehicles with td <= tau_op that wanted to park initially in the garage:
    n_want_to_garage_initially(1,k) = n_ns_s_dgp(1,k);

% %     
% %     % decision variable to finally park in garage (due to on- and off-street
% %     % parking related costs):
% %     gamma_dgp_gp(1,k) = exp((C_op(1,k) + C_drive(1,k)) - (C_gp(1,k) - C_drive(1,k)))./...
% %         (exp((C_op(1,k) + C_drive(1,k)) - (C_gp(1,k) - C_drive(1,k))) + 1);
% % 
% %     
% %     if isnan(gamma_dgp_gp(1,k))
% %         gamma_dgp_gp(1,k) = 1;
% %     end
% %    
% %     
% %     n_ns_s_dgp(1,k) = (gamcdf(timecontrol_op,a2,b2) - gamcdf(0,a2,b2)) * ...
% %         gamma_dgp_gp(1,k) * n_ns_s_dgp(1,k) ...
% %         + (1-(gamcdf(timecontrol_op,a2,b2) - gamcdf(0,a2,b2))) * n_ns_s_dgp(1,k);
    
end

% Vehicles enter garage.
if sum(n_ns_s_dgp,2) <= total_garage_capacity
    
    % The vehicles go into garage and no vehicles go into "searching state":
    n_enter_garage_total = sum(n_ns_s_dgp,2);
    
else
    % The capacity is full.
    % Only the remaining "total_garage_capacity" cars can enter the garage,
    if sum(n_ns_s_dgp,2) == 0
        n_enter_garage_total = 0;
    else
        n_enter_garage_total = total_garage_capacity;
    end
    
end

if (Ndgp_k1 + Ndgp_k2 + Ndgp_k3 + Ndgp_k4) == 0
    n_enter_garage(1,1) = 0;
    n_enter_garage(1,2) = 0;
    n_enter_garage(1,3) = 0;
    n_enter_garage(1,4) = 0;
else
    n_enter_garage(1,1) = n_enter_garage_total * Ndgp_k1/(Ndgp_k1 + Ndgp_k2 + Ndgp_k3 + Ndgp_k4);
    n_enter_garage(1,2) = n_enter_garage_total * Ndgp_k2/(Ndgp_k1 + Ndgp_k2 + Ndgp_k3 + Ndgp_k4);
    n_enter_garage(1,3) = n_enter_garage_total * Ndgp_k3/(Ndgp_k1 + Ndgp_k2 + Ndgp_k3 + Ndgp_k4);
    n_enter_garage(1,4) = n_enter_garage_total * Ndgp_k4/(Ndgp_k1 + Ndgp_k2 + Ndgp_k3 + Ndgp_k4);
end

end