function [garage_parking_pricing] = d17_gp_pricing(parking_pricing, speed, ACT, input_value_of_time, C_drive, C_walk, availableparking_3, garage_parking_pricing_past_time_slice)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of garage parking pricing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_garage_capacity = h2_getGlobal_parking_garage_capacity;
K = size(input_value_of_time,1);

[~, input_parking_price_maximum_increase] = c10_input_maximum_parking_price_increase_per_time_unit;
off_street_parking_pricing_switched_on = c16_input_switch_on_off_street_parking_pricing;

if off_street_parking_pricing_switched_on == 1
    
    if total_garage_capacity < 1
        
        garage_parking_pricing = parking_pricing + c9_input_price_per_distance * speed * ACT;
        
        for k = 1:K
            garage_parking_pricing = garage_parking_pricing + 1/K * (input_value_of_time(k) * ACT - C_drive(1,k) - C_walk(1,k));
        end
        
    else
        
        garage_parking_pricing = availableparking_3 / total_garage_capacity * (parking_pricing + c9_input_price_per_distance * speed * ACT);
        
        for k = 1:K
            garage_parking_pricing = garage_parking_pricing + 1/K * (availableparking_3 / total_garage_capacity * input_value_of_time(k) * ACT - C_drive(1,k) - C_walk(1,k));
        end
    end
    
    if garage_parking_pricing - garage_parking_pricing_past_time_slice > input_parking_price_maximum_increase
        garage_parking_pricing = garage_parking_pricing_past_time_slice + input_parking_price_maximum_increase;
        
    elseif garage_parking_pricing_past_time_slice - garage_parking_pricing > input_parking_price_maximum_increase
        garage_parking_pricing = garage_parking_pricing_past_time_slice - input_parking_price_maximum_increase;
    end
    
    if garage_parking_pricing < c23_input_parking_price_garage
        garage_parking_pricing = c23_input_parking_price_garage;
    end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif off_street_parking_pricing_switched_on == 2
   
    % Constant garage parking pricing:
    garage_parking_pricing = c23_input_parking_price_garage; 
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif off_street_parking_pricing_switched_on == 0 
    
    % No garage parking pricing:
    garage_parking_pricing = 0;
       
end    

end