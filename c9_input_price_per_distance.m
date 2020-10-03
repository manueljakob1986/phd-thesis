function price_per_distance = c9_input_price_per_distance

% price per distance (km)
% value in CHF:
 
parking_pricing_switched_on = c13_input_switch_on_parking_pricing;

if or(parking_pricing_switched_on == 1 , or(parking_pricing_switched_on == 2, parking_pricing_switched_on == 3 ))
    price_per_distance = 0.3;
elseif parking_pricing_switched_on == 0  
    price_per_distance = 0;
end   
 
end