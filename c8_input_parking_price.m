function initial_parking_price = c8_input_parking_price

% on-street parking price 
% reference value for demand-responsive pricing, fixed at interation start;
% value in CHF:
% Also first demand-responsive pricing iteration recomputes first value,
% with searching vehicles and available parking spots.

% Assumption: Fixed price is equal for all on-street parking spots:

 parking_pricing_switched_on = c13_input_switch_on_parking_pricing;
 
if parking_pricing_switched_on == 1
%   initial_parking_price = 4.5;
    initial_parking_price = h2_getGlobal_initial_parking_pricing;
elseif parking_pricing_switched_on == 0   
    initial_parking_price = 0;   
elseif parking_pricing_switched_on == 2 
    initial_parking_price = h2_getGlobal_initial_parking_pricing;
elseif parking_pricing_switched_on == 3
    initial_parking_price = h2_getGlobal_initial_parking_pricing;
end     
 
end