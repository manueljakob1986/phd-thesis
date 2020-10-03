function [use_maximum_limit, parking_price_maximum_increase] = c10_input_maximum_parking_price_increase_per_time_unit
% This function gives as an input the maximum parking price increase per
% time unit, i.e., the maximum price that is possible to increase per time
% step

parking_pricing_switched_on = c13_input_switch_on_parking_pricing;

if or(parking_pricing_switched_on == 1 , or(parking_pricing_switched_on == 2, parking_pricing_switched_on == 3))
% if use_maximum_limit is set to 1, then the maximum parking price
% increase or decrease functionality is on!
    use_maximum_limit = 1;

% define the maximum parking price increase or decrease 
% get from global variable in case script is used):
    parking_price_maximum_increase = h2_getGlobal_max_parking_price_increase;
%   parking_price_maximum_increase = 0.5;
%   parking_price_maximum_increase = 1;

elseif parking_pricing_switched_on == 0  
    use_maximum_limit = 0;
    parking_price_maximum_increase = 1000;
end

end