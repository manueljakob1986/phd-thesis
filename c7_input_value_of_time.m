function value_of_time = c7_input_value_of_time

%     value_of_time = [    
%          29.9/60
%          25.4/60
%          25.8/60
%          17.2/60];

value_of_time = h2_getGlobal_VOT;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% total number of different origins
% V = 3;

%  parking_pricing_switched_on = c13_input_switch_on_parking_pricing;

% value of time per hour/60 (= minute):
% value of times for different origins, dependent on V:

% if or(parking_pricing_switched_on == 1 , parking_pricing_switched_on == 2)
%     value_of_time = [    
%          29.9/60
%          25.4/60
%          25.8/60
%          17.2/60];
% elseif parking_pricing_switched_on == 0     
%     value_of_time = 1;
% end
 
end