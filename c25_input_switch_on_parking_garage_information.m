function switch_on_parking_garage_information = c25_input_switch_on_parking_garage_information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in August 2017 to include parking garage capacity
% information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This input parameter describes the parking garage capacity
% information, switch it on or off.
% parking garage capacity information switched on = 1
% parking garage capacity information switched off = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% switch_on_parking_garage_information = 0;
switch_on_parking_garage_information = h2_getGlobal_garage_usage_information;

end