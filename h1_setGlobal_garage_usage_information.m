function h1_setGlobal_garage_usage_information(val)
% This is a help function to set the global variable
% 'global_garage_usage_information':

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in February 2018 to include parking garage capacity
% information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This input parameter describes the parking garage capacity
% information, switch it on or off.
% parking garage capacity information switched on = 1
% parking garage capacity information switched off = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global global_garage_usage_information
global_garage_usage_information = val;

end