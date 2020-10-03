function h1_setGlobal_garage_usage_information_parameter(val)
% This is a help function to set the global variable
% 'h1_setGlobal_garage_usage_information_parameter':

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in March 2018 to include parking garage capacity
% information parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This input parameter describes the parking garage capacity
% information. Availability of garage parking information, i.e., 
% percentage of drivers who know the available garage capacity at their 
% time of decision between on- and off-street parking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global global_garage_usage_information_parameter
global_garage_usage_information_parameter = val;

end