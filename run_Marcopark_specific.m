% Run Marcopark Specific:
% Specific proportions of electric vehicles in demand and 
% parking spaces for electric vehicles in supply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes the parking space policy we would like 
% to use.
% policy = 1: Vehicle type dependent parking spaces
% policy = 2: No parking space restrictions for some vehicle types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
policy = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes whether we would like to consider
% the initial demand (demand_options = 0) or
% otherwise demand decrease or increase in % (e.g., 1.05 or 0.95)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demand_options = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes whether we would like to consider
% the initial supply (supply_options = 0) or
% otherwise supply decrease or increase in % (e.g., 1.05 or 0.95)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
supply_options = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes whether we would like to consider
% the initial parking duration (parking_duration_a_shape = 1 or 
% parking_duration_b_scale = 1) or otherwise parking duration decrease or 
% increase in % (e.g., 1.05 or 0.95)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parking_duration_a_shape = 1;
parking_duration_b_scale = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
value_of_time = [    
     29.9/60
     25.4/60
     25.8/60
     17.2/60];
h1_setGlobal_EV_infrastructure(10000);
h1_setGlobal_VOT(value_of_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set proportions of electric vehicles in demand and 
% parking spaces for electric vehicles in supply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% demand_proportion = 5;
% supply_proportion = 45;

% demand_proportion = 0; %scenario a
% supply_proportion = 0; %scenario a

% demand_proportion = 10; %scenario 2.1
% supply_proportion = 7; %scenario 2.1

% demand_proportion = 10; %scenario 2.2
% supply_proportion = 8; %scenario 2.2

% demand_proportion = 10; %scenario 2.3
% supply_proportion = 9; %scenario 2.3

demand_proportion = 10; %scenario b
supply_proportion = 10; %scenario b

% demand_proportion = 10; %scenario 4.1
% supply_proportion = 11; %scenario 4.1

% demand_proportion = 10; %scenario 4.2
% supply_proportion = 12; %scenario 4.2

% demand_proportion = 10; %scenario 4.3
% supply_proportion = 13; %scenario 4.3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[average_searching_time, average_searching_distance, final_cum_revenue, ...
    occupancy_fuel, occupancy_electric, occupancy_all, ACT,...
    optimal_target_occupancy_rate_all_vehicles, ...
    optimal_target_occupancy_rate_fuel_vehicles, ...
    optimal_target_occupancy_rate_electric_vehicles, social_impacts] = ...
    Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options, supply_options, parking_duration_a_shape, parking_duration_b_scale);

social_impacts
