% This script plots the relationship between the dumping coefficient
% (oscillation parameter) and the initial parking pricing

% The goal is to find an analytical formula for the dumping coefficient

tic
% delete cache
clear all
clc
% close all

% Switch on the damping coefficient exponential factor in delta pricing
% equation:
% value 1 = "on", value 0 = "off"
h1_setGlobal_switch_on_damp_exp_coef(1);

% set global variable for maximum price increase per time step to 0.5:  
h1_setGlobal_max_parking_price_increase(0.5);

% inputs:
% define the accuracy of the converge of the optimal value for the
% oscillation reduction in the parking pricing plot: epsilon
epsilon = 0.5; 
max_interval = 4;
i = 1;
dumping_coefficient_vector = zeros(max_interval*10,1);

for initial_parking_pricing = 0.1:0.1:max_interval
%   set global variable for initial parking pricing:  
    h1_setGlobal_initial_parking_pricing(initial_parking_pricing);

    new_dumping_coefficient_denominator = 1; %Highest oscillation reduction parameter
    dumping_coefficient = c11_input_parking_price_oscillation(1 / new_dumping_coefficient_denominator);
    [matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, penalty_distance] = a_matrix(dumping_coefficient);
    a = size(parking_pricing,1);

    while parking_pricing(a,1) > c8_input_parking_price + epsilon
%       delete cache of variable parking_pricing
        clear parking_pricing;
    
%       set new denominator of the parking pricing oscillation variable
        new_dumping_coefficient_denominator = new_dumping_coefficient_denominator + 1;

        dumping_coefficient = c11_input_parking_price_oscillation(1 / new_dumping_coefficient_denominator);
        [matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, penalty_distance] = a_matrix(dumping_coefficient);
        a = size(parking_pricing,1);
    end

    dumping_coefficient_vector(i,1) = dumping_coefficient;
    i = i + 1;
    initial_parking_pricing
end

% Initial parking pricing from 0.1 to max_interval
figure
plot(linspace(0.5,max_interval,max_interval*10),dumping_coefficient_vector(:,1))
 
title('Dumping coefficient for different initial parking pricing values')
xlabel('Initial parking pricing')
ylabel('Dumping coefficient') 

% Initial parking pricing from 0.5 to max_interval
figure
plot(linspace(0.5,max_interval,max_interval*10 - 4),dumping_coefficient_vector(5:max_interval*10,1))

title('Dumping coefficient for different initial parking pricing values')
xlabel('Initial parking pricing')
ylabel('Dumping coefficient') 
toc

save('load_a3_dump_coefficient_initial_price_script')

