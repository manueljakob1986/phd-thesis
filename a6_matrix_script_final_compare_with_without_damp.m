% Run matlab script a_matrix.m 
% Depending on switched off or on damping coefficent, run two script
% possibilities:
% 1. damping coefficent is switched off -> no exponential damping
% coefficient (e^0 = 1).
% 2. damping coefficent is switched on -> damping coefficient is computed
% based on damping ratio and undamped angular frequency (natural
% frequency).

% !!!!
% This script plots both cases (switched off and on damping
% coefficient)into one plot
% !!!!

% delete cache
clear all
clc
close all

% set global variable initial parking pricing to 2.5:
h1_setGlobal_initial_parking_pricing(2.5);

% set global variable for maximum price increase per time step to 0.5:  
h1_setGlobal_max_parking_price_increase(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parking_pricing_switched_on = c13_input_switch_on_parking_pricing;

if parking_pricing_switched_on == 1
    
    % Switch off the damping coefficient exponential factor in delta pricing
    % equation:
    % value 1 = "on", value 0 = "off"    
    h1_setGlobal_switch_on_damp_exp_coef(0);
    %   if damping coefficient exponential factor in delta pricing is switched
    %   off or on:

    [matrix, parking_pricing, ~, ~, ~, ~] = a_matrix(0);

    % plot Queuing diagram:
%     c_outputs_plots(matrix)
    % plot parking pricing:
    c_outputs_plot_parking_pricing(parking_pricing)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Switch on the damping coefficient exponential factor in delta pricing
    % equation:
    % value 1 = "on", value 0 = "off"
    h1_setGlobal_switch_on_damp_exp_coef(1);  

    % Run script with no exponential damping coefficient (e^0 = 1)   
    [~, parking_pricing, ~, ~, ~, ~] = a_matrix(0);

    % Computation of damping coefficient:  
    interpol = interp1(linspace(1,size(parking_pricing,1),size(parking_pricing,1)),parking_pricing,linspace(1,size(parking_pricing,1),size(parking_pricing,1)));
    [Wn,zeta] = damp(interpol);
    damping_coefficient = mean(zeta) * mean(Wn);
    damping_coefficient

    [matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, penalty_distance] = a_matrix(damping_coefficient);

    % plot Queuing diagram:
%     c_outputs_plots(matrix)
    % plot parking pricing:
    c_outputs_plot_parking_pricing(parking_pricing)

    % parking_pricing
    % guessed_price_vector
    % E_p_vot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif parking_pricing_switched_on == 0
    [matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, penalty_distance] = a_matrix(0);
    c_outputs_plots(matrix)
end    

save('load_a6_matrix_script_final_compare_with_without_damp')

