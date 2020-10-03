% Run matlab script a_matrix.m and increase the denominator of the parking 
% pricing oscillation variable (decrease the parking pricing oscillation
% variable) as long as an optimal value of oscillation reduction in the
% parking pricing plot exists, i.e., the graph should converge to zero
% at the end of all time steps

% delete cache
clear all
clc
% close all

% Switch on the damping coefficient exponential factor in delta pricing
% equation:
% value 1 = "on", value 0 = "off"
h1_setGlobal_switch_on_damp_exp_coef(0);

% set global variable initial parking pricing to 2.5:
h1_setGlobal_initial_parking_pricing(2.5);

% set global variable for maximum price increase per time step to 0.1:  
% h1_setGlobal_max_parking_price_increase(0.5);
h1_setGlobal_max_parking_price_increase(0.1);

parking_pricing_switched_on = c13_input_switch_on_parking_pricing;

if or(parking_pricing_switched_on == 1 , parking_pricing_switched_on == 3)

[matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, cost_travelling_next_space, penalty_distance, total_costs, bias] = a_matrix(0);
a = size(parking_pricing,1);

% plot Queuing diagram:
c_outputs_plots(matrix)
% plot parking pricing:
c_outputs_plot_parking_pricing(parking_pricing)
% plot revenue:
c_outputs_plot_revenue(matrix,parking_pricing)

% new_parking_price_oscillation_denominator
% parking_pricing(a,1)
% parking_pricing
% guessed_price_vector
% E_p_vot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif parking_pricing_switched_on == 0
    [matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, cost_travelling_next_space, penalty_distance, total_costs, bias] = a_matrix(0);
    c_outputs_plots(matrix)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif parking_pricing_switched_on == 2

[matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, cost_travelling_next_space, penalty_distance, total_costs, bias] = a_matrix(0);
a = size(parking_pricing,1);

% plot Queuing diagram:
c_outputs_plots(matrix)
% plot parking pricing:
c_outputs_plot_parking_pricing(parking_pricing)
% plot revenue:
c_outputs_plot_revenue(matrix,parking_pricing)

% new_parking_price_oscillation_denominator
% parking_pricing(a,1)
% parking_pricing
% guessed_price_vector
% E_p_vot
   
end    

b_outputs_from_the_matrix

save('load_a1_matrix_script')

