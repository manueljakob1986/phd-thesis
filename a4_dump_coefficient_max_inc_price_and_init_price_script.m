% This script plots the relationship between the dumping coefficient
% (oscillation parameter), the maximum pricing factor and the initial parking pricing

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

% inputs:
% define the accuracy of the converge of the optimal value for the
% oscillation reduction in the parking pricing plot: epsilon
epsilon = 0.5; 
% maximum interval for maximum price increase:
max_interval_x = 2.5;
% maximum interval for initial parking pricing:
max_interval_y = 3.5;
i = 1;
j = 1;
dumping_coefficient_vector = zeros(max_interval_x*10,max_interval_y*10 - 4);

for maximum_price_increase = 0.1:0.1:max_interval_x
    
%   set global variable for maximum price increase:  
    h1_setGlobal_max_parking_price_increase(maximum_price_increase);

    for initial_parking_pricing = 0.5:0.1:max_interval_y
%       set global variable for initial parking pricing:  
        h1_setGlobal_initial_parking_pricing(initial_parking_pricing);

        new_dumping_coefficient_denominator = 1; %Highest oscillation reduction parameter
        dumping_coefficient = c11_input_parking_price_oscillation(1 / new_dumping_coefficient_denominator);
        [matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, penalty_distance] = a_matrix(dumping_coefficient);
        a = size(parking_pricing,1);

        index = 0;
        while parking_pricing(a,1) > c8_input_parking_price + epsilon
%           delete cache of variable parking_pricing
            clear parking_pricing;
    
%           set new denominator of the parking pricing oscillation variable
            new_dumping_coefficient_denominator = new_dumping_coefficient_denominator + 1;

            dumping_coefficient = c11_input_parking_price_oscillation(1 / new_dumping_coefficient_denominator);
            [matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, penalty_distance] = a_matrix(dumping_coefficient);
            a = size(parking_pricing,1);
            
            
            if index > 500
                fprintf('Iterations stopped')
                fprintf('\n Parking pricing in last iteration:')
                parking_pricing(a,1)
                fprintf('\n for dumping coeeficient denominator:')
                new_dumping_coefficient_denominator
                break
            end
            index = index + 1;
        end

        dumping_coefficient_vector(i,j) = dumping_coefficient;
        j = j + 1;
        
        initial_parking_pricing
    end
    
    maximum_price_increase
    i = i + 1;
end

% [X,Y] = meshgrid(linspace(0.1,max_interval_x,max_interval_x*10), linspace(0.5,max_interval_y,max_interval_y*10 - 4));
% X
% Y
x_vector = linspace(0.1,max_interval_x,max_interval_x*10)
y_vector = linspace(0.5,max_interval_y,max_interval_y*10 - 4)
dumping_coefficient_vector

figure
% contour3(X,Y,dumping_coefficient_vector(:,:))
contour3(x_vector,y_vector,dumping_coefficient_vector(:,5:max_interval_y*10))

title('Dumping coefficient for different maximum pricing increase limits and initial parking pricing values')
xlabel('Maximum pricing increase limit')
ylabel('Initial parking pricing')
zlabel('Dumping coefficient') 

toc

save('load_a4_dump_coefficient_max_inc_price_and_init_price_script')

