% Run Marcopark ALL social impacts parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maximum_value = 11;
ratio(:,1) = zeros((maximum_value + 1)*(maximum_value + 1),1);
average_searching_time(:,1) = zeros((maximum_value + 1)*(maximum_value + 1),1);
average_searching_distance(:,1) = zeros((maximum_value + 1)*(maximum_value + 1),1);
final_cum_revenue(:,1) = zeros((maximum_value + 1)*(maximum_value + 1),1);
occupancy_fuel(:,:) = zeros((maximum_value + 1)*(maximum_value + 1),1440);
occupancy_electric(:,:) = zeros((maximum_value + 1)*(maximum_value + 1),1440);
occupancy_all(:,:) = zeros((maximum_value + 1)*(maximum_value + 1),1440);
ACT(:,:) = zeros((maximum_value + 1)*(maximum_value + 1),1440);
optimal_target_occupancy_rate_all_vehicles(:,1) = zeros((maximum_value + 1)*(maximum_value + 1),1);
optimal_target_occupancy_rate_fuel_vehicles(:,1) = zeros((maximum_value + 1)*(maximum_value + 1),1);
optimal_target_occupancy_rate_electric_vehicles(:,1) = zeros((maximum_value + 1)*(maximum_value + 1),1);
social_impacts(:,1) = zeros((maximum_value + 1)*(maximum_value + 1),1);

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
demand_proportion = 10; %scenario b
supply_proportion = 10; %scenario b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
value_of_time = [    
     29.9/60
     25.4/60
     25.8/60
     17.2/60];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 i = 0;
for EV_infrastructure = 4000:1000:15000
    for VOT_multiplier = 0.25:0.25:3
         
    clearvars -except average_searching_time average_searching_distance final_cum_revenue occupancy_fuel occupancy_electric occupancy_all ACT optimal_target_occupancy_rate_all_vehicles optimal_target_occupancy_rate_fuel_vehicles optimal_target_occupancy_rate_electric_vehicles social_impacts policy demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale i ratio maximum_value value_of_time EV_infrastructure VOT_multiplier
    
    i = i + 1;
    h1_setGlobal_EV_infrastructure(EV_infrastructure);
    h1_setGlobal_VOT(value_of_time.*VOT_multiplier);
    [average_searching_time(i,1), average_searching_distance(i,1), final_cum_revenue(i,1), ...
        occupancy_fuel(i,:), occupancy_electric(i,:), occupancy_all(i,:), ACT(i,:),...
        optimal_target_occupancy_rate_all_vehicles(i,1), ...
        optimal_target_occupancy_rate_fuel_vehicles(i,1), ...
        optimal_target_occupancy_rate_electric_vehicles(i,1), social_impacts(i,1)] = ...
        Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options, supply_options, parking_duration_a_shape, parking_duration_b_scale);   
    
    EV_infrastructure
    VOT_multiplier

    close all
    end
end

save('results_all_social_impacts_policy_1.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.) Social Impacts - Contour plot:
figure
EV_infrastructure_variable = 4000:1000:15000;
VOT_multiplier_variable = 0.25:0.25:3;
plot_max = 11;
%CHANGE FOR POLICY
policy_1_avg_VOT = 24.6678;
policy_2_avg_VOT = 24.6775;
%CHANGE FOR POLICY
for k = 1:1:(maximum_value + 1)
    x_value_all(k,1) = EV_infrastructure_variable(k);
    %CHANGE FOR POLICY
    y_value_all(k,1) = policy_1_avg_VOT*VOT_multiplier_variable(k);
    %CHANGE FOR POLICY
    z_value_all(:,k) = social_impacts(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);

end
x_value = x_value_all(1:(plot_max + 1),1);
y_value = y_value_all(1:(plot_max + 1),1);
z_value = z_value_all(1:(plot_max + 1),1:(plot_max + 1));

newpoints = 500;
[xq,yq] = meshgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints ),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints )...
          );
z_value_interp = interp2(x_value,y_value,z_value,xq,yq,'cubic');
contourf(xq,yq,z_value_interp,50)
xlabel('Cost of establishing a parking space with charging facilities for an electric vehicle, \Phi')
ylabel('Average VOT based on proportion of vehicles by user group')
% axis([0,20,0,20])
set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
% set(cb, 'Limits', [1 317])
% caxis([1 317])
grid

