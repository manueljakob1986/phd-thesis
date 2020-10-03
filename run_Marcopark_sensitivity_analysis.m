% Run Marcopark for sensitivity analysis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maximum_value = 9;
number_sensitivity_parameters = 3;
maximum_proportion_value = 20;

average_searching_time_d(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
average_searching_time_s(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
average_searching_distance_d(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
average_searching_distance_s(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
final_cum_revenue_d(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
final_cum_revenue_s(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
optimal_target_occupancy_rate_all_vehicles_d(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
optimal_target_occupancy_rate_all_vehicles_s(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
optimal_target_occupancy_rate_fuel_vehicles_d(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
optimal_target_occupancy_rate_fuel_vehicles_s(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
optimal_target_occupancy_rate_electric_vehicles_d(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
optimal_target_occupancy_rate_electric_vehicles_s(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
social_impacts_d(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);
social_impacts_s(:,:,:) = zeros(maximum_value,number_sensitivity_parameters,maximum_proportion_value + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes the parking space policy we would like 
% to use.
% policy = 1: Vehicle type dependent parking spaces
% policy = 2: No parking space restrictions for some vehicle types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
policy = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set proportions of electric vehicles in demand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demand_proportion = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set proportions of parking spaces for electric vehicles in supply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% supply_proportion = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
value_of_time = [    
     29.9/60
     25.4/60
     25.8/60
     17.2/60];
h1_setGlobal_EV_infrastructure(10000);
h1_setGlobal_VOT(value_of_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity parameter: k
% 0 = no parameter used
% 1 = demand
% 2 = supply
% 3 = average parking duration (shape parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demand entering the area (k = 1)
k = 1; %flexible demand
% sensitivity_parameter = k;
for demand_proportion = 0:1:maximum_proportion_value
    for demand_options = -0.4:0.1:0.4

    clearvars -except average_searching_time_d average_searching_time_s average_searching_distance_d average_searching_distance_s final_cum_revenue_d final_cum_revenue_s optimal_target_occupancy_rate_all_vehicles_d optimal_target_occupancy_rate_all_vehicles_s optimal_target_occupancy_rate_fuel_vehicles_d optimal_target_occupancy_rate_fuel_vehicles_s optimal_target_occupancy_rate_electric_vehicles_d optimal_target_occupancy_rate_electric_vehicles_s social_impacts_d social_impacts_s policy demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale i k maximum_value maximum_proportion_value

        supply_proportion = 10;
        supply_options = 0;
        parking_duration_a_shape = 1;
        parking_duration_b_scale = 1;
        i = int16(demand_options/0.1) + 5;
        [average_searching_time_d(i,k,demand_proportion + 1), average_searching_distance_d(i,k,demand_proportion + 1), final_cum_revenue_d(i,k,demand_proportion + 1), ...
            ~, ~, ~, ~, optimal_target_occupancy_rate_all_vehicles_d(i,k,demand_proportion + 1), ...
            optimal_target_occupancy_rate_fuel_vehicles_d(i,k,demand_proportion + 1), ...
            optimal_target_occupancy_rate_electric_vehicles_d(i,k,demand_proportion + 1),...
            social_impacts_d(i,k,demand_proportion + 1)] = ...
            Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options + 1, supply_options, parking_duration_a_shape, parking_duration_b_scale);

        close all

    end
end
clearvars demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demand entering the area (k = 1)
k = 1; %flexible supply
% sensitivity_parameter = k;
for supply_proportion = 0:1:maximum_proportion_value
    for demand_options = -0.4:0.1:0.4

    clearvars -except average_searching_time_d average_searching_time_s average_searching_distance_d average_searching_distance_s final_cum_revenue_d final_cum_revenue_s optimal_target_occupancy_rate_all_vehicles_d optimal_target_occupancy_rate_all_vehicles_s optimal_target_occupancy_rate_fuel_vehicles_d optimal_target_occupancy_rate_fuel_vehicles_s optimal_target_occupancy_rate_electric_vehicles_d optimal_target_occupancy_rate_electric_vehicles_s social_impacts_d social_impacts_s policy demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale i k maximum_value maximum_proportion_value

        demand_proportion = 10;
        supply_options = 0;
        parking_duration_a_shape = 1;
        parking_duration_b_scale = 1;
        i = int16(demand_options/0.1) + 5;
        [average_searching_time_s(i,k,supply_proportion + 1), average_searching_distance_s(i,k,supply_proportion + 1), final_cum_revenue_s(i,k,supply_proportion + 1), ...
            ~, ~, ~, ~, optimal_target_occupancy_rate_all_vehicles_s(i,k,supply_proportion + 1), ...
            optimal_target_occupancy_rate_fuel_vehicles_s(i,k,supply_proportion + 1), ...
            optimal_target_occupancy_rate_electric_vehicles_s(i,k,supply_proportion + 1),...
            social_impacts_s(i,k,supply_proportion + 1)] = ...
            Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options + 1, supply_options, parking_duration_a_shape, parking_duration_b_scale);

        close all

    end
end
clearvars demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total number of existing parking spaces in the area (k = 2)
k = 2; %flexible demand
% sensitivity_parameter = k;
for demand_proportion = 0:1:maximum_proportion_value
    for supply_options = -0.4:0.1:0.4

    clearvars -except average_searching_time_d average_searching_time_s average_searching_distance_d average_searching_distance_s final_cum_revenue_d final_cum_revenue_s optimal_target_occupancy_rate_all_vehicles_d optimal_target_occupancy_rate_all_vehicles_s optimal_target_occupancy_rate_fuel_vehicles_d optimal_target_occupancy_rate_fuel_vehicles_s optimal_target_occupancy_rate_electric_vehicles_d optimal_target_occupancy_rate_electric_vehicles_s social_impacts_d social_impacts_s policy demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale i k maximum_value maximum_proportion_value

        supply_proportion = 10;
        demand_options = 0;
        parking_duration_a_shape = 1;
        parking_duration_b_scale = 1;
        i = int16(supply_options/0.1) + 5;
        [average_searching_time_d(i,k,demand_proportion + 1), average_searching_distance_d(i,k,demand_proportion + 1), final_cum_revenue_d(i,k,demand_proportion + 1), ...
            ~, ~, ~, ~, optimal_target_occupancy_rate_all_vehicles_d(i,k,demand_proportion + 1), ...
            optimal_target_occupancy_rate_fuel_vehicles_d(i,k,demand_proportion + 1), ...
            optimal_target_occupancy_rate_electric_vehicles_d(i,k,demand_proportion + 1),...
            social_impacts_d(i,k,demand_proportion + 1)] = ...
            Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options, supply_options + 1, parking_duration_a_shape, parking_duration_b_scale);

        close all

    end
end
clearvars demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total number of existing parking spaces in the area (k = 2)
k = 2; %flexible supply
% sensitivity_parameter = k;
for supply_proportion = 0:1:maximum_proportion_value
    for supply_options = -0.4:0.1:0.4

    clearvars -except average_searching_time_d average_searching_time_s average_searching_distance_d average_searching_distance_s final_cum_revenue_d final_cum_revenue_s optimal_target_occupancy_rate_all_vehicles_d optimal_target_occupancy_rate_all_vehicles_s optimal_target_occupancy_rate_fuel_vehicles_d optimal_target_occupancy_rate_fuel_vehicles_s optimal_target_occupancy_rate_electric_vehicles_d optimal_target_occupancy_rate_electric_vehicles_s social_impacts_d social_impacts_s policy demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale i k maximum_value maximum_proportion_value

        demand_proportion = 10;
        demand_options = 0;
        parking_duration_a_shape = 1;
        parking_duration_b_scale = 1;
        i = int16(supply_options/0.1) + 5;
        [average_searching_time_s(i,k,supply_proportion + 1), average_searching_distance_s(i,k,supply_proportion + 1), final_cum_revenue_s(i,k,supply_proportion + 1), ...
            ~, ~, ~, ~, optimal_target_occupancy_rate_all_vehicles_s(i,k,supply_proportion + 1), ...
            optimal_target_occupancy_rate_fuel_vehicles_s(i,k,supply_proportion + 1), ...
            optimal_target_occupancy_rate_electric_vehicles_s(i,k,supply_proportion + 1),...
            social_impacts_s(i,k,supply_proportion + 1)] = ...
            Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options, supply_options + 1, parking_duration_a_shape, parking_duration_b_scale);

        close all

    end
end
clearvars demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average parking duration in the area (k = 3)
k = 3; %flexible demand
% sensitivity_parameter = k;
for demand_proportion = 0:1:maximum_proportion_value
    for parking_duration_a_shape = -0.4:0.1:0.4

    clearvars -except average_searching_time_d average_searching_time_s average_searching_distance_d average_searching_distance_s final_cum_revenue_d final_cum_revenue_s optimal_target_occupancy_rate_all_vehicles_d optimal_target_occupancy_rate_all_vehicles_s optimal_target_occupancy_rate_fuel_vehicles_d optimal_target_occupancy_rate_fuel_vehicles_s optimal_target_occupancy_rate_electric_vehicles_d optimal_target_occupancy_rate_electric_vehicles_s social_impacts_d social_impacts_s policy demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale i k maximum_value maximum_proportion_value

        supply_proportion = 10;
        demand_options = 0;
        supply_options = 0;
        parking_duration_b_scale = 1;
        i = int16(parking_duration_a_shape/0.1) + 5;
        [average_searching_time_d(i,k,demand_proportion + 1), average_searching_distance_d(i,k,demand_proportion + 1), final_cum_revenue_d(i,k,demand_proportion + 1), ...
            ~, ~, ~, ~, optimal_target_occupancy_rate_all_vehicles_d(i,k,demand_proportion + 1), ...
            optimal_target_occupancy_rate_fuel_vehicles_d(i,k,demand_proportion + 1), ...
            optimal_target_occupancy_rate_electric_vehicles_d(i,k,demand_proportion + 1),...
            social_impacts_d(i,k,demand_proportion + 1)] = ...
            Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options, supply_options, parking_duration_a_shape + 1, parking_duration_b_scale);

        close all

    end
end
clearvars demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average parking duration in the area (k = 3)
k = 3; %flexible supply
% sensitivity_parameter = k;
for supply_proportion = 0:1:maximum_proportion_value
    for parking_duration_a_shape = -0.4:0.1:0.4

    clearvars -except average_searching_time_d average_searching_time_s average_searching_distance_d average_searching_distance_s final_cum_revenue_d final_cum_revenue_s optimal_target_occupancy_rate_all_vehicles_d optimal_target_occupancy_rate_all_vehicles_s optimal_target_occupancy_rate_fuel_vehicles_d optimal_target_occupancy_rate_fuel_vehicles_s optimal_target_occupancy_rate_electric_vehicles_d optimal_target_occupancy_rate_electric_vehicles_s social_impacts_d social_impacts_s policy demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale i k maximum_value maximum_proportion_value

        demand_proportion = 10;
        demand_options = 0;
        supply_options = 0;
        parking_duration_b_scale = 1;
        i = int16(parking_duration_a_shape/0.1) + 5;
        [average_searching_time_s(i,k,supply_proportion + 1), average_searching_distance_s(i,k,supply_proportion + 1), final_cum_revenue_s(i,k,supply_proportion + 1), ...
            ~, ~, ~, ~, optimal_target_occupancy_rate_all_vehicles_s(i,k,supply_proportion + 1), ...
            optimal_target_occupancy_rate_fuel_vehicles_s(i,k,supply_proportion + 1), ...
            optimal_target_occupancy_rate_electric_vehicles_s(i,k,supply_proportion + 1),...
            social_impacts_s(i,k,supply_proportion + 1)] = ...
            Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options, supply_options, parking_duration_a_shape + 1, parking_duration_b_scale);

        close all

    end
end

% save ('Sensitivity_analysis_contour_policy_1.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1a.) Optimal target parking occupancy rate (actual rate, contour plot, demand proportion):
% demand entering the area (k = 1) and flexible demand
figure
for r = -0.4:0.1:0.4
    i = int16(r/0.1) + 5;
    for p = 0:1:maximum_proportion_value
        x_value(i,1) = r;
        y_value(p + 1,1) = p;
        z_value(i,p + 1) = optimal_target_occupancy_rate_all_vehicles_d(i,1,p + 1);
    end
end
z_value_mean = movmean(movmean(z_value(:,:),4,2),4,1);

newpoints = 500;
[xq,yq] = ndgrid(x_value, y_value);
F = griddedInterpolant(xq,yq,z_value_mean,'cubic');
[xi,yi] = ndgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints));
z_value_interp = F(xi,yi);

contourf(xi,yi,z_value_interp,50)
xlabel('% Change of demand entering the area')
ylabel(['Proportion of electric vehicles in demand, ', char(949)])
% axis([0,20,0,20])

xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1b.) Optimal target parking occupancy rate (actual rate, contour plot, supply proportion):
% demand entering the area (k = 1) and flexible supply
figure
clearvars x_value y_value z_value z_value_mean xq yq F xi yi z_value_interp
for r = -0.4:0.1:0.4
    i = int16(r/0.1) + 5;
    for p = 0:1:maximum_proportion_value
        x_value(i,1) = r;
        y_value(p + 1,1) = p;
        z_value(i,p + 1) = optimal_target_occupancy_rate_all_vehicles_s(i,1,p + 1);
    end
end
z_value_mean = movmean(movmean(z_value(:,:),4,2),4,1);

newpoints = 500;
[xq,yq] = ndgrid(x_value, y_value);
F = griddedInterpolant(xq,yq,z_value_mean,'cubic');
[xi,yi] = ndgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints));
z_value_interp = F(xi,yi);

contourf(xi,yi,z_value_interp,50)
xlabel('% Change of demand entering the area')
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])

xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2a.) Optimal target parking occupancy rate (actual rate, contour plot, demand proportion):
% Total number of existing parking spaces in the area (k = 2) and flexible demand
figure
clearvars x_value y_value z_value z_value_mean xq yq F xi yi z_value_interp
for r = -0.4:0.1:0.4
    i = int16(r/0.1) + 5;
    for p = 0:1:maximum_proportion_value
        x_value(i,1) = r;
        y_value(p + 1,1) = p;
        z_value(i,p + 1) = optimal_target_occupancy_rate_all_vehicles_d(i,2,p + 1);
    end
end
z_value_mean = movmean(movmean(z_value(:,:),4,2),4,1);

newpoints = 500;
[xq,yq] = ndgrid(x_value, y_value);
F = griddedInterpolant(xq,yq,z_value_mean,'cubic');
[xi,yi] = ndgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints));
z_value_interp = F(xi,yi);

contourf(xi,yi,z_value_interp,50)
xlabel('% Change of total number of existing parking spaces in the area')
ylabel(['Proportion of electric vehicles in demand, ', char(949)])
% axis([0,20,0,20])

xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2b.) Optimal target parking occupancy rate (actual rate, contour plot, supply proportion):
% Total number of existing parking spaces in the area (k = 2) and flexible supply
figure
clearvars x_value y_value z_value z_value_mean xq yq F xi yi z_value_interp
for r = -0.4:0.1:0.4
    i = int16(r/0.1) + 5;
    for p = 0:1:maximum_proportion_value
        x_value(i,1) = r;
        y_value(p + 1,1) = p;
        z_value(i,p + 1) = optimal_target_occupancy_rate_all_vehicles_s(i,2,p + 1);
    end
end
z_value_mean = movmean(movmean(z_value(:,:),4,2),4,1);

newpoints = 500;
[xq,yq] = ndgrid(x_value, y_value);
F = griddedInterpolant(xq,yq,z_value_mean,'cubic');
[xi,yi] = ndgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints));
z_value_interp = F(xi,yi);

contourf(xi,yi,z_value_interp,50)
xlabel('% Change of total number of existing parking spaces in the area')
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])

xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3a.) Optimal target parking occupancy rate (actual rate, contour plot, demand proportion):
% Average parking duration in the area (k = 3) and flexible demand
figure
clearvars x_value y_value z_value z_value_mean xq yq F xi yi z_value_interp
for r = -0.4:0.1:0.4
    i = int16(r/0.1) + 5;
    for p = 0:1:maximum_proportion_value
        x_value(i,1) = r;
        y_value(p + 1,1) = p;
        z_value(i,p + 1) = optimal_target_occupancy_rate_all_vehicles_d(i,3,p + 1);
    end
end
z_value_mean = movmean(movmean(z_value(:,:),4,2),4,1);

newpoints = 500;
[xq,yq] = ndgrid(x_value, y_value);
F = griddedInterpolant(xq,yq,z_value_mean,'cubic');
[xi,yi] = ndgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints));
z_value_interp = F(xi,yi);

contourf(xi,yi,z_value_interp,50)
xlabel('% Change of pdf of electric vehicles’ parking durations in the area')
ylabel(['Proportion of electric vehicles in demand, ', char(949)])
% axis([0,20,0,20])

xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3b.) Optimal target parking occupancy rate (actual rate, contour plot, supply proportion):
% Average parking duration in the area (k = 3) and flexible supply
figure
clearvars x_value y_value z_value z_value_mean xq yq F xi yi z_value_interp
for r = -0.4:0.1:0.4
    i = int16(r/0.1) + 5;
    for p = 0:1:maximum_proportion_value
        x_value(i,1) = r;
        y_value(p + 1,1) = p;
        z_value(i,p + 1) = optimal_target_occupancy_rate_all_vehicles_s(i,3,p + 1);
    end
end
z_value_mean = movmean(movmean(z_value(:,:),4,2),4,1);

newpoints = 500;
[xq,yq] = ndgrid(x_value, y_value);
F = griddedInterpolant(xq,yq,z_value_mean,'cubic');
[xi,yi] = ndgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints));
z_value_interp = F(xi,yi);

contourf(xi,yi,z_value_interp,50)
xlabel('% Change of pdf of electric vehicles’ parking durations in the area')
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])

xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%