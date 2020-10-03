% Run Marcopark ALL proportions of electric vehicles in demand and 
% parking spaces for electric vehicles in supply:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maximum_value = 50;
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
value_of_time = [    
     29.9/60
     25.4/60
     25.8/60
     17.2/60];
h1_setGlobal_EV_infrastructure(10000);
h1_setGlobal_VOT(value_of_time);


i = 0;
% for demand_proportion = maximum_value:1:100
%     for supply_proportion = maximum_value:1:100
for demand_proportion = 0:1:maximum_value
    for supply_proportion = 0:1:maximum_value
        
    i = i + 1;    
    clearvars -except average_searching_time average_searching_distance final_cum_revenue occupancy_fuel occupancy_electric occupancy_all ACT optimal_target_occupancy_rate_all_vehicles optimal_target_occupancy_rate_fuel_vehicles optimal_target_occupancy_rate_electric_vehicles social_impacts policy demand_proportion supply_proportion demand_options supply_options parking_duration_a_shape parking_duration_b_scale i ratio maximum_value
    
    ratio(i,1) = demand_proportion/supply_proportion;
    [average_searching_time(i,1), average_searching_distance(i,1), final_cum_revenue(i,1), ...
        occupancy_fuel(i,:), occupancy_electric(i,:), occupancy_all(i,:), ACT(i,:),...
        optimal_target_occupancy_rate_all_vehicles(i,1), ...
        optimal_target_occupancy_rate_fuel_vehicles(i,1), ...
        optimal_target_occupancy_rate_electric_vehicles(i,1), social_impacts(i,1)] = ...
        Marcopark(policy, demand_proportion/100, supply_proportion/100, demand_options, supply_options, parking_duration_a_shape, parking_duration_b_scale);
    
    if average_searching_time(i,1) == 0 && average_searching_distance(i,1) == 0 && final_cum_revenue(i,1) == 0 && occupancy_fuel(i,1) == 0 && occupancy_electric(i,1) == 0 && occupancy_all(i,1) == 0 && ACT(i,1) == 0 && optimal_target_occupancy_rate_all_vehicles(i,1) == 0 && optimal_target_occupancy_rate_fuel_vehicles(i,1) == 0 && optimal_target_occupancy_rate_electric_vehicles(i,1) == 0
        ratio(i,1) = ratio(i-1,1);
        average_searching_time(i,1) = average_searching_time(i-1,1);
        average_searching_distance(i,1) = average_searching_distance(i-1,1);
        final_cum_revenue(i,1) = final_cum_revenue(i-1,1);
        occupancy_fuel(i,:) = occupancy_fuel(i-1,:);
        occupancy_electric(i,:) = occupancy_electric(i-1,:);
        occupancy_all(i,:) = occupancy_all(i-1,:);
        ACT(i,:) = ACT(i-1,:);
        optimal_target_occupancy_rate_all_vehicles(i,1) = optimal_target_occupancy_rate_all_vehicles(i-1,1);
        optimal_target_occupancy_rate_fuel_vehicles(i,1) = optimal_target_occupancy_rate_fuel_vehicles(i-1,1);
        optimal_target_occupancy_rate_electric_vehicles(i,1) = optimal_target_occupancy_rate_electric_vehicles(i-1,1);
        social_impacts(i,1) = social_impacts(i-1,1);
    end
    demand_proportion
    supply_proportion

    close all
    end
end

% save('results_trade_offs_demand_and_supply_for_electric_vehicles_policy_1.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 1.) Average searching time:
% plot_start = 8;
% plot_end = 14;
% plot_line{1} = [':'];
% plot_line{2} = ['--'];
% plot_line{3} = ['-.'];
% plot_line{4} = ['-'];
% plot_line{5} = ['-.'];
% plot_line{6} = ['--'];
% plot_line{7} = [':'];
% plot_marker{1} = ['+'];
% plot_marker{2} = ['*'];
% plot_marker{3} = ['d'];
% plot_marker{4} = ['x'];
% plot_marker{5} = ['s'];
% plot_marker{6} = ['o'];
% plot_marker{7} = ['h'];
% % black blue green red magenta cyan orange
% cmap = [0 0 0; 0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 0.9100 0.4100 0.1700];
% % cmap = colormap(parula(5));
% % cmap = colormap(jet(5));
% for k = plot_start:plot_end
%     hold on
%     x_value = ratio(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
%     y_value = average_searching_time(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
% 
%     scatter(x_value, y_value, 75, cmap(k - (plot_start - 1),:), 'Marker', plot_marker{k - (plot_start - 1)}, 'LineWidth', 2)
%     hold on
%     leg{k - (plot_start - 1)} = [num2str(k-1), '% electric vehicles in demand'];
%     relevant_leg(k - (plot_start - 1)) = plot(x_value, y_value, 'Color', cmap(k - (plot_start - 1),:), 'LineStyle', plot_line{k - (plot_start - 1)}, 'LineWidth', 2, 'DisplayName', leg{k - (plot_start - 1)});
%     
% end
% xlabel('Proportion of electric vehicles in demand / Proportion of parking spaces for only electric vehicles')
% ylabel('Average searching time (in min/veh)')
% set(gca,'FontSize',17)
% hold off
% leg = legend(relevant_leg);
% leg.ItemTokenSize = [50,18];
% axis([0,13,0,75])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1a.) Average searching time - Contour plot:
figure
plot_max = 50;
for k = 1:1:(maximum_value + 1)
    x_value_all(k,1) = k - 1;
    y_value_all(k,1) = k - 1;
    z_value_all(:,k) = average_searching_time(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
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
xlabel(['Proportion of electric vehicles in demand, ', char(949)])
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])
set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [1 317])
caxis([1 317])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 2.) Average searching distance:
% clear leg
% figure
% for k = plot_start:plot_end
%     hold on
%     x_value = ratio(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
%     y_value = average_searching_distance(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
% 
%     scatter(x_value, y_value, 75, cmap(k - (plot_start - 1),:), 'Marker', plot_marker{k - (plot_start - 1)}, 'LineWidth', 2)
%     hold on
%     leg{k - (plot_start - 1)} = [num2str(k-1), '% electric vehicles in demand'];
%     relevant_leg(k - (plot_start - 1)) = plot(x_value, y_value, 'Color', cmap(k - (plot_start - 1),:), 'LineStyle', plot_line{k - (plot_start - 1)}, 'LineWidth', 2, 'DisplayName', leg{k - (plot_start - 1)});
%     
% end
% xlabel('Proportion of electric vehicles in demand / Proportion of parking spaces for only electric vehicles')
% ylabel('Average searching distance (in km/veh)')
% set(gca,'FontSize',17)
% hold off
% leg = legend(relevant_leg);
% leg.ItemTokenSize = [50,18];
% axis([0,13,0,16])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2a.) Average searching distance - Contour plot:
figure
plot_max = 50;
for k = 1:1:(maximum_value + 1)
    x_value_all(k,1) = k - 1;
    y_value_all(k,1) = k - 1;    
    z_value_all(:,k) = average_searching_distance(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
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
xlabel(['Proportion of electric vehicles in demand, ', char(949)])
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])
set(gca,'FontSize',21)
hold off
colormap(jet)
colorbar
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 3.) Total revenue from parking pricing:
% clear leg
% figure
% for k = plot_start:plot_end
%     hold on
%     x_value = ratio(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
%     y_value = final_cum_revenue(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
% 
%     scatter(x_value, y_value, 75, cmap(k - (plot_start - 1),:), 'Marker', plot_marker{k - (plot_start - 1)}, 'LineWidth', 2)
%     hold on
%     leg{k - (plot_start - 1)} = [num2str(k-1), '% electric vehicles in demand'];
%     relevant_leg(k - (plot_start - 1)) = plot(x_value, y_value, 'Color', cmap(k - (plot_start - 1),:), 'LineStyle', plot_line{k - (plot_start - 1)}, 'LineWidth', 2, 'DisplayName', leg{k - (plot_start - 1)});
%     
% end
% xlabel('Proportion of electric vehicles in demand / Proportion of parking spaces for only electric vehicles')
% ylabel('Total revenue from parking pricing (in CHF)')
% set(gca,'FontSize',17)
% hold off
% leg = legend(relevant_leg);
% leg.ItemTokenSize = [50,18];
% axis([0,13,15600,17700])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3a.) Total revenue from parking pricing - Contour plot:
figure
plot_max = 50;
for k = 1:1:(maximum_value + 1)
    x_value_all(k,1) = k - 1;
    y_value_all(k,1) = k - 1;    
    z_value_all(:,k) = final_cum_revenue(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
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
xlabel(['Proportion of electric vehicles in demand, ', char(949)])
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])
set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [8793 17630])
caxis([8793 17630])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 4.) Optimal target parking occupancy rate:
% clear leg
% figure
% for k = plot_start:plot_end
%     hold on
%     x_value = ratio(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
%     y_value = optimal_target_occupancy_rate_all_vehicles(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
% 
%     scatter(x_value, y_value, 75, cmap(k - (plot_start - 1),:), 'Marker', plot_marker{k - (plot_start - 1)}, 'LineWidth', 2)
%     hold on
%     leg{k - (plot_start - 1)} = [num2str(k-1), '% electric vehicles in demand'];
%     relevant_leg(k - (plot_start - 1)) = plot(x_value, y_value, 'Color', cmap(k - (plot_start - 1),:), 'LineStyle', plot_line{k - (plot_start - 1)}, 'LineWidth', 2, 'DisplayName', leg{k - (plot_start - 1)});
%     
% end
% xlabel('Proportion of electric vehicles in demand / Proportion of parking spaces for only electric vehicles')
% ylabel('Optimal target parking occupancy rate (all parking spaces)')
% set(gca,'FontSize',17)
% hold off
% ay = gca;
% ay.YTick = 0:0.1:1;
% yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
% percentsy = repmat('%', length(yticks),1);  %  equal to the size
% yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
% set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
% leg = legend(relevant_leg);
% leg.ItemTokenSize = [50,18];
% axis([0,13,0,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4a.) Optimal target parking occupancy rate - Contour plot:
figure
plot_max = 50;
for k = 1:1:(maximum_value + 1)
    x_value_all(k,1) = k - 1;
    y_value_all(k,1) = k - 1;    
    z_value_all(:,k) = optimal_target_occupancy_rate_all_vehicles(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
end
x_value = x_value_all(1:(plot_max + 1),1);
y_value = y_value_all(1:(plot_max + 1),1);
z_value = movmean(movmean(z_value_all(1:(plot_max + 1),1:(plot_max + 1)),4,2),4,1);

newpoints = 500;
[xq,yq] = meshgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints ),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints )...
          );
z_value_interp = interp2(x_value,y_value,z_value,xq,yq,'cubic');
for m = 1:size(z_value_interp,1)
    for n = 1:size(z_value_interp,2)
        if z_value_interp(m,n) < 0
            z_value_interp(m,n) = 0;
        end
        if z_value_interp(m,n) > 1
            z_value_interp(m,n) = 1;
        end
    end
end     
contourf(xq,yq,z_value_interp,50)
xlabel(['Proportion of electric vehicles in demand, ', char(949)])
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])
set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 5.) Optimal target parking occupancy rate (fuel vehicles):
% clear leg
% figure
% for k = plot_start:plot_end
%     hold on
%     x_value = ratio(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
%     y_value = optimal_target_occupancy_rate_fuel_vehicles(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
% 
%     scatter(x_value, y_value, 75, cmap(k - (plot_start - 1),:), 'Marker', plot_marker{k - (plot_start - 1)}, 'LineWidth', 2)
%     hold on
%     leg{k - (plot_start - 1)} = [num2str(k-1), '% electric vehicles in demand'];
%     relevant_leg(k - (plot_start - 1)) = plot(x_value, y_value, 'Color', cmap(k - (plot_start - 1),:), 'LineStyle', plot_line{k - (plot_start - 1)}, 'LineWidth', 2, 'DisplayName', leg{k - (plot_start - 1)});
%     
% end
% xlabel('Proportion of electric vehicles in demand / Proportion of parking spaces for only electric vehicles')
% ylabel('Optimal target parking occupancy rate (parking spaces for fuel vehicles)')
% set(gca,'FontSize',17)
% hold off
% ay = gca;
% ay.YTick = 0:0.1:1;
% yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
% percentsy = repmat('%', length(yticks),1);  %  equal to the size
% yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
% set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
% leg = legend(relevant_leg);
% leg.ItemTokenSize = [50,18];
% axis([0,13,0,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5a.) Optimal target parking occupancy rate (fuel vehicles) - Contour plot:
figure
plot_max = 50;
for k = 1:1:(maximum_value + 1)
    x_value_all(k,1) = k - 1;
    y_value_all(k,1) = k - 1;    
    z_value_all(:,k) = optimal_target_occupancy_rate_fuel_vehicles(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
end
x_value = x_value_all(1:(plot_max + 1),1);
y_value = y_value_all(1:(plot_max + 1),1);
% z_value = z_value_all(1:(plot_max + 1),1:(plot_max + 1));
z_value = movmean(movmean(z_value_all(1:(plot_max + 1),1:(plot_max + 1)),4,2),4,1);

newpoints = 500;
[xq,yq] = meshgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints ),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints )...
          );
z_value_interp = interp2(x_value,y_value,z_value,xq,yq,'cubic');
for m = 1:size(z_value_interp,1)
    for n = 1:size(z_value_interp,2)
        if z_value_interp(m,n) < 0
            z_value_interp(m,n) = 0;
        end
        if z_value_interp(m,n) > 1
            z_value_interp(m,n) = 1;
        end
    end
end  
contourf(xq,yq,z_value_interp,50)
% contourf(x_value,y_value,z_value,50)
xlabel(['Proportion of electric vehicles in demand, ', char(949)])
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])
set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 6.) Optimal target parking occupancy rate (electric vehicles):
% clear leg
% figure
% for k = plot_start:plot_end
%     hold on
%     x_value = ratio(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
%     y_value = optimal_target_occupancy_rate_electric_vehicles(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
% 
%     scatter(x_value, y_value, 75, cmap(k - (plot_start - 1),:), 'Marker', plot_marker{k - (plot_start - 1)}, 'LineWidth', 2)
%     hold on
%     leg{k - (plot_start - 1)} = [num2str(k-1), '% electric vehicles in demand'];
%     relevant_leg(k - (plot_start - 1)) = plot(x_value, y_value, 'Color', cmap(k - (plot_start - 1),:), 'LineStyle', plot_line{k - (plot_start - 1)}, 'LineWidth', 2, 'DisplayName', leg{k - (plot_start - 1)});
%     
% end
% xlabel('Proportion of electric vehicles in demand / Proportion of parking spaces for only electric vehicles')
% ylabel('Optimal target parking occupancy rate (parking spaces for electric vehicles)')
% set(gca,'FontSize',17)
% hold off
% ay = gca;
% ay.YTick = 0:0.1:1;
% yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
% percentsy = repmat('%', length(yticks),1);  %  equal to the size
% yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
% set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
% leg = legend(relevant_leg);
% leg.ItemTokenSize = [50,18];
% axis([0,13,0,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6a.) Optimal target parking occupancy rate (electric vehicles) - Contour plot:
figure
plot_max = 50;
for k = 1:1:(maximum_value + 1)
    x_value_all(k,1) = k - 1;
    y_value_all(k,1) = k - 1;    
    z_value_all(:,k) = optimal_target_occupancy_rate_electric_vehicles(((maximum_value + 1)*(k-1)+1):1:((maximum_value + 1)*k),1);
end
x_value = x_value_all(2:(plot_max + 1),1);
y_value = y_value_all(2:(plot_max + 1),1);
% z_value = z_value_all(2:(plot_max + 1),2:(plot_max + 1));
z_value = movmean(movmean(z_value_all(2:(plot_max + 1),2:(plot_max + 1)),4,2),4,1);

newpoints = 500;
[xq,yq] = meshgrid(...
            linspace(min(min(x_value,[],2)),max(max(x_value,[],2)),newpoints ),...
            linspace(min(min(y_value,[],1)),max(max(y_value,[],1)),newpoints )...
          );
z_value_interp = interp2(x_value,y_value,z_value,xq,yq,'cubic');
for m = 1:size(z_value_interp,1)
    for n = 1:size(z_value_interp,2)
        if z_value_interp(m,n) < 0
            z_value_interp(m,n) = 0;
        end
        if z_value_interp(m,n) > 1
            z_value_interp(m,n) = 1;
        end
    end
end  
contourf(xq,yq,z_value_interp,50)
% contourf(x_value,y_value,z_value,50)
xlabel(['Proportion of electric vehicles in demand, ', char(949)])
ylabel(['Proportion of parking spaces for only electric vehicles, ', char(950)])
% axis([0,20,0,20])
set(gca,'FontSize',21)
hold off
colormap(jet)
cb = colorbar();
set(cb, 'Limits', [0 1], 'Ticks', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], 'TickLabels', {'10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'})
caxis([0 1])
grid
