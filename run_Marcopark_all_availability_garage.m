% Run Marcopark ALL Availability Garage Usage:
% Multiple availability garage usage parameter scenarios vs. avg searching time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % This input parameter describes the parking garage capacity
    % information, switch it on or off.
    % parking garage capacity information switched on = 1
    % parking garage capacity information switched off = 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h1_setGlobal_initial_parking_pricing(1.5);
    h1_setGlobal_initial_garage_parking_pricing(3);
    h1_setGlobal_number_of_parking_garages(2);
    h1_setGlobal_capacity_garage(166);
    h1_setGlobal_number_of_on_street_parking_spaces(207);
    h1_setGlobal_garage_usage_information(1);
    h1_setGlobal_s_dgp_penalty(2);
            
average_searching_time(:,1) = zeros(21,1); 
total_searching_time(:,1) = zeros(21,1); 
average_non_searching_time_parkers(:,1) = zeros(21,1); 
total_non_searching_time_parkers(:,1) = zeros(21,1); 
average_searching_distance(:,1) = zeros(21,1); 
total_searching_distance(:,1) = zeros(21,1); 
average_non_searching_distance_parkers(:,1) = zeros(21,1); 
total_non_searching_distance_parkers(:,1) = zeros(21,1);
average_deciding_gp_time(:,1) = zeros(21,1); 
total_deciding_gp_time(:,1) = zeros(21,1); 
average_deciding_gp_distance(:,1) = zeros(21,1); 
total_deciding_gp_distance(:,1) = zeros(21,1);
avg_total_time(:,1) = zeros(21,1); 
tot_time(:,1) = zeros(21,1); 
avg_total_distance(:,1) = zeros(21,1); 
tot_distance(:,1) = zeros(21,1); 
final_cum_revenue(:,1) = zeros(21,1); 
final_cum_revenue_garage(:,1) = zeros(21,1); 
total_revenue(:,1) = zeros(21,1);
total_on_street_revenue(:,1) = zeros(21,1);
total_garage_revenue(:,1) = zeros(21,1);
enterthearea_1(:,:,:) = zeros(1440,1,21);
enterthearea_2(:,:,:) = zeros(1440,1,21);
enterthearea_3(:,:,:) = zeros(1440,1,21);
enterthearea(:,:,:) = zeros(1440,1,21); 
n_vehicles_ns_dgp_VOT(:,:,:) = zeros(1440,4,21);
starttosearch_2(:,:,:) = zeros(1440,1,21);
starttosearch_3(:,:,:) = zeros(1440,1,21);
starttosearch(:,:,:) = zeros(1440,1,21);
n_vehicles_s_dgp_VOT(:,:,:) = zeros(1440,4,21);
findparking_2(:,:,:) = zeros(1440,1,21);
findparking_3(:,:,:) = zeros(1440,1,21);
findparking(:,:,:) = zeros(1440,1,21);
n_enter_garage(:,:,:) = zeros(1440,4,21);
n_dgp_searching(:,:,:) = zeros(1440,4,21); 
departparking_2(:,:,:) = zeros(1440,1,21); 
departparking_3(:,:,:) = zeros(1440,1,21); 
departparking(:,:,:) = zeros(1440,1,21); 
ndepart_garage(:,:,:) = zeros(1440,1,21); 
leavethearea_1(:,:,:) = zeros(1440,1,21); 
leavethearea_2(:,:,:) = zeros(1440,1,21); 
leavethearea_3(:,:,:) = zeros(1440,1,21); 
leavethearea(:,:,:) = zeros(1440,1,21);

% set global variable initial on-street parking pricing:
for j = 0:0.05:1
       
    h1_setGlobal_garage_usage_information_parameter(j);
    
    clearvars -except global_initial_garage_parking_pricing k i j m average_searching_time total_searching_time average_non_searching_time_parkers total_non_searching_time_parkers...
    average_searching_distance total_searching_distance average_non_searching_distance_parkers total_non_searching_distance_parkers ...
    average_deciding_gp_time total_deciding_gp_time average_deciding_gp_distance total_deciding_gp_distance ...
    avg_total_time tot_time avg_total_distance tot_distance final_cum_revenue final_cum_revenue_garage total_revenue total_on_street_revenue total_garage_revenue ...
    enterthearea_1 enterthearea_2 enterthearea_3 enterthearea n_vehicles_ns_dgp_VOT starttosearch_2 starttosearch_3 starttosearch ...
    n_vehicles_s_dgp_VOT findparking_2 findparking_3 findparking n_enter_garage n_dgp_searching ...
    departparking_2 departparking_3 departparking ndepart_garage ...
    leavethearea_1 leavethearea_2 leavethearea_3 leavethearea set_scenario

    i = int8(j/0.05 + 1);
        
    [average_searching_time(i,1), total_searching_time(i,1), average_non_searching_time_parkers(i,1), total_non_searching_time_parkers(i,1),...
    average_searching_distance(i,1), total_searching_distance(i,1), average_non_searching_distance_parkers(i,1), total_non_searching_distance_parkers(i,1), ...
    average_deciding_gp_time(i,1), total_deciding_gp_time(i,1), average_deciding_gp_distance(i,1), total_deciding_gp_distance(i,1),...    
    avg_total_time(i,1), tot_time(i,1), avg_total_distance(i,1), tot_distance(i,1), final_cum_revenue(i,1), final_cum_revenue_garage(i,1), total_revenue(i,1), total_on_street_revenue(i,1), total_garage_revenue(i,1), ...
    enterthearea_1(:,:,i), enterthearea_2(:,:,i), enterthearea_3(:,:,i), enterthearea(:,:,i), n_vehicles_ns_dgp_VOT(:,:,i), starttosearch_2(:,:,i), starttosearch_3(:,:,i), starttosearch(:,:,i),...
    n_vehicles_s_dgp_VOT(:,:,i), findparking_2(:,:,i), findparking_3(:,:,i), findparking(:,:,i), n_enter_garage(:,:,i), n_dgp_searching(:,:,i),...
    departparking_2(:,:,i), departparking_3(:,:,i), departparking(:,:,i), ndepart_garage(:,:,i), ...
    leavethearea_1(:,:,i), leavethearea_2(:,:,i), leavethearea_3(:,:,i), leavethearea(:,:,i)] = Marcopark();
        
    i
    
end

save('results_scenario_garage_usage_information.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1a.) Average searching time:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = average_searching_time(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Average searching time (in min/veh)')
set(gca,'FontSize',24)

axis([0 1 1.7 5.3])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Avg_searching_time_garage_usage_information.jpg')
saveas(gcf,'Avg_searching_time_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b.) Average dgp time:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = average_deciding_gp_time(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Average dgp time (in min/veh)')
set(gca,'FontSize',24)

axis([0 1 1.7 5.3])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Avg_dgp_time_garage_usage_information.jpg')
saveas(gcf,'Avg_dgp_time_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1c.) Average time:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = average_searching_time(:,1);
plot(x_value, y_value, 'b--', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value2 = 0:0.05:1;
y_value2 = average_deciding_gp_time(:,1);
plot(x_value2, y_value2, 'r-', 'MarkerSize', 10, 'LineWidth', 4)

legend('Vehicles searching for on-street parking','Vehicles deciding for garage parking')
xlabel('Influence parameter for garage usage information')
ylabel('Average time (in min/veh)')
set(gca,'FontSize',24)

axis([0 1 0 3.5])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1d.) Average searching and dgp time:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = average_searching_time(:,1) + average_deciding_gp_time(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Average searching and dgp time (in min/veh)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Avg_searching_and_dgp_time_garage_usage_information.jpg')
saveas(gcf,'Avg_searching_and_dgp_time_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.) Total searching time:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = total_searching_time(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Total searching time (in min/veh)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Total_searching_time_garage_usage_information.jpg')
saveas(gcf,'Total_searching_time_garage_usage_information.fig')

hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3a.) Average searching distance:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = average_searching_distance(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Average searching distance (in km/veh)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Avg_searching_distance_garage_usage_information.jpg')
saveas(gcf,'Avg_searching_distance_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3b.) Average dgp distance:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = average_deciding_gp_distance(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Average dgp distance (in km/veh)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Avg_dgp_distance_garage_usage_information.jpg')
saveas(gcf,'Avg_dgp_distance_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3c.) Average searching and dgp distance:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = average_searching_distance(:,1) + average_deciding_gp_distance(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Average searching and dgp distance (in km/veh)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Avg_searching_and_dgp_distance_garage_usage_information.jpg')
saveas(gcf,'Avg_searching_and_dgp_distance_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.) Total searching distance:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = total_searching_distance(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Total searching distance (in km/veh)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Total_searching_distance_garage_usage_information.jpg')
saveas(gcf,'Total_searching_distance_garage_usage_information.fig')

hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5.) Total on-street parking revenue:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = total_on_street_revenue(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Total on-street parking revenue (in CHF)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Total_on_street_revenue_garage_usage_information.jpg')
saveas(gcf,'Total_on_street_revenue_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6.) Total garage parking revenue:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = total_garage_revenue(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Total garage parking revenue (in CHF)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Total_garage_revenue_garage_usage_information.jpg')
saveas(gcf,'Total_garage_revenue_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7.) Total revenue:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:1;
y_value = total_revenue(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Influence parameter for garage usage information')
ylabel('Total revenue (in CHF)')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Total_revenue_both_garage_usage_information.jpg')
saveas(gcf,'Total_revenue_both_garage_usage_information.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 8.) ALL revenues:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

x_value = 0:0.05:1;
y_value = total_on_street_revenue(:,1);
plot(x_value, y_value, 'b--', 'MarkerSize', 10, 'LineWidth', 4)
hold on

x_value = 0:0.05:1;
y_value = total_garage_revenue(:,1);
plot(x_value, y_value, 'g:', 'MarkerSize', 10, 'LineWidth', 4)
hold on

x_value = 0:0.05:1;
y_value = total_revenue(:,1);
plot(x_value, y_value, 'r-', 'MarkerSize', 10, 'LineWidth', 4)
hold on

xlabel('Influence parameter for garage usage information')
ylabel('Total revenue (in CHF)')
legend('Total on-street parking revenue','Total garage parking revenue','Total revenue created by both on-street and garage parking')
set(gca,'FontSize',24)
% t = get(gca, 'ytick');
% set(gca, 'yticklabel', num2str(t'))

axis([0 1 0 17100])
ax = gca;
ax.XTick = 0:0.1:1;
set(gca,'FontSize',24)

saveas(gcf,'Total_revenue_garage_usage_information.jpg')
saveas(gcf,'Total_revenue_garage_usage_information.fig')

hold off

