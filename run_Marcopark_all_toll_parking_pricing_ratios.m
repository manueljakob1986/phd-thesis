% Run Marcopark ALL congestion toll and parking pricing ratios:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maximum_value = 60;

average_searching_time(:,:) = zeros(maximum_value,1);
average_searching_distance(:,:) = zeros(maximum_value,1); 
final_cum_revenue(:,:) = zeros(maximum_value,1); 
cum_revenue_PR(:,:) = zeros(maximum_value,1); 
cum_revenue_congestion_toll(:,:) = zeros(maximum_value,1); 
cum_revenue_parking_pricing(:,:) = zeros(maximum_value,1);
share_veh_pr_decision_mean_all_times(:,:) = zeros(maximum_value,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes the parking and P+R capacity real-time
% information, switch it on or off.
% parking and P+R capacity real-time information switched on = 1
% parking and P+R capacity real-time information switched off = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1_setGlobal_parking_pr_capacity_information(0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for all ratios between 0.5 and maximum_value (here 30)
for j = 0.5:0.5:(maximum_value/2)
    
    % set global variable initial parking pricing:
    h1_setGlobal_initial_parking_pricing(c16_input_toll/j);
          
    clearvars -except average_searching_time average_searching_distance final_cum_revenue cum_revenue_PR cum_revenue_congestion_toll cum_revenue_parking_pricing share_veh_pr_decision_mean_all_times i j maximum_value

    i = j/0.5;
    [average_searching_time(i,1), average_searching_distance(i,1), final_cum_revenue(i,1), cum_revenue_PR(i,1), cum_revenue_congestion_toll(i,1), ...
       cum_revenue_parking_pricing(i,1), share_veh_pr_decision_mean_all_times(i,1)] = Marcopark(0,0);

    close all

end

% save ('Finding_optimal_congestion_toll_parking_price_9_3DMFD.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.) Average searching time:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = average_searching_time(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Congestion toll / Parking price')
ylabel('Average searching time (in min/veh)')
set(gca,'FontSize',24)
% axis([0,maximum_value,5,18])
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1a.) Average searching time (Comparison - High congestion toll):
% figure('units','normalized','outerposition',[0 0 1 1])
% hold on
% x_value = 0.5:0.5:(maximum_value/2);
% plot(x_value, average_searching_time_16, 'k--', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_18, 'b:', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_20, 'r-', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_22, 'm-.', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_24, 'g--', 'LineWidth', 4)
% hold off
% xlabel('Congestion toll / Parking price')
% ylabel('Average searching time (in min/veh)')
% leg = legend('Congestion toll = 16 CHF', 'Congestion toll = 18 CHF', 'Congestion toll = 20 CHF', 'Congestion toll = 22 CHF', 'Congestion toll = 24 CHF');
% leg.ItemTokenSize = [50,18];
% set(gca,'FontSize',24)
% % axis([0,maximum_value,5,18])
% hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1b.) Average searching time (Comparison - Low congestion toll):
% figure('units','normalized','outerposition',[0 0 1 1])
% hold on
% x_value = 0.5:0.5:(maximum_value/2);
% plot(x_value, average_searching_time_9, 'k--', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_10, 'b:', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_11, 'c-.', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_12, 'r-', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_13, 'm-.', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_14, 'g:', 'LineWidth', 4)
% hold on 
% plot(x_value, average_searching_time_15, 'b--', 'LineWidth', 4)
% hold off
% xlabel('Congestion toll / Parking price')
% ylabel('Average searching time (in min/veh)')
% leg = legend('Congestion toll = 9 CHF', 'Congestion toll = 10 CHF', 'Congestion toll = 11 CHF', 'Congestion toll = 12 CHF', 'Congestion toll = 13 CHF', 'Congestion toll = 14 CHF', 'Congestion toll = 15 CHF');
% leg.ItemTokenSize = [50,18];
% set(gca,'FontSize',24)
% % axis([0,maximum_value,5,18])
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.) Average searching distance:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = average_searching_distance(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Congestion toll / Parking price')
ylabel('Average searching distance (in km/veh)')
set(gca,'FontSize',24)
% axis([0,maximum_value,5,18])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3.) Total revenue from P+R and from congestion tolls and/or parking pricing:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = final_cum_revenue(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Congestion toll / Parking price')
ylabel('Total revenue (in CHF)')
set(gca,'FontSize',24)
% axis([0,maximum_value,5,18])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.) Total revenue from P+R and PT:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = cum_revenue_PR(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Congestion toll / Parking price')
ylabel('Total revenue from P+R facilities and PT (in CHF)')
set(gca,'FontSize',24)
% axis([0,maximum_value,5,18])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5.) Total revenue congestion tolls:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = cum_revenue_congestion_toll(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Congestion toll / Parking price')
ylabel('Total revenue from congestion tolls (in CHF)')
set(gca,'FontSize',24)
% axis([0,maximum_value,5,18])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6.) Total revenue from parking pricing:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = cum_revenue_parking_pricing(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Congestion toll / Parking price')
ylabel('Total revenue from parking pricing (in CHF)')
set(gca,'FontSize',24)
% axis([0,maximum_value,5,18])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7.) Choice of drivers for entering the network by vehicle:
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = share_veh_pr_decision_mean_all_times(:,1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Congestion toll / Parking price')
ylabel('Choice of drivers for entering the network by vehicle')
set(gca,'FontSize',24)
% axis([0,maximum_value,5,18])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 8.) All revenues:
figure('units','normalized','outerposition',[0 0 1 1])

hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = final_cum_revenue(:,1);
plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = cum_revenue_PR(:,1);
plot(x_value, y_value, 'r--', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = cum_revenue_congestion_toll(:,1);
plot(x_value, y_value, 'g:', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = 0.5:0.5:(maximum_value/2);
y_value = cum_revenue_parking_pricing(:,1);
plot(x_value, y_value, 'm-.', 'MarkerSize', 10, 'LineWidth', 4)

xlabel('Congestion toll / Parking price')
ylabel('Revenue (in CHF)')
leg = legend('Total revenue','Total revenue from P+R facilities and PT','Total revenue from congestion tolls', 'Total revenue from parking pricing');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
axis([0,30,0,inf])
hold off
