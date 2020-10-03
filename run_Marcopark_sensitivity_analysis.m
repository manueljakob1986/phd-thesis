% Run Marcopark for sensitivity analysis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maximum_value = 11;
% maximum_value = 19;
number_of_sensitivity_analysis = 6;

average_searching_time(:,:) = zeros(maximum_value,number_of_sensitivity_analysis);
average_searching_distance(:,:) = zeros(maximum_value,number_of_sensitivity_analysis); 
final_cum_revenue(:,:) = zeros(maximum_value,number_of_sensitivity_analysis); 
cum_revenue_PR(:,:) = zeros(maximum_value,number_of_sensitivity_analysis); 
cum_revenue_congestion_toll(:,:) = zeros(maximum_value,number_of_sensitivity_analysis); 
cum_revenue_parking_pricing(:,:) = zeros(maximum_value,number_of_sensitivity_analysis);
share_veh_pr_decision_mean_all_times(:,:) = zeros(maximum_value,number_of_sensitivity_analysis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set global variable initial parking pricing to 4.5 CHF.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1_setGlobal_initial_parking_pricing(4.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes the parking and P+R capacity real-time
% information, switch it on or off.
% parking and P+R capacity real-time information switched on = 1
% parking and P+R capacity real-time information switched off = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1_setGlobal_parking_pr_capacity_information(0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity parameter: k
% 0 = no parameter used
% 1 = VOT
% 2 = number of existing parking spaces inside the area
% 3 = number of existing P+R spaces outside the area
% 4 = hourly P+R fee
% 5 = round-trip PT fee
% 6 = number of existing parking spaces inside the area (A and P dependent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:number_of_sensitivity_analysis
%     sensitivity_parameter = k;

    for j = -0.5:0.1:0.5

        clearvars -except average_searching_time average_searching_distance final_cum_revenue cum_revenue_PR cum_revenue_congestion_toll cum_revenue_parking_pricing share_veh_pr_decision_mean_all_times i j k number_of_sensitivity_analysis maximum_value

        i = int16(j/0.1) + 6;
        [average_searching_time(i,k), average_searching_distance(i,k), final_cum_revenue(i,k), cum_revenue_PR(i,k), cum_revenue_congestion_toll(i,k), ...
           cum_revenue_parking_pricing(i,k), share_veh_pr_decision_mean_all_times(i,k)] = Marcopark(k, j);

        close all

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6 = number of existing parking spaces inside the area (A and P dependent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k = 6:number_of_sensitivity_analysis
% %     sensitivity_parameter = k;
% 
%     for j = 0:0.05:0.9
% 
%         clearvars -except average_searching_time average_searching_distance final_cum_revenue cum_revenue_PR cum_revenue_congestion_toll cum_revenue_parking_pricing share_veh_pr_decision_mean_all_times i j k number_of_sensitivity_analysis maximum_value
% 
%         i = int16(j/0.05) + 1;
%         [average_searching_time(i,k), average_searching_distance(i,k), final_cum_revenue(i,k), cum_revenue_PR(i,k), cum_revenue_congestion_toll(i,k), ...
%            cum_revenue_parking_pricing(i,k), share_veh_pr_decision_mean_all_times(i,k)] = Marcopark(k, j);
% 
%         close all
% 
%     end
% end

save ('Sensitivity_analysis_A_P_independent_peak_3DMFD.mat')
% save ('Sensitivity_analysis_A_P_dependent_peak_3DMFD.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.) Average searching time:

figure('units','normalized','outerposition',[0 0 1 1])
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = average_searching_time(:,1)/average_searching_time(6,1) - 1;
% plot(x_value, y_value, 'k:', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = average_searching_time(:,2)/average_searching_time(6,2) - 1;
% plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = average_searching_time(:,3)/average_searching_time(6,3) - 1;
% plot(x_value, y_value, 'r--', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = average_searching_time(:,4)/average_searching_time(6,4) - 1;
% plot(x_value, y_value, 'g-.', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = average_searching_time(:,1);
plot(x_value, y_value, 'k:', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = average_searching_time(:,2);
plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = average_searching_time(:,3);
plot(x_value, y_value, 'r--', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = average_searching_time(:,4);
plot(x_value, y_value, 'g-.', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)
symlog('y')

xlabel('% Change of influencing factor')
% ylabel('% Change of average searching time')
ylabel('Average searching time (in min)')
leg = legend('VOT for all user groups','Total number of existing parking spaces inside the area','Total number of existing P+R spaces outside the area','Hourly P+R and round-trip PT fee');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
% axis([-0.5,0.5,-0.3,50])
hold off

ax.XTick = -0.5:0.1:0.5;
% ax.YTick = -0.3:0.1:0.3;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
% yticks = [-2; -1; -0.3; 0; 0.3; 1; 2];
percentsx = repmat('%', length(xticks),1);  %  equal to the size
% percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% yticklabel = [num2str(yticks)]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.) Average searching distance:

figure('units','normalized','outerposition',[0 0 1 1])
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = average_searching_distance(:,1)/average_searching_distance(6,1) - 1;
% plot(x_value, y_value, 'k:', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = average_searching_distance(:,2)/average_searching_distance(6,2) - 1;
% plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = average_searching_distance(:,3)/average_searching_distance(6,3) - 1;
% plot(x_value, y_value, 'r--', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = average_searching_distance(:,4)/average_searching_distance(6,4) - 1;
% plot(x_value, y_value, 'g-.', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = average_searching_distance(:,1);
plot(x_value, y_value, 'k:', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = average_searching_distance(:,2);
plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = average_searching_distance(:,3);
plot(x_value, y_value, 'r--', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = average_searching_distance(:,4);
plot(x_value, y_value, 'g-.', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)
symlog('y')

xlabel('% Change of influencing factor')
% ylabel('% Change of average searching distance')
ylabel('Average searching distance (in km)')
leg = legend('VOT for all user groups','Total number of existing parking spaces inside the area','Total number of existing P+R spaces outside the area','Hourly P+R and round-trip PT fee');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
% axis([-0.5,0.5,-1,45])
hold off

ax.XTick = -0.5:0.1:0.5;
% ax.YTick = -0.3:0.1:0.3;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
% percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3.) Total revenue from P+R and from congestion tolls and/or parking pricing:

figure('units','normalized','outerposition',[0 0 1 1])
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = final_cum_revenue(:,1)/final_cum_revenue(6,1) - 1;
% plot(x_value, y_value, 'k:', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = final_cum_revenue(:,2)/final_cum_revenue(6,2) - 1;
% plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = final_cum_revenue(:,3)/final_cum_revenue(6,3) - 1;
% plot(x_value, y_value, 'r--', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 4)
% 
% hold on
% x_value = -0.5:0.1:0.5;
% y_value = final_cum_revenue(:,4)/final_cum_revenue(6,4) - 1;
% plot(x_value, y_value, 'g-.', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = final_cum_revenue(:,1);
plot(x_value, y_value, 'k:', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = final_cum_revenue(:,2);
plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = final_cum_revenue(:,3);
plot(x_value, y_value, 'r--', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = -0.5:0.1:0.5;
y_value = final_cum_revenue(:,4);
plot(x_value, y_value, 'g-.', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)
% symlog('y')

xlabel('% Change of influencing factor')
ylabel('Total revenue (in CHF)')
leg = legend('VOT for all user groups','Total number of existing parking spaces inside the area','Total number of existing P+R spaces outside the area','Hourly P+R and round-trip PT fee');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%,d';
% axis([-0.5,0.5,0,32000])
hold off

ax.XTick = -0.5:0.1:0.5;
% ax.YTick = -0.3:0.1:0.3;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
% percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.) Average searching time (A und P are dependent):

figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:0.5;
y_value = average_searching_time(:,6)/average_searching_time(1,6) - 1;
plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

xlabel('% of A converted to P')
ylabel('% Change of average searching time')
set(gca,'FontSize',24)
axis([0,0.5,-0.35,0])
hold off

ax.XTick = 0:0.05:0.5;
ax.YTick = -0.35:-0.05:0;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5.) Total revenue from P+R and from congestion tolls and/or parking pricing (A und P are dependent):

figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:0.5;
y_value = final_cum_revenue(:,6)/final_cum_revenue(1,6) - 1;
plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

xlabel('% of A converted to P')
ylabel('% Change of total revenue')
set(gca,'FontSize',24)
axis([0,0.5,0,0.22])
hold off

ax.XTick = 0:0.1:0.5;
ax.YTick = -0.3:0.1:0.3;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6.) Average searching time and total revenue from P+R and from congestion tolls and/or parking pricing (A und P are dependent):

figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0:0.05:0.9;
y_value = average_searching_time(:,6)/average_searching_time(1,6) - 1;
plot(x_value, y_value, 'b-', 'Marker', '+', 'MarkerSize', 10, 'LineWidth', 4)

hold on
x_value = 0:0.05:0.9;
y_value = final_cum_revenue(:,6)/final_cum_revenue(1,6) - 1;
plot(x_value, y_value, 'r--', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 4)

xlabel('% of A converted to P')
ylabel('% Change')
leg = legend('Average searching time','Total revenue');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
% axis([0,0.5,-0.35,0])
hold off

ax.XTick = 0:0.05:0.9;
ax.YTick = -0.35:-0.05:0;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
percentsy = repmat('%', length(yticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
