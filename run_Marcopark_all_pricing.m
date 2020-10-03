% Run Marcopark ALL Pricing:
% Multiple on-street and garage parking prices vs. avg searching time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for scenario = 1:1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % This input parameter describes the parking garage capacity
    % information, switch it on or off.
    % parking garage capacity information switched on = 1
    % parking garage capacity information switched off = 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if scenario == 1
        h1_setGlobal_number_of_parking_garages(2)
        h1_setGlobal_capacity_garage(166)
        h1_setGlobal_number_of_on_street_parking_spaces(207)
        h1_setGlobal_garage_usage_information(0)
        h1_setGlobal_garage_usage_information_parameter(0.5)
        h1_setGlobal_s_dgp_penalty(2);
        set_scenario = 1;
        
%     elseif scenario == 2
%         h1_setGlobal_number_of_parking_garages(7)
%         h1_setGlobal_capacity_garage(296)
%         h1_setGlobal_number_of_on_street_parking_spaces(207)
%         h1_setGlobal_garage_usage_information(0)
%         h1_setGlobal_garage_usage_information_parameter(0.5)
%         h1_setGlobal_s_dgp_penalty(2);
%         set_scenario = 2;
%         
%     elseif scenario == 3
%         h1_setGlobal_number_of_parking_garages(2)
%         h1_setGlobal_capacity_garage(166)
%         h1_setGlobal_number_of_on_street_parking_spaces(207)
%         h1_setGlobal_garage_usage_information(1)
%         h1_setGlobal_s_dgp_penalty(2);
%         set_scenario = 3;
%        
% %         for subscenario = 1:10
% %             
% %             if subscenario == 1
% %                 h1_setGlobal_garage_usage_information_parameter(0.1)
% %                 set_subscenario = 1;
% %                 
% %             elseif subscenario == 2
% %                 h1_setGlobal_garage_usage_information_parameter(0.2)
% %                 set_subscenario = 2;
% %                 
% %             elseif subscenario == 3
% %                 h1_setGlobal_garage_usage_information_parameter(0.3)
% %                 set_subscenario = 3;
% %                 
% %             elseif subscenario == 4
% %                 h1_setGlobal_garage_usage_information_parameter(0.4)
% %                 set_subscenario = 4;
% %                 
% %             elseif subscenario == 5
% %                 h1_setGlobal_garage_usage_information_parameter(0.5)
% %                 set_subscenario = 5;
% %                 
% %             elseif subscenario == 6
% %                 h1_setGlobal_garage_usage_information_parameter(0.6)
% %                 set_subscenario = 6;
% %                 
% %             elseif subscenario == 7
% %                 h1_setGlobal_garage_usage_information_parameter(0.7)
% %                 set_subscenario = 7;
% %                 
% %             elseif subscenario == 8
% %                 h1_setGlobal_garage_usage_information_parameter(0.8)
% %                 set_subscenario = 8;
% %                 
% %             elseif subscenario == 9
% %                 h1_setGlobal_garage_usage_information_parameter(0.9)
% %                 set_subscenario = 9;
% %                 
% %             elseif subscenario == 10
% %                 h1_setGlobal_garage_usage_information_parameter(1)
% %                 set_subscenario = 10;
% %                 
% %             end
% %         end               
            
%     elseif scenario == 4
% %       25% of on-street parking spaces (207) = 52 to garage parking spaces:
%         h1_setGlobal_number_of_parking_garages(2)
%         h1_setGlobal_capacity_garage(192) % 2 x 26 = 52 garage parking spaces
%   
%         h1_setGlobal_number_of_on_street_parking_spaces(155) % 207 - 52 = 155
%         h1_setGlobal_garage_usage_information(0)
%         h1_setGlobal_garage_usage_information_parameter(0.5)
%         h1_setGlobal_s_dgp_penalty(2);
%         set_scenario = 4;
        
    end


average_searching_time(:,:) = zeros(10,10); 
total_searching_time(:,:) = zeros(10,10); 
average_non_searching_time_parkers(:,:) = zeros(10,10); 
total_non_searching_time_parkers(:,:) = zeros(10,10); 
average_searching_distance(:,:) = zeros(10,10); 
total_searching_distance(:,:) = zeros(10,10); 
average_non_searching_distance_parkers(:,:) = zeros(10,10); 
total_non_searching_distance_parkers(:,:) = zeros(10,10);
average_deciding_gp_time(:,:) = zeros(10,10); 
total_deciding_gp_time(:,:) = zeros(10,10);
average_deciding_gp_distance(:,:) = zeros(10,10);
total_deciding_gp_distance(:,:) = zeros(10,10);
avg_total_time(:,:) = zeros(10,10); 
tot_time(:,:) = zeros(10,10); 
avg_total_distance(:,:) = zeros(10,10); 
tot_distance(:,:) = zeros(10,10); 
final_cum_revenue(:,:) = zeros(10,10); 
final_cum_revenue_garage(:,:) = zeros(10,10); 
total_revenue(:,:) = zeros(10,10);
total_on_street_revenue(:,:) = zeros(10,10);
total_garage_revenue(:,:) = zeros(10,10);
enterthearea_1(:,:,:,:) = zeros(1440,1,10,10);
enterthearea_2(:,:,:,:) = zeros(1440,1,10,10);
enterthearea_3(:,:,:,:) = zeros(1440,1,10,10);
enterthearea(:,:,:,:) = zeros(1440,1,10,10); 
n_vehicles_ns_dgp_VOT(:,:,:,:) = zeros(1440,4,10,10);
starttosearch_2(:,:,:,:) = zeros(1440,1,10,10);
starttosearch_3(:,:,:,:) = zeros(1440,1,10,10);
starttosearch(:,:,:,:) = zeros(1440,1,10,10);
n_vehicles_s_dgp_VOT(:,:,:,:) = zeros(1440,4,10,10);
findparking_2(:,:,:,:) = zeros(1440,1,10,10);
findparking_3(:,:,:,:) = zeros(1440,1,10,10);
findparking(:,:,:,:) = zeros(1440,1,10,10);
n_enter_garage(:,:,:,:) = zeros(1440,4,10,10);
n_dgp_searching(:,:,:,:) = zeros(1440,4,10,10); 
departparking_2(:,:,:,:) = zeros(1440,1,10,10); 
departparking_3(:,:,:,:) = zeros(1440,1,10,10); 
departparking(:,:,:,:) = zeros(1440,1,10,10); 
ndepart_garage(:,:,:,:) = zeros(1440,1,10,10); 
leavethearea_1(:,:,:,:) = zeros(1440,1,10,10); 
leavethearea_2(:,:,:,:) = zeros(1440,1,10,10); 
leavethearea_3(:,:,:,:) = zeros(1440,1,10,10); 
leavethearea(:,:,:,:) = zeros(1440,1,10,10);

% set global variable initial on-street parking pricing:
for j = 0.5:0.5:5
    
    h1_setGlobal_initial_parking_pricing(j);
    m = j/0.5;
    
    % set global variable initial garage parking pricing:
    for k = 0.5:0.5:5
    
        clearvars -except global_initial_garage_parking_pricing k i j m average_searching_time total_searching_time average_non_searching_time_parkers total_non_searching_time_parkers...
        average_searching_distance total_searching_distance average_non_searching_distance_parkers total_non_searching_distance_parkers ...
        average_deciding_gp_time total_deciding_gp_time average_deciding_gp_distance total_deciding_gp_distance ...
        avg_total_time tot_time avg_total_distance tot_distance final_cum_revenue final_cum_revenue_garage total_revenue total_on_street_revenue total_garage_revenue ...
        enterthearea_1 enterthearea_2 enterthearea_3 enterthearea n_vehicles_ns_dgp_VOT starttosearch_2 starttosearch_3 starttosearch ...
        n_vehicles_s_dgp_VOT findparking_2 findparking_3 findparking n_enter_garage n_dgp_searching ...
        departparking_2 departparking_3 departparking ndepart_garage ...
        leavethearea_1 leavethearea_2 leavethearea_3 leavethearea set_scenario
        
        h1_setGlobal_initial_garage_parking_pricing(k);
        i = k/0.5;
        [average_searching_time(m,i), total_searching_time(m,i), average_non_searching_time_parkers(m,i), total_non_searching_time_parkers(m,i),...
        average_searching_distance(m,i), total_searching_distance(m,i), average_non_searching_distance_parkers(m,i), total_non_searching_distance_parkers(m,i), ...
        average_deciding_gp_time(m,i), total_deciding_gp_time(m,i), average_deciding_gp_distance(m,i), total_deciding_gp_distance(m,i),...
        avg_total_time(m,i), tot_time(m,i), avg_total_distance(m,i), tot_distance(m,i), final_cum_revenue(m,i), final_cum_revenue_garage(m,i), total_revenue(m,i), total_on_street_revenue(m,i), total_garage_revenue(m,i), ...
        enterthearea_1(:,:,m,i), enterthearea_2(:,:,m,i), enterthearea_3(:,:,m,i), enterthearea(:,:,m,i), n_vehicles_ns_dgp_VOT(:,:,m,i), starttosearch_2(:,:,m,i), starttosearch_3(:,:,m,i), starttosearch(:,:,m,i),...
        n_vehicles_s_dgp_VOT(:,:,m,i), findparking_2(:,:,m,i), findparking_3(:,:,m,i), findparking(:,:,m,i), n_enter_garage(:,:,m,i), n_dgp_searching(:,:,m,i),...
        departparking_2(:,:,m,i), departparking_3(:,:,m,i), departparking(:,:,m,i), ndepart_garage(:,:,m,i), ...
        leavethearea_1(:,:,m,i), leavethearea_2(:,:,m,i), leavethearea_3(:,:,m,i), leavethearea(:,:,m,i)] = Marcopark();
        
        j
        k
    end
    
end


if set_scenario == 1
    save('results_scenario_1_capacity_limited.mat')
elseif set_scenario == 2
    save('results_scenario_2_capacity_larger_demand.mat')
elseif set_scenario == 3
    save('results_scenario_3_availability_garage_information.mat')
    
%     if set_subscenario == 1
%         save('results_scenario_3_1_availability_garage_information.mat')
%     elseif set_subscenario == 2
%         save('results_scenario_3_2_availability_garage_information.mat')
%     elseif set_subscenario == 3
%         save('results_scenario_3_3_availability_garage_information.mat')
%     elseif set_subscenario == 4
%         save('results_scenario_3_4_availability_garage_information.mat')
%     elseif set_subscenario == 5
%         save('results_scenario_3_5_availability_garage_information.mat')
%     elseif set_subscenario == 6
%         save('results_scenario_3_6_availability_garage_information.mat')
%     elseif set_subscenario == 7
%         save('results_scenario_3_7_availability_garage_information.mat')
%     elseif set_subscenario == 8
%         save('results_scenario_3_8_availability_garage_information.mat')
%     elseif set_subscenario == 9
%         save('results_scenario_3_9_availability_garage_information.mat')
%     elseif set_subscenario == 10
%         save('results_scenario_3_10_availability_garage_information.mat')
%     end
    
elseif set_scenario == 4
    save('results_scenario_4_converting_on_street_to_garage.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.) Average searching time:
% set on-street parking price to 6.5 (constant):
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:10;
y_value = average_searching_time(13,:);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
% title('Garage parking price vs. avg searching time')
xlabel('Garage parking price (in CHF)')
ylabel('Average searching time (in min/veh)')
set(gca,'FontSize',24)
% axis([0,10,5,18])

if set_scenario == 1
    saveas(gcf,'(1) Avg_searching_time - On_Street_Constant.jpg')
    saveas(gcf,'(1) Avg_searching_time - On_Street_Constant.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Avg_searching_time - On_Street_Constant.jpg')
    saveas(gcf,'(2) Avg_searching_time - On_Street_Constant.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Avg_searching_time - On_Street_Constant.jpg')
    saveas(gcf,'(3) Avg_searching_time - On_Street_Constant.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_1) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_2) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_3) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_4) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_5) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_6) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_7) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_8) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_9) Avg_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Avg_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_10) Avg_searching_time - On_Street_Constant.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Avg_searching_time - On_Street_Constant.jpg')
    saveas(gcf,'(4) Avg_searching_time - On_Street_Constant.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.) Total searching time:
% set on-street parking price to 6.5 (constant):
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x2_value = 0.5:0.5:10;
y2_value = total_searching_time(13,:);

plot(x2_value, y2_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
% title('Garage parking price vs. total searching time')
xlabel('Garage parking price (in CHF)')
ylabel('Total searching time (in min)')
set(gca,'FontSize',24)
% axis([0,10,0,100])

if set_scenario == 1
    saveas(gcf,'(1) Total_searching_time - On_Street_Constant.jpg')
    saveas(gcf,'(1) Total_searching_time - On_Street_Constant.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Total_searching_time - On_Street_Constant.jpg')
    saveas(gcf,'(2) Total_searching_time - On_Street_Constant.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Total_searching_time - On_Street_Constant.jpg')
    saveas(gcf,'(3) Total_searching_time - On_Street_Constant.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_1) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_2) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_3) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_4) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_5) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_6) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_7) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_8) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_9) Total_searching_time - On_Street_Constant.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Total_searching_time - On_Street_Constant.jpg')
%         saveas(gcf,'(3_10) Total_searching_time - On_Street_Constant.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Total_searching_time - On_Street_Constant.jpg')
    saveas(gcf,'(4) Total_searching_time - On_Street_Constant.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.) Average searching time:
% set garage parking price to 2.5 (constant):
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x3_value = 0.5:0.5:10;
y3_value = average_searching_time(:,5);

plot(x3_value, y3_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
% title('Garage parking price vs. avg searching time')
xlabel('On-street parking price (in CHF)')
ylabel('Average searching time (in min/veh)')
set(gca,'FontSize',24)
% axis([0,10,5.5,7])

if set_scenario == 1
    saveas(gcf,'(1) Avg_searching_time - Garage_Constant.jpg')
    saveas(gcf,'(1) Avg_searching_time - Garage_Constant.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Avg_searching_time - Garage_Constant.jpg')
    saveas(gcf,'(2) Avg_searching_time - Garage_Constant.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Avg_searching_time - Garage_Constant.jpg')
    saveas(gcf,'(3) Avg_searching_time - Garage_Constant.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_1) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_2) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_3) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_4) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_5) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_6) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_7) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_8) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_9) Avg_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Avg_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_10) Avg_searching_time - Garage_Constant.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Avg_searching_time - Garage_Constant.jpg')
    saveas(gcf,'(4) Avg_searching_time - Garage_Constant.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.) Total searching time:
% set garage parking price to 2.5 (constant):
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x4_value = 0.5:0.5:10;
y4_value = total_searching_time(:,5);

plot(x4_value, y4_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
% title('Garage parking price vs. total searching time')
xlabel('On-street parking price (in CHF)')
ylabel('Total searching time (in min)')
set(gca,'FontSize',24)
% axis([0,10,0,100])

if set_scenario == 1
    saveas(gcf,'(1) Total_searching_time - Garage_Constant.jpg')
    saveas(gcf,'(1) Total_searching_time - Garage_Constant.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Total_searching_time - Garage_Constant.jpg')
    saveas(gcf,'(2) Total_searching_time - Garage_Constant.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Total_searching_time - Garage_Constant.jpg')
    saveas(gcf,'(3) Total_searching_time - Garage_Constant.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_1) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_2) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_3) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_4) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_5) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_6) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_7) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_8) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_9) Total_searching_time - Garage_Constant.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Total_searching_time - Garage_Constant.jpg')
%         saveas(gcf,'(3_10) Total_searching_time - Garage_Constant.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Total_searching_time - Garage_Constant.jpg')
    saveas(gcf,'(4) Total_searching_time - Garage_Constant.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5a.) Average searching time:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

op_price = 0.5:0.5:5;
gp_price = 0.5:0.5:5;
x5_value = zeros(100,1);
y5_value = zeros(100,1);
for r = 1:10
    for s = 1:10
        x5_value(r+s-1 + 9*(r-1),1) = op_price(s)./gp_price(r);
        y5_value(r+s-1 + 9*(r-1),1) = average_searching_time(s,r);
    end
end

[x5_sorted, x5_order] = sort(x5_value);
y5_sorted = y5_value(x5_order,:);

x5_sorted_final(1,1) = x5_sorted(1,1);
y5_sorted_final(1,1) = y5_sorted(1,1);
v = 2;
for u = 2:size(x5_sorted,1)
    if x5_sorted(u,1) == x5_sorted(u-1,1)
       x5_sorted_final(v-1,1) = x5_sorted(u,1);
       y5_sorted_final(v-1,1) =  y5_sorted(u,1);
    else
        x5_sorted_final(v,1) =  x5_sorted(u,1);
        y5_sorted_final(v,1) =  y5_sorted(u,1);
        v = v + 1;
    end
end

plot(x5_sorted_final(:,1), y5_sorted_final(:,1), 'b-', 'MarkerSize', 10, 'LineWidth', 4)

xlabel('On-street parking price / Garage parking price')
ylabel('Average searching time (in min/veh)')
% zlabel('Average searching time (in min/veh)')
% axis([0,10,0,10,5.7,21])
% caxis([5.7 21]);
set(gca,'FontSize',24)
% c = colorbar;
% c.Label.String = 'Average searching time (in min)';
% set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%16.f')) )

% saveas(gcf,'Avg_searching_time - All_values.jpg')
% saveas(gcf,'Avg_searching_time - All_values.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5b.) Average dgp time:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

op_price = 0.5:0.5:5;
gp_price = 0.5:0.5:5;
x5b_value = zeros(100,1);
y5b_value = zeros(100,1);
for r = 1:10
    for s = 1:10
        x5b_value(r+s-1 + 9*(r-1),1) = op_price(s)./gp_price(r);
        y5b_value(r+s-1 + 9*(r-1),1) = average_deciding_gp_time(s,r);
    end
end

[x5b_sorted, x5b_order] = sort(x5b_value);
y5b_sorted = y5b_value(x5b_order,:);

x5b_sorted_final(1,1) = x5b_sorted(1,1);
y5b_sorted_final(1,1) = y5b_sorted(1,1);
v = 2;
for u = 2:size(x5b_sorted,1)
    if x5b_sorted(u,1) == x5b_sorted(u-1,1)
       x5b_sorted_final(v-1,1) = x5b_sorted(u,1);
       y5b_sorted_final(v-1,1) =  y5b_sorted(u,1);
    else
        x5b_sorted_final(v,1) =  x5b_sorted(u,1);
        y5b_sorted_final(v,1) =  y5b_sorted(u,1);
        v = v + 1;
    end
end

plot(x5b_sorted_final(:,1), y5b_sorted_final(:,1), 'b-', 'MarkerSize', 10, 'LineWidth', 4)

xlabel('On-street parking price / Garage parking price')
ylabel('Average dgp time (in min/veh)')
% zlabel('Average searching time (in min/veh)')
% axis([0,10,0,10,5.7,21])
% caxis([5.7 21]);
set(gca,'FontSize',24)
% c = colorbar;
% c.Label.String = 'Average searching time (in min)';
% set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%16.f')) )

% saveas(gcf,'Avg_searching_time - All_values.jpg')
% saveas(gcf,'Avg_searching_time - All_values.fig')

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5c.) Average searching and dgp time:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

plot(x5_sorted_final(:,1), y5_sorted_final(:,1), 'b--', 'MarkerSize', 10, 'LineWidth', 4)
hold on
plot(x5b_sorted_final(:,1), y5b_sorted_final(:,1), 'r-', 'MarkerSize', 10, 'LineWidth', 4)
hold on
y5c_sorted_final(:,1) = y5_sorted_final(:,1) + y5b_sorted_final(:,1);
plot(x5_sorted_final(:,1), y5c_sorted_final(:,1), 'g:', 'MarkerSize', 10, 'LineWidth', 4)

legend('Vehicles searching for on-street parking','Vehicles driving to garage parking','All vehicles searching for on-street parking and driving to garage parking')
% [lgd,icons,plots,txt] = legend('Vehicles searching for on-street parking','Vehicles driving to garage parking','All vehicles searching for on-street parking and driving to garage parking');

xlabel('On-street parking price / Garage parking price')
ylabel('Average time (in min/veh)')
axis([0,10,0,8])
set(gca,'FontSize',24)

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.) Total searching time:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

surf(0.5:0.5:10,0.5:0.5:10,total_searching_time(:,:))

% title('On-street and garage parking price vs. avg searching time')
xlabel('Garage parking price (in CHF)')
ylabel('On-street parking price (in CHF)')
zlabel('Total searching time (in min)')
set(gca,'FontSize',24)
c = colorbar;
c.Label.String = 'Total searching time (in min)';
set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%16.f')) )

if set_scenario == 1
    saveas(gcf,'(1) Total_searching_time - All_values.jpg')
    saveas(gcf,'(1) Total_searching_time - All_values.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Total_searching_time - All_values.jpg')
    saveas(gcf,'(2) Total_searching_time - All_values.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Total_searching_time - All_values.jpg')
    saveas(gcf,'(3) Total_searching_time - All_values.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_1) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_2) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_3) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_4) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_5) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_6) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_7) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_8) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_9) Total_searching_time - All_values.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Total_searching_time - All_values.jpg')
%         saveas(gcf,'(3_10) Total_searching_time - All_values.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Total_searching_time - All_values.jpg')
    saveas(gcf,'(4) Total_searching_time - All_values.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7.) Average searching distance:
% set on-street parking price to 6.5 (constant):
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x_value = 0.5:0.5:10;
y_value = average_searching_distance(13,:);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('Garage parking price (in CHF)')
ylabel('Average searching distance (in km/veh)')
set(gca,'FontSize',24)
% axis([0,10,5,18])

if set_scenario == 1
    saveas(gcf,'(1) Avg_searching_distance - On_Street_Constant.jpg')
    saveas(gcf,'(1) Avg_searching_distance - On_Street_Constant.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Avg_searching_distance - On_Street_Constant.jpg')
    saveas(gcf,'(2) Avg_searching_distance - On_Street_Constant.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Avg_searching_distance - On_Street_Constant.jpg')
    saveas(gcf,'(3) Avg_searching_distance - On_Street_Constant.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_1) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_2) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_3) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_4) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_5) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_6) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_7) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_8) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_9) Avg_searching_distance - On_Street_Constant.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Avg_searching_distance - On_Street_Constant.jpg')
%         saveas(gcf,'(3_10) Avg_searching_distance - On_Street_Constant.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Avg_searching_distance - On_Street_Constant.jpg')
    saveas(gcf,'(4) Avg_searching_distance - On_Street_Constant.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.) Average searching distance:
% set garage parking price to 2.5 (constant):
figure('units','normalized','outerposition',[0 0 1 1])
hold on
x3_value = 0.5:0.5:10;
y3_value = average_searching_distance(:,5);

plot(x3_value, y3_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)
xlabel('On-street parking price (in CHF)')
ylabel('Average searching distance (in km/veh)')
set(gca,'FontSize',24)
% axis([0,10,0,10])

if set_scenario == 1
    saveas(gcf,'(1) Avg_searching_distance - Garage_Constant.jpg')
    saveas(gcf,'(1) Avg_searching_distance - Garage_Constant.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Avg_searching_distance - Garage_Constant.jpg')
    saveas(gcf,'(2) Avg_searching_distance - Garage_Constant.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Avg_searching_distance - Garage_Constant.jpg')
    saveas(gcf,'(3) Avg_searching_distance - Garage_Constant.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_1) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_2) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_3) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_4) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_5) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_6) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_7) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_8) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_9) Avg_searching_distance - Garage_Constant.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Avg_searching_distance - Garage_Constant.jpg')
%         saveas(gcf,'(3_10) Avg_searching_distance - Garage_Constant.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Avg_searching_distance - Garage_Constant.jpg')
    saveas(gcf,'(4) Avg_searching_distance - Garage_Constant.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 9.) Average searching distance:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

surf(0.5:0.5:10,0.5:0.5:10,average_searching_distance(:,:))

xlabel('Garage parking price (in CHF)')
ylabel('On-street parking price (in CHF)')
zlabel('Average searching distance (in km/veh)')
% caxis([1.2 4.4]);
set(gca,'FontSize',24)
c = colorbar;
c.Label.String = 'Average searching distance (in km/veh)';
set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%16.f')) )

if set_scenario == 1
    saveas(gcf,'(1) Avg_searching_distance - All_values.jpg')
    saveas(gcf,'(1) Avg_searching_distance - All_values.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Avg_searching_distance - All_values.jpg')
    saveas(gcf,'(2) Avg_searching_distance - All_values.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Avg_searching_distance - All_values.jpg')
    saveas(gcf,'(3) Avg_searching_distance - All_values.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_1) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_2) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_3) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_4) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_5) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_6) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_7) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_8) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_9) Avg_searching_distance - All_values.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Avg_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_10) Avg_searching_distance - All_values.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Avg_searching_distance - All_values.jpg')
    saveas(gcf,'(4) Avg_searching_distance - All_values.fig')
end


hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9b.) Average dgp distance:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

op_price = 0.5:0.5:5;
gp_price = 0.5:0.5:5;
x9b_value = zeros(100,1);
y9b_value = zeros(100,1);
for r = 1:10
    for s = 1:10
        x9b_value(r+s-1 + 19*(r-1),1) = op_price(s)./gp_price(r);
        y9b_value(r+s-1 + 19*(r-1),1) = average_deciding_gp_distance(s,r);
    end
end

[x9b_sorted, x9b_order] = sort(x9b_value);
y9b_sorted = y9b_value(x9b_order,:);

x9b_sorted_final(1,1) = x9b_sorted(1,1);
y9b_sorted_final(1,1) = y9b_sorted(1,1);
v = 2;
for u = 2:size(x9b_sorted,1)
    if x9b_sorted(u,1) == x9b_sorted(u-1,1)
       x9b_sorted_final(v-1,1) = x9b_sorted(u,1);
       y9b_sorted_final(v-1,1) =  y9b_sorted(u,1);
    else
        x9b_sorted_final(v,1) =  x9b_sorted(u,1);
        y9b_sorted_final(v,1) =  y9b_sorted(u,1);
        v = v + 1;
    end
end

plot(x9b_sorted_final(:,1), y9b_sorted_final(:,1), 'b-', 'MarkerSize', 10, 'LineWidth', 4)

xlabel('On-street parking price / Garage parking price')
ylabel('Average dgp distance (in km/veh)')
set(gca,'FontSize',24)

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10.) Total searching distance:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

surf(0.5:0.5:10,0.5:0.5:10,total_searching_distance(:,:))

xlabel('Garage parking price (in CHF)')
ylabel('On-street parking price (in CHF)')
zlabel('Total searching distance (in km)')
set(gca,'FontSize',24)
c = colorbar;
c.Label.String = 'Total searching distance (in km/veh)';
set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%16.f')) )

if set_scenario == 1
    saveas(gcf,'(1) Total_searching_distance - All_values.jpg')
    saveas(gcf,'(1) Total_searching_distance - All_values.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Total_searching_distance - All_values.jpg')
    saveas(gcf,'(2) Total_searching_distance - All_values.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Total_searching_distance - All_values.jpg')
    saveas(gcf,'(3) Total_searching_distance - All_values.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_1) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_2) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_3) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_4) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_5) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_6) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_7) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_8) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_9) Total_searching_distance - All_values.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Total_searching_distance - All_values.jpg')
%         saveas(gcf,'(3_10) Total_searching_distance - All_values.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Total_searching_distance - All_values.jpg')
    saveas(gcf,'(4) Total_searching_distance - All_values.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 11.) Total on-street parking revenue:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

surf(0.5:0.5:10,0.5:0.5:10,total_on_street_revenue(:,:))

xlabel('Garage parking price (in CHF)')
ylabel('On-street parking price (in CHF)')
zlabel('Total on-street parking revenue (in CHF)')
% caxis([457 9131]);
set(gca,'FontSize',24)
c = colorbar;
c.Label.String = 'Total on-street parking revenue (in CHF)';
set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%16.f')) )

if set_scenario == 1
    saveas(gcf,'(1) Total_On_Street_Parking_Revenue.jpg')
    saveas(gcf,'(1) Total_On_Street_Parking_Revenue.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Total_On_Street_Parking_Revenue.jpg')
    saveas(gcf,'(2) Total_On_Street_Parking_Revenue.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Total_On_Street_Parking_Revenue.jpg')
    saveas(gcf,'(3) Total_On_Street_Parking_Revenue.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_1) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_2) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_3) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_4) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_5) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_6) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_7) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_8) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_9) Total_On_Street_Parking_Revenue.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Total_On_Street_Parking_Revenue.jpg')
%         saveas(gcf,'(3_10) Total_On_Street_Parking_Revenue.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Total_On_Street_Parking_Revenue.jpg')
    saveas(gcf,'(4) Total_On_Street_Parking_Revenue.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 12.) Total garage parking revenue:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

surf(0.5:0.5:10,0.5:0.5:10,total_garage_revenue(:,:))

xlabel('Garage parking price (in CHF)')
ylabel('On-street parking price (in CHF)')
zlabel('Total garage parking revenue (in CHF)')
% caxis([578 11930]);
set(gca,'FontSize',24)
c = colorbar;
c.Label.String = 'Total garage parking revenue (in CHF)';
set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%16.f')) )

if set_scenario == 1
    saveas(gcf,'(1) Total_Garage_Parking_Revenue.jpg')
    saveas(gcf,'(1) Total_Garage_Parking_Revenue.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Total_Garage_Parking_Revenue.jpg')
    saveas(gcf,'(2) Total_Garage_Parking_Revenue.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Total_Garage_Parking_Revenue.jpg')
    saveas(gcf,'(3) Total_Garage_Parking_Revenue.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_1) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_2) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_3) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_4) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_5) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_6) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_7) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_8) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_9) Total_Garage_Parking_Revenue.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Total_Garage_Parking_Revenue.jpg')
%         saveas(gcf,'(3_10) Total_Garage_Parking_Revenue.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Total_Garage_Parking_Revenue.jpg')
    saveas(gcf,'(4) Total_Garage_Parking_Revenue.fig')
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 13.) Total on- and off-street parking revenue:
% set on-street and garage parking pricing as variable:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

surf(0.5:0.5:5,0.5:0.5:5,total_revenue(:,:))

xlabel('Garage parking price (in CHF)')
ylabel('On-street parking price (in CHF)')
zlabel('Total on-street and garage parking revenue (in CHF)')
% caxis([1035 20690]);
set(gca,'FontSize',24)
c = colorbar;
c.Label.String = 'Total on-street and garage parking revenue (in CHF)';
set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%16.f')) )

if set_scenario == 1
    saveas(gcf,'(1) Total_Revenue.jpg')
    saveas(gcf,'(1) Total_Revenue.fig')
elseif set_scenario == 2
    saveas(gcf,'(2) Total_Revenue.jpg')
    saveas(gcf,'(2) Total_Revenue.fig')
elseif set_scenario == 3
    saveas(gcf,'(3) Total_Revenue.jpg')
    saveas(gcf,'(3) Total_Revenue.fig')
    
%     if set_subscenario == 1
%         saveas(gcf,'(3_1) Total_Revenue.jpg')
%         saveas(gcf,'(3_1) Total_Revenue.fig')
%     elseif set_subscenario == 2
%         saveas(gcf,'(3_2) Total_Revenue.jpg')
%         saveas(gcf,'(3_2) Total_Revenue.fig')
%     elseif set_subscenario == 3
%         saveas(gcf,'(3_3) Total_Revenue.jpg')
%         saveas(gcf,'(3_3) Total_Revenue.fig')
%     elseif set_subscenario == 4
%         saveas(gcf,'(3_4) Total_Revenue.jpg')
%         saveas(gcf,'(3_4) Total_Revenue.fig')
%     elseif set_subscenario == 5
%         saveas(gcf,'(3_5) Total_Revenue.jpg')
%         saveas(gcf,'(3_5) Total_Revenue.fig')
%     elseif set_subscenario == 6
%         saveas(gcf,'(3_6) Total_Revenue.jpg')
%         saveas(gcf,'(3_6) Total_Revenue.fig')
%     elseif set_subscenario == 7
%         saveas(gcf,'(3_7) Total_Revenue.jpg')
%         saveas(gcf,'(3_7) Total_Revenue.fig')
%     elseif set_subscenario == 8
%         saveas(gcf,'(3_8) Total_Revenue.jpg')
%         saveas(gcf,'(3_8) Total_Revenue.fig')
%     elseif set_subscenario == 9
%         saveas(gcf,'(3_9) Total_Revenue.jpg')
%         saveas(gcf,'(3_9) Total_Revenue.fig')
%     elseif set_subscenario == 10
%         saveas(gcf,'(3_10) Total_Revenue.jpg')
%         saveas(gcf,'(3_10) Total_Revenue.fig')
%     end
    
elseif set_scenario == 4
    saveas(gcf,'(4) Total_Revenue.jpg')
    saveas(gcf,'(4) Total_Revenue.fig')
end

hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 14.) ALL revenues:
figure('units','normalized','outerposition',[0 0 1 1])
hold on

op_price = 0.5:0.5:5;
gp_price = 0.5:0.5:5;
x14a_value = zeros(100,1);
y14a_value = zeros(100,1);
for r = 1:10
    for s = 1:10
        x14a_value(r+s-1 + 19*(r-1),1) = op_price(s)./gp_price(r);
        y14a_value(r+s-1 + 19*(r-1),1) = total_on_street_revenue(s,r);
    end
end

[x14a_sorted, x14a_order] = sort(x14a_value);
y14a_sorted = y14a_value(x14a_order,:);

x14a_sorted_final(1,1) = x14a_sorted(1,1);
y14a_sorted_final(1,1) = y14a_sorted(1,1);
v = 2;
% for u = 2:size(x14a_sorted,1)
%     if x14a_sorted(u,1) == x14a_sorted(u-1,1)
%        x14a_sorted_final(v-1,1) = x14a_sorted(u,1);
%        y14a_sorted_final(v-1,1) =  y14a_sorted(u,1);
%     else
%         x14a_sorted_final(v,1) =  x14a_sorted(u,1);
%         y14a_sorted_final(v,1) =  y14a_sorted(u,1);
%         v = v + 1;
%     end
% end

% x_value = 0:0.05:1;
% y_value = total_on_street_revenue(:,1);
plot(x14a_sorted_final, y14a_sorted_final, 'b--', 'MarkerSize', 10, 'LineWidth', 4)
hold on


op_price = 0.5:0.5:5;
gp_price = 0.5:0.5:5;
x14b_value = zeros(100,1);
y14b_value = zeros(100,1);
for r = 1:10
    for s = 1:10
        x14b_value(r+s-1 + 19*(r-1),1) = op_price(s)./gp_price(r);
        y14b_value(r+s-1 + 19*(r-1),1) = total_garage_revenue(s,r);
    end
end

[x14b_sorted, x14b_order] = sort(x14b_value);
y14b_sorted = y14b_value(x14b_order,:);

x14b_sorted_final(1,1) = x14b_sorted(1,1);
y14b_sorted_final(1,1) = y14b_sorted(1,1);
v = 2;
% for u = 2:size(x14b_sorted,1)
%     if x14b_sorted(u,1) == x14b_sorted(u-1,1)
%        x14b_sorted_final(v-1,1) = x14b_sorted(u,1);
%        y14b_sorted_final(v-1,1) =  y14b_sorted(u,1);
%     else
%         x14b_sorted_final(v,1) =  x14b_sorted(u,1);
%         y14b_sorted_final(v,1) =  y14b_sorted(u,1);
%         v = v + 1;
%     end
% end

% x_value = 0:0.05:1;
% y_value = total_garage_revenue(:,1);
plot(x14b_sorted_final, y14b_sorted_final, 'g:', 'MarkerSize', 10, 'LineWidth', 4)
hold on

op_price = 0.5:0.5:5;
gp_price = 0.5:0.5:5;
x14c_value = zeros(100,1);
y14c_value = zeros(100,1);
for r = 1:10
    for s = 1:10
        x14c_value(r+s-1 + 19*(r-1),1) = op_price(s)./gp_price(r);
        y14c_value(r+s-1 + 19*(r-1),1) = total_revenue(s,r);
    end
end

[x14c_sorted, x14c_order] = sort(x14c_value);
y14c_sorted = y14c_value(x14c_order,:);

x14c_sorted_final(1,1) = x14c_sorted(1,1);
y14c_sorted_final(1,1) = y14c_sorted(1,1);
v = 2;
% for u = 2:size(x14c_sorted,1)
%     if x14c_sorted(u,1) == x14c_sorted(u-1,1)
%        x14c_sorted_final(v-1,1) = x14c_sorted(u,1);
%        y14c_sorted_final(v-1,1) =  y14c_sorted(u,1);
%     else
%         x14c_sorted_final(v,1) =  x14c_sorted(u,1);
%         y14c_sorted_final(v,1) =  y14c_sorted(u,1);
%         v = v + 1;
%     end
% end

% x_value = 0:0.05:1;
% y_value = total_revenue(:,1);
plot(x14c_sorted_final, y14c_sorted_final, 'r-', 'MarkerSize', 10, 'LineWidth', 4)
hold on

xlabel('On-street parking price / Garage parking price')
ylabel('Total revenue (in CHF)')
legend('Total on-street parking revenue','Total garage parking revenue','Total revenue created by both on-street and garage parking')
set(gca,'FontSize',24)

axis([0 1 0 inf])
ax = gca;
ax.XTick = 0:0.1:1;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(xticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsy]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
set(gca,'FontSize',24)

% saveas(gcf,'Total_revenue_converting_on_street_to_garage.jpg')
% saveas(gcf,'Total_revenue_converting_on_street_to_garage.fig')

hold off
