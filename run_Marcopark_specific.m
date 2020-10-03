% Run Marcopark Specific:
% Specific on-street and garage parking prices vs. avg searching time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% h1_setGlobal_initial_parking_pricing(1.5);
% h1_setGlobal_initial_garage_parking_pricing(3);
% h1_setGlobal_number_of_parking_garages(2);
% h1_setGlobal_capacity_garage(166);
% h1_setGlobal_number_of_on_street_parking_spaces(207);
% h1_setGlobal_garage_usage_information(0);
% h1_setGlobal_garage_usage_information_parameter(0);
% h1_setGlobal_s_dgp_penalty(2);
% 
% [average_searching_time, total_searching_time, average_non_searching_time_parkers, total_non_searching_time_parkers,...
% average_searching_distance, total_searching_distance, average_non_searching_distance_parkers, total_non_searching_distance_parkers, ...
% average_deciding_gp_time, total_deciding_gp_time, average_deciding_gp_distance, total_deciding_gp_distance,...
% avg_total_time, tot_time, avg_total_distance, tot_distance, final_cum_revenue, final_cum_revenue_garage, total_revenue, total_on_street_revenue, total_garage_revenue, ...
% enterthearea_1, enterthearea_2, enterthearea_3, enterthearea, n_vehicles_ns_dgp_VOT, starttosearch_2, starttosearch_3, starttosearch,...
% n_vehicles_s_dgp_VOT, findparking_2, findparking_3, findparking, n_enter_garage, n_dgp_searching,...
% departparking_2, departparking_3, departparking, ndepart_garage, leavethearea_1, leavethearea_2, leavethearea_3, leavethearea] = Marcopark();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1_setGlobal_initial_parking_pricing(1.5);
h1_setGlobal_initial_garage_parking_pricing(3);
h1_setGlobal_number_of_parking_garages(2);
h1_setGlobal_capacity_garage(166);
h1_setGlobal_number_of_on_street_parking_spaces(207);
h1_setGlobal_garage_usage_information(0);
h1_setGlobal_garage_usage_information_parameter(0.5);
h1_setGlobal_s_dgp_penalty(2);

[average_searching_time, total_searching_time, average_non_searching_time_parkers, total_non_searching_time_parkers,...
average_searching_distance, total_searching_distance, average_non_searching_distance_parkers, total_non_searching_distance_parkers, ...
average_deciding_gp_time, total_deciding_gp_time, average_deciding_gp_distance, total_deciding_gp_distance,...
avg_total_time, tot_time, avg_total_distance, tot_distance, final_cum_revenue, final_cum_revenue_garage, total_revenue, total_on_street_revenue, total_garage_revenue, ...
enterthearea_1, enterthearea_2, enterthearea_3, enterthearea, n_vehicles_ns_dgp_VOT, starttosearch_2, starttosearch_3, starttosearch,...
n_vehicles_s_dgp_VOT, findparking_2, findparking_3, findparking, n_enter_garage, n_dgp_searching,...
departparking_2, departparking_3, departparking, ndepart_garage, leavethearea_1, leavethearea_2, leavethearea_3, leavethearea] = Marcopark();
