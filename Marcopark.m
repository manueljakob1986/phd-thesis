function [average_searching_time, total_searching_time, average_non_searching_time_parkers, total_non_searching_time_parkers,...
    average_searching_distance, total_searching_distance, average_non_searching_distance_parkers, total_non_searching_distance_parkers, ...
    average_deciding_gp_time, total_deciding_gp_time, average_deciding_gp_distance, total_deciding_gp_distance,...
    avg_total_time, tot_time, avg_total_distance, tot_distance, final_cum_revenue, final_cum_revenue_garage, total_revenue, total_on_street_revenue, total_garage_revenue, ...
    enterthearea_1, enterthearea_2, enterthearea_3, enterthearea, n_vehicles_ns_dgp_VOT, starttosearch_2, starttosearch_3, starttosearch,...
    n_vehicles_s_dgp_VOT, findparking_2, findparking_3, findparking, n_enter_garage, n_dgp_searching,...
    departparking_2, departparking_3, departparking, ndepart_garage, leavethearea_1, leavethearea_2, leavethearea_3, leavethearea] = Marcopark()

% clear all% SMALL CIRCLE
% clc

global R t L A_2 A_3 v Lnetwork kc kj Qmax;
    R=0.1; % unit: km
    t=1/60; % STEP SLICE UNIT: hours
    L=7.7; %KM
    A_2=0; % unit: parking spaces in garages
%     A_3=207; % unit: parking spaces on-street
    A_3 = h2_getGlobal_number_of_on_street_parking_spaces; % unit: parking spaces on-street
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % supply decrease by 4%:
%     A_2 = A_2 * 0.96;
%     A_3 = A_3 * 0.96;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % supply increase by 4%:
%     A_2 = A_2 * 1.04;
%     A_3 = A_3 * 1.04;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    v=12.5; % unit: km/hr
    Lnetwork=15.4; %km
    kc=20; % unit: veh/km
    kj=55; % unit: veh/km
    Qmax= v*kc; %unit: veh/hour
% load('ouput_enter.mat')
load('load_traffic_demand')
[demand,~] = hist(act_start_time_seconds./60, linspace(1,1440,1440));
enterarea = demand';
initial_garage_capacity = c19_input_number_of_parking_garages * c20_input_capacity_garage;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % demand decrease by 4%:
% enterarea = enterarea * 0.96;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % demand increase by 4%:
% enterarea = enterarea * 1.04;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('parameter')%parking duration parameters
load('startinginarea') % number of available on-street parking spaces at the beginning
% startinginarea = 183;
load('startinginarea_garage') % number of available garage parking spaces at the beginning
% startinginarea_garage = 113;

T=size(enterarea,1);
% throughtraffic_ratio(1:size(smallcircle_in,1),1)=0; % tt is the percentage of throughtraffic corresponding to each time slice.
% garage_ratio(1:size(smallcircle_in,1),1)=0; % gg is the percentage of garage corresponding to each time slice.

% % parking_garage_information_level = c24_input_parking_garage_information;
switch_on_parking_garage_information = c25_input_switch_on_parking_garage_information;
%%
%ENTER THE AREA--------------------------------
%1st category: throughtraffic
%2nd category: vehicles to garage
%3rd category: vehicles to on-street parking
    enterthearea_1(1:T,1)=enterarea*0.23; 
    enterthearea_2(1:T,1)=0; 
    enterthearea_3(1:T,1)=enterarea*0.77; 
%total:
    enterthearea(:,1)=enterarea; 
%%
    density=0;
    speed=v;  
    availableparking_2 = A_2;
%%%%%%%%%%%%%%%%%
%     startinginarea = 0;
%%%%%%%%%%%%%%%%%
    availableparking_3 = A_3 - startinginarea;
%     availableparking_gp = initial_garage_capacity - startinginarea_garage;
    total_garage_capacity(1,1) = initial_garage_capacity - startinginarea_garage;
    
    Nns_1(1,1)=0;   Nns_2(1,1)=0;    Nns_3(1,1)=0;
    Np_2(1,1)=0;    Np_3(1,1)=startinginarea;
    Ns_2(1,1)=0;    Ns_3(1,1)=0;
    Ns_3_k1(1,1)=0;
    Ns_3_k2(1,1)=0;
    Ns_3_k3(1,1)=0;
    Ns_3_k4(1,1)=0;
    
    Ndgp_k1(1,1)=0;
    Ndgp_k2(1,1)=0;
    Ndgp_k3(1,1)=0;
    Ndgp_k4(1,1)=0;
    Ndgp(1,1) = 0;
    
    Ngp(1,1)=startinginarea_garage;
    
    Nns(1,1)= Nns_1(1,1)+Nns_2(1,1)+Nns_3(1,1); %non-searching can be all, throughtraffic, garage and on-street parkers
    Np(1,1)= Np_2(1,1)+Np_3(1,1); %parking can be both garage and on-street parkers
    Ns(1,1)= Ns_2(1,1)+Ns_3(1,1); %searching can only be on-street parkers
%%
%the first row of every event--------------------------------
    leavethearea_1(1,1) =0;

    starttosearch_2(1,1)=0;
    findparking_2(1,1)  =0; 
    departparking_2(1,1)=0; 
    leavethearea_2(1,1) =0;

    starttosearch_3(1,1)=0; 
    findparking_3(1,1)  =0; 
    decideforparking_3(1,1) =0;
    departparking_3(1,1)=0; 
    leavethearea_3(1,1) =0;
    
    starttosearch(1,1)  = starttosearch_2(1,1)+starttosearch_3(1,1);
    findparking(1,1)    = findparking_2(1,1)+findparking_3(1,1);
    decideforparking(1,1) = decideforparking_3(1,1);
    departparking(1,1)  = departparking_2(1,1)+departparking_3(1,1);
    leavethearea(1,1)   = leavethearea_1(1,1)+ leavethearea_2(1,1)+leavethearea_3(1,1);
    
    parking_probability(1,:) = ones(1,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1_setGlobal_switch_on_damp_exp_coef(0);

% set global variable initial parking pricing to 2.5:
% % h1_setGlobal_initial_parking_pricing(2.5);

% set global variable for maximum price increase per time step to 0.1:  
h1_setGlobal_max_parking_price_increase(0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,n_demand_VOT]= c1_input_n_demand(1);%/
        
 % demand matrix for VOT origins:
 matrix_demand_VOT(1,:) = n_demand_VOT'.*0.77; 
        
% initial price:
 parking_pricing(1,1) = c8_input_parking_price;
 
% initial delta searching:
parking_pricing_switched_on = c13_input_switch_on_parking_pricing;
if parking_pricing_switched_on == 1
    delta_searching_veh_available_parking_spots(1,1) = Ns(1,1)/availableparking_3(1,1);
    delta_searching_veh_available_parking_spots_i_plus_one(1,1) = Ns(1,1)/availableparking_3(1,1);
elseif parking_pricing_switched_on == 3
    delta_searching_veh_available_parking_spots(1,1) = 1/availableparking_3(1,1);
    delta_searching_veh_available_parking_spots_i_plus_one(1,1) = 1/availableparking_3(1,1);
else
    delta_searching_veh_available_parking_spots(1,1) = Ns(1,1)/availableparking_3(1,1);
    delta_searching_veh_available_parking_spots_i_plus_one(1,1) = Ns(1,1)/availableparking_3(1,1);
end

% ns/s matrix for VOT origins:
V = size(n_demand_VOT,1);
matrix_ns_s_VOT(1,:) = zeros(1,V); 

% Garage parking initial transition events:
decideforparking_3_VOT(1,:) = zeros(1,V);
n_vehicles_ns_dgp_VOT(1,:) = zeros(1,V);
n_vehicles_s_dgp_VOT(1,:) = zeros(1,V);
n_enter_garage(1,:) = zeros(1,V);
n_dgp_searching(1,:) = zeros(1,V);
ndepart_garage(1,1) = 0;

% initial bias:
bias = zeros(1440,1);
bias(1,1) = (speed*t)/2 - speed*t;
bias_tot(1,1) = findparking_3(1,1) * bias(1,1);
driven_distance_per_time_slice(1,1) = speed*t*(Nns(1,1) + Ns(1,1));
E_driven_distance(1,1) = driven_distance_per_time_slice(1,1) - bias_tot(1,1);

% set initial value for parking garage capacity:
h1_setGlobal_parking_garage_capacity(c19_input_number_of_parking_garages * c20_input_capacity_garage - Ngp(1,1));

% VOT per user group:
input_value_of_time = c7_input_value_of_time;

% Garage parking initial cost values:
ACT(1) = 0;
C_op(1,:) = zeros(1,4);
C_drive(1,:) = zeros(1,4);
C_walk(1,:) = zeros(1,4);
C_gp(1,:) = zeros(1,4);
share_off_street_parking_decision(1,:) = zeros(1,4);

% Garage parking pricing:
% off_street_parking_pricing_switched_on = c16_input_switch_on_off_street_parking_pricing;
% if off_street_parking_pricing_switched_on == 2
%     garage_parking_pricing(1,1) = c23_input_parking_price_garage;
% elseif off_street_parking_pricing_switched_on == 0
%     garage_parking_pricing(1,1) = 0;
% elseif off_street_parking_pricing_switched_on == 1
%     garage_parking_pricing(1,1) = c23_input_parking_price_garage;
% end
garage_parking_pricing(1,1) = h2_getGlobal_initial_garage_parking_pricing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%
i=1;
while i<1440
    %this part is about the number of vehicles in each state
    Nns_1(i+1,1) = Nns_1(i,1)+enterthearea_1(i,1)-leavethearea_1(i,1); % only throughtraffic
%     Nns_2(i+1,1) = Nns_2(i,1)+enterthearea_2(i,1)-starttosearch_2(i,1)+departparking_2(i,1)-leavethearea_2(i,1); % only garage parked vehicles
%     Nns_3(i+1,1) = Nns_3(i,1)+enterthearea_3(i,1)-starttosearch_3(i,1)+departparking_3(i,1)-leavethearea_3(i,1); % only on-street parked vehicles
    Nns_2(i+1,1) = Nns_2(i,1)+enterthearea_2(i,1)-starttosearch_2(i,1)+departparking_2(i,1)-leavethearea_2(i,1);
    Nns_3(i+1,1) = Nns_3(i,1)+enterthearea_3(i,1)-starttosearch_3(i,1)+departparking_3(i,1)-leavethearea_3(i,1) + ndepart_garage(i,1) - sum(n_vehicles_ns_dgp_VOT(i,:),2); % only on-street parked vehicles
    
    Np_2(i+1,1)  = Np_2(i,1)+findparking_2(i,1)-departparking_2(i,1);% only garage parked vehicles
    Np_3(i+1,1)  = Np_3(i,1)+decideforparking_3(i,1)-departparking_3(i,1);% only on-street parked vehicles
    
    Ns_2(i+1,1)  = Ns_2(i,1)+starttosearch_2(i,1)-findparking_2(i,1);
    
    if Ns_3(i,1) == 0
%         Ns_3_k1(i+1,1)  = Ns_3_k1(i,1) + matrix_ns_s_VOT(i,1) - findparking_3(i,1) * 1/4 * parking_probability(i,1);
%         Ns_3_k2(i+1,1)  = Ns_3_k2(i,1) + matrix_ns_s_VOT(i,2) - findparking_3(i,1) * 1/4 * parking_probability(i,2);
%         Ns_3_k3(i+1,1)  = Ns_3_k3(i,1) + matrix_ns_s_VOT(i,3) - findparking_3(i,1) * 1/4 * parking_probability(i,3);
%         Ns_3_k4(i+1,1)  = Ns_3_k4(i,1) + matrix_ns_s_VOT(i,4) - findparking_3(i,1) * 1/4 * parking_probability(i,4);     

        Ns_3_k1(i+1,1)  = Ns_3_k1(i,1) + matrix_ns_s_VOT(i,1) + n_dgp_searching(i,1) - n_vehicles_s_dgp_VOT(i,1) - findparking_3(i,1) * 1/4 * parking_probability(i,1);
        Ns_3_k2(i+1,1)  = Ns_3_k2(i,1) + matrix_ns_s_VOT(i,2) + n_dgp_searching(i,2) - n_vehicles_s_dgp_VOT(i,2) - findparking_3(i,1) * 1/4 * parking_probability(i,2);
        Ns_3_k3(i+1,1)  = Ns_3_k3(i,1) + matrix_ns_s_VOT(i,3) + n_dgp_searching(i,3) - n_vehicles_s_dgp_VOT(i,3) - findparking_3(i,1) * 1/4 * parking_probability(i,3);
        Ns_3_k4(i+1,1)  = Ns_3_k4(i,1) + matrix_ns_s_VOT(i,4) + n_dgp_searching(i,4) - n_vehicles_s_dgp_VOT(i,4) - findparking_3(i,1) * 1/4 * parking_probability(i,4);
        
    else
%         Ns_3_k1(i+1,1)  = Ns_3_k1(i,1) + matrix_ns_s_VOT(i,1) - findparking_3(i,1) * Ns_3_k1(i,1) ./ Ns_3(i,1) * parking_probability(i,1);
%         Ns_3_k2(i+1,1)  = Ns_3_k2(i,1) + matrix_ns_s_VOT(i,2) - findparking_3(i,1) * Ns_3_k2(i,1) ./ Ns_3(i,1) * parking_probability(i,2);
%         Ns_3_k3(i+1,1)  = Ns_3_k3(i,1) + matrix_ns_s_VOT(i,3) - findparking_3(i,1) * Ns_3_k3(i,1) ./ Ns_3(i,1) * parking_probability(i,3);
%         Ns_3_k4(i+1,1)  = Ns_3_k4(i,1) + matrix_ns_s_VOT(i,4) - findparking_3(i,1) * Ns_3_k4(i,1) ./ Ns_3(i,1) * parking_probability(i,4);
 
        Ns_3_k1(i+1,1)  = Ns_3_k1(i,1) + matrix_ns_s_VOT(i,1) + n_dgp_searching(i,1) - n_vehicles_s_dgp_VOT(i,1) - findparking_3(i,1) * Ns_3_k1(i,1) ./ Ns_3(i,1) * parking_probability(i,1);
        Ns_3_k2(i+1,1)  = Ns_3_k2(i,1) + matrix_ns_s_VOT(i,2) + n_dgp_searching(i,2) - n_vehicles_s_dgp_VOT(i,2) - findparking_3(i,1) * Ns_3_k2(i,1) ./ Ns_3(i,1) * parking_probability(i,2);
        Ns_3_k3(i+1,1)  = Ns_3_k3(i,1) + matrix_ns_s_VOT(i,3) + n_dgp_searching(i,3) - n_vehicles_s_dgp_VOT(i,3) - findparking_3(i,1) * Ns_3_k3(i,1) ./ Ns_3(i,1) * parking_probability(i,3);
        Ns_3_k4(i+1,1)  = Ns_3_k4(i,1) + matrix_ns_s_VOT(i,4) + n_dgp_searching(i,4) - n_vehicles_s_dgp_VOT(i,4) - findparking_3(i,1) * Ns_3_k4(i,1) ./ Ns_3(i,1) * parking_probability(i,4);
    end
    
 
    Ns_3(i+1,1)  = Ns_3_k1(i+1,1) + Ns_3_k2(i+1,1) + Ns_3_k3(i+1,1) + Ns_3_k4(i+1,1);
%     Ns_3_alt(i+1,1)  = Ns_3(i,1)+starttosearch_3(i,1)-findparking_3(i,1);
    
    Nns(i+1,1)   = Nns_1(i+1,1)+Nns_2(i+1,1)+Nns_3(i+1,1);
    Np(i+1,1)    = Np_2(i+1,1)+Np_3(i+1,1);
    Ns(i+1,1)    = Ns_2(i+1,1)+Ns_3(i+1,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Garage implementation:

    Ndgp_k1(i+1,1) = Ndgp_k1(i,1) + n_vehicles_ns_dgp_VOT(i,1) + n_vehicles_s_dgp_VOT(i,1) - n_enter_garage(i,1) - n_dgp_searching(i,1);
    Ndgp_k2(i+1,1) = Ndgp_k2(i,1) + n_vehicles_ns_dgp_VOT(i,2) + n_vehicles_s_dgp_VOT(i,2) - n_enter_garage(i,2) - n_dgp_searching(i,2);
    Ndgp_k3(i+1,1) = Ndgp_k3(i,1) + n_vehicles_ns_dgp_VOT(i,3) + n_vehicles_s_dgp_VOT(i,3) - n_enter_garage(i,3) - n_dgp_searching(i,3);
    Ndgp_k4(i+1,1) = Ndgp_k4(i,1) + n_vehicles_ns_dgp_VOT(i,4) + n_vehicles_s_dgp_VOT(i,4) - n_enter_garage(i,4) - n_dgp_searching(i,4);
    Ndgp(i+1,1) = Ndgp_k1(i+1,1) + Ndgp_k2(i+1,1) + Ndgp_k3(i+1,1) + Ndgp_k4(i+1,1);
    
    Ngp(i+1,1) = Ngp(i,1) + sum(n_enter_garage(i,:),2) - ndepart_garage(i,1);
    total_garage_capacity(i+1,1) = h2_getGlobal_parking_garage_capacity;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%     density(i+1,1)=(Nns(i+1,1)+Ns(i+1,1))/Lnetwork;
    density(i+1,1)=(Nns(i+1,1)+Ns(i+1,1)+Ndgp(i+1,1))/Lnetwork;
    if density(i+1,1)<=kc
        speed(i+1,1)=v;
    elseif density(i+1,1)<=kj
        speed(i+1,1)=Qmax/(kc-kj)*(1-kj/density(i+1,1));
    else
        speed(i+1,1)=0;
    end
    availableparking_2(i+1,1) = A_2 - Np_2(i+1,1);
    availableparking_3(i+1,1) = A_3 - Np_3(i+1,1);
%     availableparking_gp(i+1,1) = initial_garage_capacity - Ngp(i+1,1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine bias for maximum driven distance (for each vehicle)
        bias(i+1,1) = (speed(i+1,1)*t)/2 - speed(i+1,1)*t;
        driven_distance_per_time_slice(i+1,1) = speed(i+1,1)*t*(Nns(i+1,1) + Ns(i+1,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%this part is about the number of vehicles experiencing each event----------------------------------
%1st throughtraffic i.e., leavethearea_1(i+1,1)
         dist_now = flipud(cumsum(flipud(speed(1:i,1))))*t;    
         dist_last= dist_now-(speed(i)*t);
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);% UNIFORM DISTRIBUTION: 
         % in case of FIXED VALUE of z: prob = (dist_now>=z).*(dist_last<z);
         % in case of GAMMA DISTRIBUTION: prob = gamcdf(dist_now,2,2.5)-gamcdf(dist_last,2,2.5);
     leavethearea_1(i+1,1) = dot(enterthearea_1(1:i,1),prob); %throughtraffic
         clear prob;
%2nd garage
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
     starttosearch_2(i+1,1)= dot(enterthearea_2(1:i,1), prob);
         
     clear prob
         if Ns_2(i+1)<availableparking_2(i+1)
             findparking_2(i+1,1)  = Ns_2(i+1);
         else
             findparking_2(i+1,1)  = availableparking_2(i+1);
         end
         clear prob;
         prob =fliplr(diff(gamcdf((0:i)*60*t,parameter(1),parameter(2))));
         departparking_2(i+1,1)= dot(findparking_2(1:i,1), prob(1:end));clear prob;
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         leavethearea_2(i+1,1) = dot(departparking_2(1:i,1), prob);clear prob;
%3rd on-street
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         starttosearch_3(i+1,1)= dot(enterthearea_3(1:i,1), prob);
%          clear prob;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         [~,n_demand_VOT]= c1_input_n_demand(i+1);   %| n_demand
         
         % demand matrix for VOT origins:
         matrix_demand_VOT(i+1,:) = n_demand_VOT'.*0.77;   %| /ns for VOT, (i x V)-matrix
         
         for j= 1:V
             matrix_ns_s_VOT(i+1,j) = dot(matrix_demand_VOT(1:i,j), prob);
         end
         clear prob;
         
         %         [~,n_demand_VOT]= c1_input_n_demand(i+1);   %| n_demand
         %
         %         % demand matrix for VOT origins:
         %         matrix_demand_VOT(i+1,:) = n_demand_VOT';   %| /ns for VOT, (i x V)-matrix
         %
         %         [~,n_ns_s_VOT] = d1_transitions_n_ns_s(enterthearea(1:i,1),speed(1:i,1),t,matrix_demand_VOT);      %| ns/s
         %
         %         % ns/s matrix for VOT origins:
         %         matrix_ns_s_VOT(i+1,:) = n_ns_s_VOT';       %| ns/s for VOT, (i x V)-matrix
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         [~, parking_probability(i+1,:), guessed_price_vector(i,1), tau(i,1), VOT_k1(i,1), VOT_k2(i,1), VOT_k3(i,1), VOT_k4(i,1), ...
             E_p_vot(i,1), travel_distance_cost(i,1), penalty_distance(i,:), ...
             total_costs(i,:), parking_pricing(i,1), delta_searching_veh_available_parking_spots(i,1),...
             delta_searching_veh_available_parking_spots_i_plus_one(i,1)] = ...
             d2_transitions_n_s_p(availableparking_3(i+1,1),Ns_3(i+1,1),L,speed(i+1,1),t,matrix_ns_s_VOT,availableparking_3(1:i,1),...
             Ns_3(1:i,1),speed(1:i,1),0,A_3, parking_pricing, delta_searching_veh_available_parking_spots,...
             Ns_3_k1(i+1,1),Ns_3_k2(i+1,1),Ns_3_k3(i+1,1),Ns_3_k4(i+1,1),decideforparking_3(1:i,1),starttosearch_3(1:i,1));              %| s/p
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         findparking_3(i+1,1) = d2_transitions_n_s_p_new(availableparking_3(i+1,1),Ns_3(i+1,1),speed(i+1,1));
%          probability(i+1,1) = d2_transitions_n_s_p_new(availableparking_3(i+1,1),Ns_3(i+1,1),speed(i+1,1))/Ns_3(i+1,1);
         
         if Ns_3(i+1,1) == 0
             decideforparking_3(i+1,1) = ...
                 findparking_3(i+1,1) * 1/4 * parking_probability(i+1,1) + ...
                 findparking_3(i+1,1) * 1/4 * parking_probability(i+1,2) + ...
                 findparking_3(i+1,1) * 1/4 * parking_probability(i+1,3) + ...
                 findparking_3(i+1,1) * 1/4 * parking_probability(i+1,4);

             decideforparking_3_VOT(i+1,1) = findparking_3(i+1,1) * 1/4 * parking_probability(i+1,1);
             decideforparking_3_VOT(i+1,2) = findparking_3(i+1,1) * 1/4 * parking_probability(i+1,2);
             decideforparking_3_VOT(i+1,3) = findparking_3(i+1,1) * 1/4 * parking_probability(i+1,3);
             decideforparking_3_VOT(i+1,4) = findparking_3(i+1,1) * 1/4 * parking_probability(i+1,4);
         else
             decideforparking_3(i+1,1) = ...
                 findparking_3(i+1,1) * Ns_3_k1(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,1) + ...
                 findparking_3(i+1,1) * Ns_3_k2(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,2) + ...
                 findparking_3(i+1,1) * Ns_3_k3(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,3) + ...
                 findparking_3(i+1,1) * Ns_3_k4(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,4);
             
             decideforparking_3_VOT(i+1,1) = findparking_3(i+1,1) * Ns_3_k1(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,1);
             decideforparking_3_VOT(i+1,2) = findparking_3(i+1,1) * Ns_3_k2(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,2);
             decideforparking_3_VOT(i+1,3) = findparking_3(i+1,1) * Ns_3_k3(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,3);
             decideforparking_3_VOT(i+1,4) = findparking_3(i+1,1) * Ns_3_k4(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,4);
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Determine total bias for all vehicles
         
%          bias_tot(i+1,1) = findparking_3(i+1,1) * bias(i+1,1);
         bias_tot(i+1,1) = decideforparking_3(i+1,1) * bias(i+1,1);
         E_driven_distance(i+1,1) = driven_distance_per_time_slice(i+1,1) + bias_tot(i+1,1);
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Garage implementation:
         
         number_of_parking_garages = c19_input_number_of_parking_garages;
         avg_driving_to_next_garage = L / (2 * number_of_parking_garages);
         
% %          if switch_on_parking_garage_information == 1
% %              garage_information_penalty =  ((initial_garage_capacity/total_garage_capacity(i,1))^parking_garage_information_level) - 1;
% %                          
% %          end
         
         for k = 1:4
             
             [ACT(i+1)] = d15_avg_cruising_time(starttosearch_3(1:i,1), decideforparking_3(1:i,1), n_vehicles_s_dgp_VOT(1:i,:), n_dgp_searching(1:i,:), t);
             
             [walking_distance_gp, walking_distance_op] = d14_walking_distance(L);
             
             C_op(i+1,k) = parking_pricing(i,1) + c9_input_price_per_distance*speed(i+1)*ACT(i+1) + input_value_of_time(k)*ACT(i+1) ...
                 + input_value_of_time(k)*walking_distance_op/c22_input_walking_speed;
             
             avg_driving_to_next_garage = d16_avg_driving_to_next_garage(L);
             C_drive(i+1,k) = c9_input_price_per_distance*avg_driving_to_next_garage + input_value_of_time(k)*avg_driving_to_next_garage/speed(i+1);
             
             C_walk(i+1,k) = input_value_of_time(k)*walking_distance_gp/c22_input_walking_speed;
                       
         end
             
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Off-street parking pricing:
         % garage_parking_pricing(i+1,1) = d17_gp_pricing(parking_pricing(i,1), speed(i+1), ACT(i+1), input_value_of_time, C_drive(i+1,:), C_walk(i+1,:), availableparking_3(i+1,1), garage_parking_pricing(i,1));
         garage_parking_pricing(i+1,1) = h2_getGlobal_initial_garage_parking_pricing;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         for k = 1:4
             
             C_gp(i+1,k) = garage_parking_pricing(i+1,1) + C_drive(i+1,k) + C_walk(i+1,k);
             
             if switch_on_parking_garage_information == 1  
                 C_gp(i+1,k) = C_gp(i+1,k) * A_3 / (total_garage_capacity(i,1) + A_3);                
                 C_op(i+1,k) = C_op(i+1,k) * total_garage_capacity(i,1) / (total_garage_capacity(i,1) + A_3);                
             else
                 C_gp(i+1,k) = C_gp(i+1,k) * A_3 / (initial_garage_capacity + A_3);                 
                 C_op(i+1,k) = C_op(i+1,k) * initial_garage_capacity / (initial_garage_capacity + A_3);
             end
             
             share_off_street_parking_decision(i+1,k) = exp( (C_op(i+1,k) - C_gp(i+1,k)) / min(C_op(i+1,k),C_gp(i+1,k)) ) /...
                 (exp( (C_op(i+1,k) - C_gp(i+1,k)) / min(C_op(i+1,k),C_gp(i+1,k)) ) + 1);
         end
         
         [n_vehicles_ns_dgp(i+1,1), n_vehicles_ns_dgp_VOT(i+1,:)] = d9_transitions_n_ns_dgp(matrix_ns_s_VOT(i+1,:), share_off_street_parking_decision(i+1,:));
         starttosearch_3(i+1,1) = starttosearch_3(i+1,1) - n_vehicles_ns_dgp(i+1,1);
         for j= 1:V
             matrix_ns_s_VOT(i+1,j) = matrix_ns_s_VOT(i+1,j) - n_vehicles_ns_dgp_VOT(i+1,j);
         end
         
         [n_vehicles_s_dgp_VOT(i+1,:)] = d10_transitions_n_s_dgp(Ns_3_k1(i+1,1), Ns_3_k2(i+1,1), Ns_3_k3(i+1,1), Ns_3_k4(i+1,1), share_off_street_parking_decision(i+1,:), decideforparking_3_VOT(i+1,:));
         
         [n_enter_garage(i+1,:), n_enter_garage_total, n_want_to_garage_initially] = d11_transitions_n_dgp_gp(speed, t, avg_driving_to_next_garage, n_vehicles_ns_dgp_VOT, ...
             n_vehicles_s_dgp_VOT, C_op(i+1,:), C_gp(i+1,:), C_drive(i+1,:), Ndgp_k1(i+1,1), Ndgp_k2(i+1,1), Ndgp_k3(i+1,1), Ndgp_k4(i+1,1));
         [n_dgp_searching(i+1,:)] = d12_transitions_n_dgp_s(n_enter_garage_total, n_want_to_garage_initially, Ndgp_k1(i+1,1), Ndgp_k2(i+1,1), Ndgp_k3(i+1,1), Ndgp_k4(i+1,1), A_3, initial_garage_capacity);        
         
         ndepart_garage(i+1,1)= startinginarea_garage*gampdf(i*60*t,parameter(1),parameter(2));
         [ndepart_gar] = d13_transitions_n_gp_ns(n_enter_garage, t);
         ndepart_garage(i+1,1)= ndepart_garage(i+1,1) + ndepart_gar;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         prob =fliplr(diff(gamcdf((0:i)*60*t,parameter(1),parameter(2))));
         departparking_3(i+1,1)= startinginarea*gampdf(i*60*t,parameter(1),parameter(2));
         ndepart = d3_transitions_n_p_ns(decideforparking_3,t);
%          departparking_3(i+1,1)= departparking_3(i+1,1) + dot(decideforparking_3(1:i,1), prob(1:end));
         departparking_3(i+1,1)= departparking_3(i+1,1) + ndepart;
         clear prob;
         
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         leavethearea_3(i+1,1) = 0;
         for j = 1:i
             leavethearea_3(i+1,1) = leavethearea_3(i+1,1) + (departparking_3(j,1) + ndepart_garage(j,1)) * prob(j);
         end
%          leavethearea_3(i+1,1) = dot(departparking_3(1:i,1),prob);
         clear prob;
         %total
         starttosearch(i+1,1)  = starttosearch_2(i+1,1)+starttosearch_3(i+1,1);
         findparking(i+1,1)    = findparking_2(i+1,1)+findparking_3(i+1,1);
         decideforparking(i+1,1) = decideforparking_3(i+1,1);
         departparking(i+1,1)  = departparking_2(i+1,1)+departparking_3(i+1,1);
         leavethearea(i+1,1)   = leavethearea_1(i+1,1)+leavethearea_2(i+1,1)+leavethearea_3(i+1,1);
         %hard part of modelling each individual transition event.
%%
    clc
    disp(['== Coding line: ' num2str(i+1) ' ==']);
    i = i + 1;
    if i == 1000
        i = 1000;
    end
    
    if speed(i,1)==0
        display('== G R I D L O C K ==')
        break
    end
end
%%      
if speed(i,1)~=0
    disp('===  FINAL RESULT  ===');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine total revenue:
[total_revenue, total_on_street_revenue, total_garage_revenue] = d18_total_revenue(departparking, ndepart_garage);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine total bias for all vehicles over all time slices
sum_bias_tot = sum(bias_tot(1:1440,1));
sum_E_driven_distance = sum(E_driven_distance(1:1440,1));
total_error = sum_bias_tot/sum_E_driven_distance;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
occupancy_2(:,1)=1-availableparking_2(:,1)/A_2;
occupancy_3(:,1)=1-availableparking_3(:,1)/A_3;
% occupancy_gp(:,1)=1-availableparking_gp(:,1)/initial_garage_capacity;
occupancy_gp(:,1)=1-total_garage_capacity(:,1)/initial_garage_capacity;
% occupancy_all(:,1) = 1 - (availableparking_3(:,1) + availableparking_gp(:,1)) / (A_3 + initial_garage_capacity);
occupancy_all(:,1) = 1 - (availableparking_3(:,1) + total_garage_capacity(:,1)) / (A_3 + initial_garage_capacity);

%cumulative
for j=1:i
    cum_enterthearea_1(j) = sum(enterthearea_1(1:j));
    cum_leavethearea_1(j) = sum(leavethearea_1(1:j));
    
    cum_enterthearea_2(j) = sum(enterthearea_2(1:j));
    cum_starttosearch_2(j)= sum(starttosearch_2(1:j));
    cum_findparking_2(j)  = sum(findparking_2(1:j));
    cum_departparking_2(j)= sum(departparking_2(1:j));
    cum_leavethearea_2(j) = sum(leavethearea_2(1:j));
    
    cum_enterthearea_3(j)  = startinginarea + startinginarea_garage + sum(enterthearea_3(1:j));
    cum_starttosearch_3(j) = startinginarea + startinginarea_garage + sum(starttosearch_3(1:j));
    cum_findparking_3(j)   = startinginarea + startinginarea_garage + sum(findparking_3(1:j));
    cum_decideforparking_3(j) = startinginarea + startinginarea_garage + sum(decideforparking_3(1:j));
    cum_departparking_3(j) = sum(departparking_3(1:j));
    cum_leavethearea_3(j)  = sum(leavethearea_3(1:j));
    cum_vehicles_ns_dgp(j) = startinginarea + startinginarea_garage + sum(n_vehicles_ns_dgp(1:j,1));
    cum_vehicles_s_dgp(j)  = startinginarea + startinginarea_garage + sum(sum(n_vehicles_s_dgp_VOT(1:j,:)));
    cum_dgp_searching(j)   = startinginarea + startinginarea_garage + sum(sum(n_dgp_searching(1:j,:)));
    cum_enter_garage(j)    = startinginarea + startinginarea_garage + sum(sum(n_enter_garage(1:j,:)));
    cum_depart_garage(j)   = sum(ndepart_garage(1:j,1));

    cum_enterthearea(j)   = startinginarea + startinginarea_garage + sum(enterthearea(1:j));
    cum_starttosearch(j)  = startinginarea + startinginarea_garage + sum(starttosearch(1:j));
    cum_findparking(j)    = startinginarea + startinginarea_garage + sum(findparking(1:j));
    cum_decideforparking(j) = startinginarea + startinginarea_garage + sum(decideforparking(1:j));
    cum_departparking(j)  = sum(departparking(1:j));
    cum_leavethearea(j)   = sum(leavethearea(1:j));
end
% % total_searchers_onstreet = sum(enterthearea_3(1:T,1))
% % % total_searchtime_onstreet=sum(cum_starttosearch_3-cum_findparking_3)*t % unit in hours.
% % total_searchtime_onstreet = sum((cum_starttosearch_3 + cum_dgp_searching)-(cum_decideforparking_3 + cum_vehicles_s_dgp))*t % unit in hours.
% % % total_searchtime_onstreet_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_findparking_3(1+11*60:16*60))*t % unit in hours.
% % total_searchtime_onstreet_worsthours = sum(cum_starttosearch_3(1+11*60:16*60) + cum_dgp_searching(1+11*60:16*60))*t - ...
% %     sum(cum_decideforparking_3(1+11*60:16*60) + cum_vehicles_s_dgp(1+11*60:16*60))*t % unit in hours.
% % total_searchtime_onstreet_deduct1min = total_searchtime_onstreet - total_searchers_onstreet*t % unit in hours.
% % % total_searchtime_onstreet_deduct1min_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_findparking_3(1+11*60:16*60))*t-total_searchers_onstreet*t % unit in hours.
% % total_searchtime_onstreet_deduct1min_worsthours = sum(cum_starttosearch_3(1+11*60:16*60) + cum_dgp_searching(1+11*60:16*60))*t ...
% %     - sum(cum_decideforparking_3(1+11*60:16*60) + cum_vehicles_s_dgp(1+11*60:16*60))*t-total_searchers_onstreet*t % unit in hours.
% % 
% % % average_searchtime_onstreet=total_searchtime_onstreet/sum(findparking_3)*60 % average searching time per person unit in minutes
% % average_searchtime_onstreet = total_searchtime_onstreet/...
% %     (sum(decideforparking_3) + sum(sum(n_vehicles_s_dgp_VOT(:,:),1),2))*60 % average searching time per person unit in minutes
% % % total_searchtime_garage = sum(cum_starttosearch_2-cum_findparking_2)*t % unit in hours.
%%
time=1:1440;

% % 1ST PLOT ----------------------------------------------------------------input (parking demand enterthearea)
% % figure
% % plot(time/60,cum_enterthearea_3)
% % xlabel('Hour of the day')
% % ylabel('Cumulative number of vehicles entering the area')
% % axis([0 24 0 2500])
% % ax = gca;
% % ax.XTick = 0:1:24;
% % ax.YTick = 0:500:2500;
% 
% % 2ND PLOT ----------------------------------------------------------------cumulative plot for first 3 curves 
% % figure
% % a=2.2;
% % a1=1+10*60;
% % a2=1+18*60;
% % background=1:a:1+1439*a;
% % plot(time(a1:a2)/60,cum_enterthearea_3(a1:a2)-background(a1:a2));hold on;
% % plot(time(a1:a2)/60,cum_starttosearch_3(a1:a2)-background(a1:a2),'--');
% % plot(time(a1:a2)/60,cum_findparking_3(a1:a2)-background(a1:a2),'k:');
% % % plot(time(a1:a2)/60,cum_departparking_3(a1:a2)-background(a1:a2));
% % % plot(time(a1:a2)/60,cum_leavethearea_3(a1:a2)-background(a1:a2));
% % hold off;   
% % axis([10 18 -380 -280])
% % ax = gca;
% % ax.XTick = 10:1:18;
% % ax.YTick = -380:50:-280;
% % xlabel('Hour of the day')
% % ylabel('Transformed cumulative number of vehicles (only for parkers)')
% % legend('enter area','start search','found parking')
% 
% %
% -------------------------------------------------------------------occupancy on-street
figure
plot(time/60,occupancy_3,'LineWidth',2)
xlabel('Hour of the day')
axis([0 24 0 1])
ylabel('On-street parking occupancy')
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
set(gca,'FontSize',24)

% -------------------------------------------------------------------occupancy garage
figure
plot(time/60,occupancy_gp,'LineWidth',2)
xlabel('Hour of the day')
axis([0 24 0 1])
ylabel('Garage parking occupancy')
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
set(gca,'FontSize',24)

% -------------------------------------------------------------------occupancy ALL
figure
plot(time/60,occupancy_all,'LineWidth',2)
xlabel('Hour of the day')
axis([0 24 0 1])
ylabel('Parking occupancy')
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
set(gca,'FontSize',24)

% -------------------------------------------------------------------------number of on-street parking searchers
figure
plot(time/60,Ns,'--.','LineWidth', 2)
hold on
plot(time/60,availableparking_3, 'LineWidth', 2)
xlabel('Hour of the day')
ylabel('Number')
legend('Searching vehicles on-street','Number of available on-street parking')
axis([0 24 0 inf])
ax = gca;
ax.XTick = 0:1:24;
set(gca,'FontSize',24)
% ax.YTick = 0:5:40;

% -------------------------------------------------------------------------number of garage parking searchers
figure
plot(time/60,Ndgp,'--.','LineWidth', 2)
hold on
plot(time/60,total_garage_capacity, 'LineWidth', 2)
xlabel('Hour of the day')
ylabel('Number')
legend('Vehicles driving to garage parking','Number of available garage parking')
axis([0 24 0 inf])
ax = gca;
ax.XTick = 0:1:24;
set(gca,'FontSize',24)
% ax.YTick = 0:5:40;


% PLOT --------------------------------------------------------------------share of traffic
figure
area(time/60, [movmean(Ns_3./(Ns+Nns+Ndgp),10) movmean(Ndgp./(Ns+Nns+Ndgp),10) movmean(Nns./(Ns+Nns+Ndgp),10)], 'LineStyle','none')
xlabel('Hour of the day') % x-axis label
ylabel('Share of Traffic') % y-axis label
legend('Share of traffic searching for on-street parking','Share of traffic driving to garage parking','Share of traffic non-searching')
axis([0 24 0 1])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
set(gca,'FontSize',24)

% PLOT
% ----------------------------------------------------------------s/dgp and dgp/s
figure
plot(time/60,movmean(sum(n_vehicles_ns_dgp_VOT,2),10),'g-.','LineWidth',3)
hold on
plot(time/60,movmean(sum(n_vehicles_s_dgp_VOT,2),10),'r-','LineWidth',3)
hold on
n_dgp_s_total = movmean(sum(n_dgp_searching,2),10);
for j = 1:1:564
    n_dgp_s_total(j) = 0;
end
test = [0.0319571747338214;0.159290231207619;0.244142412180423;0.311546118187965;0.382567850253119;0.502243109148401;0.657291684754052;0.827067652521479;0.991342955696998;1.15408383139618;1.27801165408423;1.27651681841215;1.30204539679821;1.32495533959009;1.36254470906936;1.40081220663529;1.43319998437154;1.46189454879937;1.48861896501626;1.52992797220949];
for j = 565:1:584
    n_dgp_s_total(j) = test(j-564);
end

plot(time/60,n_dgp_s_total,'b:','LineWidth',3)
xlabel('Hour of the day')
axis([0 24 0 inf])
ylabel('Number of vehicles in transition event')
hold off
legend('Go to garage parking (ns/dgp)','Switch to garage parking (s/dgp)','Not access garage parking (dgp/s)')
ax = gca;
ax.XTick = 0:1:24;
set(gca,'FontSize',24)
 
% % -------------------------------------------------------------------------probability  
% figure
% plot(time/60,findparking_3)
% xlabel('Hour of the day')
% ylabel('Probability of finding parking')
% axis([0 24 0 1])
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:0.1:1;
% yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
% percentsy = repmat('%', length(yticks),1);  %  equal to the size
% yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
% set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis


% -------------------------------------------------------------------------average cruising time
avg_cruising_time(:,1) = zeros(1440,1);
for i=1:1:1440 
    avg_cruising_time(i,1) = ACT(i) / t;
end    

figure
plot(time/60,avg_cruising_time, 'b', 'LineWidth', 2)
xlabel('Hour of the day')
ylabel('Average cruising time (minutes)')
axis([8 15 0 8])
ax = gca;
ax.XTick = 8:1:15;
ax.YTick = 0:1:8;
set(gca,'FontSize',24)

% -------------------------------------------------------------------------parking choice
share_off_street_parking_decision_average(:,1) = zeros(1440,1);
for i=1:1:1440 
    share_off_street_parking_decision_average(i,1) = sum(share_off_street_parking_decision(i,:),2) / 4;
end    

figure
plot(time/60,share_off_street_parking_decision_average, 'b', 'LineWidth', 2)
xlabel('Hour of the day')
ylabel('Parking choice for garage over on-street parking')
axis([8 15 0.39 0.48])
ax = gca;
ax.XTick = 8:1:15;
ax.YTick = 0.39:0.05:0.48;
set(gca,'FontSize',24)

% -----------------------------------------------------------parking choice incl garage usage plot
% figure
% plot(time_test,share_off_street_parking_decision_average_test, 'b-', 'LineWidth', 2)
% hold on
% 
% plot(time_test,share_off_street_parking_decision_garage_usage, 'r--', 'LineWidth', 2)
% hold on
% 
% xlabel('Hour of the day')
% ylabel('Parking choice for garage over on-street parking')
% legend('No availability of garage usage information for all drivers','Availability of garage usage information for all drivers')
% axis([0 24 0 0.5])
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:0.05:0.5;
% set(gca,'FontSize',24)

% --------------------------------------------------------------------parking choice test
time_test(:,1) = zeros(288,1);
share_off_street_parking_decision_average_test(:,1) = zeros(288,1);
for i=5:5:1440
    time_test(i/5,1) = time(i)/60;
    share_off_street_parking_decision_average_test(i/5,1) = sum(share_off_street_parking_decision(i,:),2) / 4;
end

figure
plot(time_test,share_off_street_parking_decision_average_test, 'b', 'LineWidth', 2)
xlabel('Hour of the day')
ylabel('Parking choice for garage over on-street parking')
axis([0 24 0 0.5])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.05:0.5;
set(gca,'FontSize',24)


% figure
% plot(time/60,share_off_street_parking_decision(:,1), 'LineWidth', 2)
% xlabel('Hour of the day')
% ylabel('1Parking choice for garage over on-street parking')
% axis([0 24 0 inf])
% ax = gca;
% ax.XTick = 0:1:24;
% % ax.YTick = 0:1:15;
% set(gca,'FontSize',24)
% 
% figure
% plot(time/60,share_off_street_parking_decision(:,2), 'LineWidth', 2)
% xlabel('Hour of the day')
% ylabel('2Parking choice for garage over on-street parking')
% axis([0 24 0 inf])
% ax = gca;
% ax.XTick = 0:1:24;
% % ax.YTick = 0:1:15;
% set(gca,'FontSize',24)
% 
% figure
% plot(time/60,share_off_street_parking_decision(:,3), 'LineWidth', 2)
% xlabel('Hour of the day')
% ylabel('3Parking choice for garage over on-street parking')
% axis([0 24 0 inf])
% ax = gca;
% ax.XTick = 0:1:24;
% % ax.YTick = 0:1:15;
% set(gca,'FontSize',24)
% 
% figure
% plot(time/60,share_off_street_parking_decision(:,4), 'LineWidth', 2)
% xlabel('Hour of the day')
% ylabel('4Parking choice for garage over on-street parking')
% axis([0 24 0 inf])
% ax = gca;
% ax.XTick = 0:1:24;
% % ax.YTick = 0:1:15;
% set(gca,'FontSize',24)

% % -------------------------------------------------------------------------average searching time
% for i=1:1:size(cum_findparking_3,2) %time 
%     % the ith row corresponds to the round(cum_starttosearch_3(i),3)
% for x=i:1:size(cum_findparking_3,2) %time
%        if round(cum_findparking_3(x))==round(cum_starttosearch_3(i))
%             searchingtime(i)=x-i;
%             break
%         end
% end
% end
% figure
% plot(time/60,1./probability-1)
% xlabel('Hour of the day')
% ylabel('Average search time before parking (minutes)')
% axis([0 24 0 15])
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:1:15;
% % -------------------------------------------------------------------------------------------------------
% 
% plot(time/60,availableparking_3,'LineWidth',2)
% legend();
% % axis([0 24 0 539])
% xlabel('Time (hr)')
% ylabel('Number of available parking spaces')
% ax = gca;
% ax.XTick = 0:1:24;
% 
% figure
% plot(time/60,Ns_3,'LineWidth',2)
% % axis([0 24 0 35])
% xlabel('Time (hr)')
% ylabel('Number of searching vehicles')
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:5:35;
% 
% 
% figure
% plot(time/60,cum_enterthearea_3,'DisplayName','cum_enterthearea_3','LineWidth',2);hold on;
% plot(time/60,cum_starttosearch_3,'DisplayName','cum_starttosearch_3','LineWidth',2);
% plot(time/60,cum_findparking_3,'DisplayName','cum_findparking_3','LineWidth',2);
% plot(time/60,cum_departparking_3,'DisplayName','cum_departparking_3','LineWidth',2);
% plot(time/60,cum_leavethearea_3,'DisplayName','cum_leavethearea_3','LineWidth',2);hold off;
% xlabel('Time (hr)')
% ylabel('Cumulative number of vehicles')
% legend('cum enterthearea','cum starttosearch','cum accessparking','cum departparking','cum leavethearea')
% axis([0 24 0 2500])
% ax = gca;
% ax.XTick = 0:1:24;
% 
% % figure
% % plot(occupancy_3,1./probability-1)
% % xlabel('Parking occupancy') % x-axis label
% % ylabel('Average searching time') % y-axis label
% % xlim([0.8 1])
% % 
% % figure
% % plot(Ns,probability)
% % xlabel('# searching vehicles') % x-axis label
% % ylabel('Probability') % y-axis label
% % 
% % figure
% % plot(availableparking_3,probability)
% % xlabel('# available parking spaces') % x-axis label
% % ylabel('Probability of finding parking') % y-axis label
% % xlim([0 30])
% 
% 
% % figure
% % plot(occupancy_3,(1./probability-1).*Ns/60)
% % xlabel('Parking occupancy') % x-axis label
% % ylabel('Total searching time (hrs)') % y-axis label
% % xlim([0.9 1])
% 
% figure
% plot(occupancy_3(1:750),Ns(1:750))
% hold on
% plot(occupancy_3(751:1440),Ns(751:1440))
% xlabel('Parking occupancy') % x-axis label
% ylabel('# searching vehicles') % y-axis label
% xlim([0 1])
% legend('00:00-12:00 loading','12:00-24:00 unloading')
% 
% figure
% plot(availableparking_3(1:750),Ns(1:750))
% hold on
% plot(availableparking_3(751:1440),Ns(751:1440))
% xlabel('# available parking space') % x-axis label
% ylabel('# searching vehicles') % y-axis label
% legend('00:00-12:00 loading','12:00-24:00 unloading')
% xlim([0 30])
% 
% figure 
% plot(time/60,cum_enterthearea,'DisplayName','cum_enterthearea');hold on;
% plot(time/60,cum_starttosearch,'DisplayName','cum_starttosearch');
% plot(time/60,cum_findparking,'DisplayName','cum_findparking');
% plot(time/60,cum_departparking,'DisplayName','cum_findparking');
% plot(time/60,cum_leavethearea,'DisplayName','cum_leavethearea');
% hold off;
% axis([0 24 0 3000])
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:500:3000;
% xlabel('Time (hr)')
% ylabel('Cumulative number of vehicles')
% legend('cum enterthearea','cum starttosearch','cum findparking','cum departparking','cum leavethearea')
% 
% % % a=1.3;background=1:a:1+1439*a;
% % figure 
% % plot(time/60,cum_enterthearea_3-background,'DisplayName','cum_enterthearea');hold on;
% % plot(time/60,cum_starttosearch_3-background,'DisplayName','cum_starttosearch');
% % plot(time/60,cum_findparking_3-background,'DisplayName','cum_findparking');
% % plot(time/60,cum_departparking_3-background,'DisplayName','cum_findparking');
% % plot(time/60,cum_leavethearea_3-background,'DisplayName','cum_leavethearea');
% % hold off;
% % xlabel('Time (hr)')
% % ylabel('Cumulative number of vehicles (only for parkers)')
% % legend('cum enterthearea','cum starttosearch','cum findparking','cum departparking','cum leavethearea')
% 
% figure 
% plot(time/60,Nns,'DisplayName','Non-searching vehicles');hold on;
% plot(time/60,Ns,'DisplayName','Searching vehicles');
% hold off;
% axis([0 24 0 40])
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:10:40;
% xlabel('Time (hr)')
% ylabel('Number of vehicles')
% legend('Non-searching vehicles','Searching vehicles')
 
[average_searching_time, total_searching_time, average_non_searching_time_parkers, total_non_searching_time_parkers,...
    average_searching_distance, total_searching_distance, average_non_searching_distance_parkers, total_non_searching_distance_parkers, ...
    average_deciding_gp_time, total_deciding_gp_time, average_deciding_gp_distance, total_deciding_gp_distance,...
    avg_total_time, tot_time, avg_total_distance, tot_distance, final_cum_revenue, final_cum_revenue_garage] = ...
    b_new_outputs_from_the_matrix(Nns, Nns_1, Nns_2, Nns_3, ...
    Ns, Np, Ndgp, Ngp, speed, cum_enterthearea, cum_enterthearea_1, cum_enterthearea_2, cum_enterthearea_3, cum_departparking_3, cum_vehicles_ns_dgp, ...
    cum_vehicles_s_dgp, findparking_3, findparking_3, parking_pricing, n_enter_garage, garage_parking_pricing, startinginarea);

save('results_and_validation_all_simulation_outputs.mat')

end