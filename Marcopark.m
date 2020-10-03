function [average_searching_time, average_searching_distance, final_cum_revenue, occupancy_fuel, occupancy_electric, occupancy_3_all, ACT, optimal_target_occupancy_rate_all_vehicles, optimal_target_occupancy_rate_fuel_vehicles, optimal_target_occupancy_rate_electric_vehicles, social_impacts] = ...
    Marcopark(policy, demand_proportion, supply_proportion, demand_options, supply_options, parking_duration_a_shape, parking_duration_b_scale)

% clear all
% clc

global R t L A_2 A_3 v Lnetwork kc kj Qmax;
    R=0.1; % unit: km
    t=1/60; % STEP SLICE UNIT: hours
    L=7.7; %KM
    A_2=0; % unit: parking spaces in garages
    A_3=539; % unit: parking spaces on-street
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
enterarea(:,1) = demand'*(1 - demand_proportion); %fuel vehicles
enterarea(:,2) = demand'*demand_proportion; %electric vehicles

if demand_options == 0
    % initial demand
else
    % demand increase or decrease:
    enterarea = demand_options .* enterarea; %more or less demand
end

if supply_options == 0
    % initial supply
else
    % supply increase or decrease:
    A_3 = supply_options .* A_3; %more or less supply
end

load('parameter')%parking duration parameters
load('startinginarea') % number of available parking spaces at the beginning

T=size(enterarea,1);
% throughtraffic_ratio(1:size(smallcircle_in,1),1)=0; % tt is the percentage of throughtraffic corresponding to each time slice.
% garage_ratio(1:size(smallcircle_in,1),1)=0; % gg is the percentage of garage corresponding to each time slice.
%%
%ENTER THE AREA--------------------------------
%1st category: throughtraffic
%2nd category: vehicles to garage
%3rd category: vehicles to on-street parking
    enterthearea_1(1:T,1:2)=enterarea*0.23; 
    enterthearea_2(1:T,1:2)=0; 
    enterthearea_3(1:T,1:2)=enterarea*0.77; 
%total:
    enterthearea(:,1) = enterarea(:,1) + enterarea(:,2);
%%
    density = 0;
    speed = v;
    availableparking_2 = A_2;
    availableparking_3(:,1) = (A_3 - startinginarea)*(1 - supply_proportion); %fuel vehicle parking spaces
    availableparking_3(:,2) = (A_3 - startinginarea)*supply_proportion; %electric vehicle parking spaces
    availableparking_3_total(:,1) = availableparking_3(:,1) + availableparking_3(:,2);
    
    Nns_1(1,1)=0;   Nns_2(1,1)=0;    Nns_3(1,1)=0;
    Nns_1(1,2)=0;   Nns_2(1,2)=0;    Nns_3(1,2)=0;
    
    Np_2(1,1)=0;    Np_3(1,1)=startinginarea*(1 - supply_proportion);
    Np_2(1,2)=0;    Np_3(1,2)=startinginarea*supply_proportion;
    
    Ns_2(1,1)=0;    Ns_3(1,1)=0;
    Ns_2(1,2)=0;    Ns_3(1,2)=0;
    
    Ns_3_k1(1,1)=0;     Ns_3_k1(1,2)=0;
    Ns_3_k2(1,1)=0;     Ns_3_k2(1,2)=0;
    Ns_3_k3(1,1)=0;     Ns_3_k3(1,2)=0;
    Ns_3_k4(1,1)=0;     Ns_3_k4(1,2)=0;
    
    Nns(1,1)= Nns_1(1,1)+Nns_2(1,1)+Nns_3(1,1)+Nns_1(1,2)+Nns_2(1,2)+Nns_3(1,2); %non-searching can be all, throughtraffic, garage and on-street parkers
    Np(1,1)= Np_2(1,1)+Np_3(1,1)+Np_2(1,2)+Np_3(1,2); %parking can be both garage and on-street parkers
    Ns(1,1)= Ns_2(1,1)+Ns_3(1,1)+Ns_2(1,2)+Ns_3(1,2); %searching can only be on-street parkers
%%
%the first row of every event--------------------------------
    leavethearea_1(1,1) =0;     leavethearea_1(1,2) =0;    
    starttosearch_2(1,1)=0;     starttosearch_2(1,2)=0;
    findparking_2(1,1)  =0;     findparking_2(1,2)  =0; 
    departparking_2(1,1)=0;     departparking_2(1,2)=0;
    leavethearea_2(1,1) =0;     leavethearea_2(1,2) =0;

    starttosearch_3(1,1)=0;     starttosearch_3(1,2)=0;
    findparking_3(1,1)  =0;     findparking_3(1,2)  =0;
    findparking_3_total(1,1)  =0;
    decideforparking_3(1,1)=0;  decideforparking_3(1,2)=0;
    departparking_3(1,1)=0;     departparking_3(1,2)=0;
    leavethearea_3(1,1) =0;     leavethearea_3(1,2) =0;
    
    starttosearch(1,1)  = starttosearch_2(1,1)+starttosearch_3(1,1)+starttosearch_2(1,2)+starttosearch_3(1,2);
    findparking(1,1)    = findparking_2(1,1)+findparking_3(1,1)+findparking_2(1,2)+findparking_3(1,2);
    decideforparking(1,1) = decideforparking_3(1,1)+decideforparking_3(1,2);
    departparking(1,1)  = departparking_2(1,1)+departparking_3(1,1)+departparking_2(1,2)+departparking_3(1,2);
    leavethearea(1,1)   = leavethearea_1(1,1)+ leavethearea_2(1,1)+leavethearea_3(1,1)+leavethearea_1(1,2)+ leavethearea_2(1,2)+leavethearea_3(1,2);
    
    parking_probability(1,:) = ones(1,4);
    
    
    [~, parkingduration_expectation_fuel, parkingduration_expectation_electric,~,~] = c3_input_parkingduration(1,1,parking_duration_a_shape,parking_duration_b_scale); %the numbers (1,1,...) are irrelevant since only expecation value is needed
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1_setGlobal_switch_on_damp_exp_coef(0);

% set global variable initial parking pricing to 2.25:
h1_setGlobal_initial_parking_pricing(2.25);

% set global variable for maximum price increase per time step to 0.1:  
h1_setGlobal_max_parking_price_increase(0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,n_demand_VOT_ini]= c1_input_n_demand(1,demand_options);%/
n_demand_VOT(:,1) = n_demand_VOT_ini(:,1)'.*(1 - demand_proportion); %fuel vehicles
n_demand_VOT(:,2) = n_demand_VOT_ini(:,1)'.*demand_proportion; %electric vehicles
        
 % demand matrix for VOT origins:
matrix_demand_VOT(1,:,1) = n_demand_VOT(:,1)'.*0.77; 
matrix_demand_VOT(1,:,2) = n_demand_VOT(:,2)'.*0.77;
        
% initial price:
parking_pricing(1,1) = c8_input_parking_price;
 
parking_revenue(1,1) = sum(decideforparking_3(1,1)*parking_pricing(1,1))*parkingduration_expectation_fuel/60;
parking_revenue(1,2) = sum(decideforparking_3(1,2)*parking_pricing(1,1))*parkingduration_expectation_electric/60;
 
% initial delta searching:
parking_pricing_switched_on = c13_input_switch_on_parking_pricing;
if parking_pricing_switched_on == 1
    delta_searching_veh_available_parking_spots(1,1) = Ns(1,1)/availableparking_3_total(1,1);
    delta_searching_veh_available_parking_spots_i_plus_one(1,1) = Ns(1,1)/availableparking_3_total(1,1);
elseif parking_pricing_switched_on == 3
    delta_searching_veh_available_parking_spots(1,1) = 1/availableparking_3_total(1,1);
    delta_searching_veh_available_parking_spots_i_plus_one(1,1) = 1/availableparking_3_total(1,1);
else
    delta_searching_veh_available_parking_spots(1,1) = Ns(1,1)/availableparking_3_total(1,1);
    delta_searching_veh_available_parking_spots_i_plus_one(1,1) = Ns(1,1)/availableparking_3_total(1,1);
end

% ns/s matrix for VOT origins:
V = size(n_demand_VOT,1);
matrix_ns_s_VOT(1,:,1:3) = zeros(1,V,3); 

% initial bias:
bias = zeros(1440,1);
bias(1,1) = (speed*t)/2 - speed*t;
bias_tot(1,1) = (findparking_3(1,1) + findparking_3(1,2)) * bias(1,1);
driven_distance_per_time_slice(1,1) = speed*t*(Nns(1,1) + Ns(1,1));
E_driven_distance(1,1) = driven_distance_per_time_slice(1,1) - bias_tot(1,1);
ACT(:) = zeros(1,1440);
ACT(1) = 1/60;

optimal_target_occupancy_rate_all_vehicles = 0;
optimal_target_occupancy_rate_fuel_vehicles = 0;
optimal_target_occupancy_rate_electric_vehicles = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%
i=1;
while i<1440
    %this part is about the number of vehicles in each state
    Nns_1(i+1,1) = Nns_1(i,1)+enterthearea_1(i,1)-leavethearea_1(i,1); % only throughtraffic
    Nns_1(i+1,2) = Nns_1(i,2)+enterthearea_1(i,2)-leavethearea_1(i,2); % only throughtraffic
    Nns_2(i+1,1) = Nns_2(i,1)+enterthearea_2(i,1)-starttosearch_2(i,1)+departparking_2(i,1)-leavethearea_2(i,1); % only garage parked vehicles
    Nns_2(i+1,2) = Nns_2(i,2)+enterthearea_2(i,2)-starttosearch_2(i,2)+departparking_2(i,2)-leavethearea_2(i,2); % only garage parked vehicles
    Nns_3(i+1,1) = Nns_3(i,1)+enterthearea_3(i,1)-starttosearch_3(i,1)+departparking_3(i,1)-leavethearea_3(i,1); % only on-street parked vehicles
    Nns_3(i+1,2) = Nns_3(i,2)+enterthearea_3(i,2)-starttosearch_3(i,2)+departparking_3(i,2)-leavethearea_3(i,2); % only on-street parked vehicles
    
    Np_2(i+1,1)  = Np_2(i,1)+findparking_2(i,1)-departparking_2(i,1);% only garage parked vehicles
    Np_2(i+1,2)  = Np_2(i,2)+findparking_2(i,2)-departparking_2(i,2);% only garage parked vehicles
%     Np_3(i+1,1)  = Np_3(i,1)+findparking_3(i,1)-departparking_3(i,1);% only on-street parked vehicles
    Np_3(i+1,1)  = Np_3(i,1)+decideforparking_3(i,1)-departparking_3(i,1);% only on-street parked vehicles
    Np_3(i+1,2)  = Np_3(i,2)+decideforparking_3(i,2)-departparking_3(i,2);% only on-street parked vehicles
    
    Ns_2(i+1,1)  = Ns_2(i,1)+starttosearch_2(i,1)-findparking_2(i,1);
    Ns_2(i+1,2)  = Ns_2(i,2)+starttosearch_2(i,2)-findparking_2(i,2);
    
    if Ns_3(i,1) == 0
        Ns_3_k1(i+1,1)  = Ns_3_k1(i,1) + matrix_ns_s_VOT(i,1,1) - findparking_3(i,1) * 1/4 * parking_probability(i,1);
        Ns_3_k2(i+1,1)  = Ns_3_k2(i,1) + matrix_ns_s_VOT(i,2,1) - findparking_3(i,1) * 1/4 * parking_probability(i,2);
        Ns_3_k3(i+1,1)  = Ns_3_k3(i,1) + matrix_ns_s_VOT(i,3,1) - findparking_3(i,1) * 1/4 * parking_probability(i,3);
        Ns_3_k4(i+1,1)  = Ns_3_k4(i,1) + matrix_ns_s_VOT(i,4,1) - findparking_3(i,1) * 1/4 * parking_probability(i,4);
    else
        Ns_3_k1(i+1,1)  = Ns_3_k1(i,1) + matrix_ns_s_VOT(i,1,1) - findparking_3(i,1) * Ns_3_k1(i,1) ./ Ns_3(i,1) * parking_probability(i,1);
        Ns_3_k2(i+1,1)  = Ns_3_k2(i,1) + matrix_ns_s_VOT(i,2,1) - findparking_3(i,1) * Ns_3_k2(i,1) ./ Ns_3(i,1) * parking_probability(i,2);
        Ns_3_k3(i+1,1)  = Ns_3_k3(i,1) + matrix_ns_s_VOT(i,3,1) - findparking_3(i,1) * Ns_3_k3(i,1) ./ Ns_3(i,1) * parking_probability(i,3);
        Ns_3_k4(i+1,1)  = Ns_3_k4(i,1) + matrix_ns_s_VOT(i,4,1) - findparking_3(i,1) * Ns_3_k4(i,1) ./ Ns_3(i,1) * parking_probability(i,4);        
    end
    
    if Ns_3(i,2) == 0
        Ns_3_k1(i+1,2)  = Ns_3_k1(i,2) + matrix_ns_s_VOT(i,1,2) - findparking_3(i,2) * 1/4 * parking_probability(i,1);
        Ns_3_k2(i+1,2)  = Ns_3_k2(i,2) + matrix_ns_s_VOT(i,2,2) - findparking_3(i,2) * 1/4 * parking_probability(i,2);
        Ns_3_k3(i+1,2)  = Ns_3_k3(i,2) + matrix_ns_s_VOT(i,3,2) - findparking_3(i,2) * 1/4 * parking_probability(i,3);
        Ns_3_k4(i+1,2)  = Ns_3_k4(i,2) + matrix_ns_s_VOT(i,4,2) - findparking_3(i,2) * 1/4 * parking_probability(i,4);
    else
        Ns_3_k1(i+1,2)  = Ns_3_k1(i,2) + matrix_ns_s_VOT(i,1,2) - findparking_3(i,2) * Ns_3_k1(i,2) ./ Ns_3(i,2) * parking_probability(i,1);
        Ns_3_k2(i+1,2)  = Ns_3_k2(i,2) + matrix_ns_s_VOT(i,2,2) - findparking_3(i,2) * Ns_3_k2(i,2) ./ Ns_3(i,2) * parking_probability(i,2);
        Ns_3_k3(i+1,2)  = Ns_3_k3(i,2) + matrix_ns_s_VOT(i,3,2) - findparking_3(i,2) * Ns_3_k3(i,2) ./ Ns_3(i,2) * parking_probability(i,3);
        Ns_3_k4(i+1,2)  = Ns_3_k4(i,2) + matrix_ns_s_VOT(i,4,2) - findparking_3(i,2) * Ns_3_k4(i,2) ./ Ns_3(i,2) * parking_probability(i,4);        
    end

    Ns_3(i+1,1)  = Ns_3_k1(i+1,1) + Ns_3_k2(i+1,1) + Ns_3_k3(i+1,1) + Ns_3_k4(i+1,1);
    Ns_3(i+1,2)  = Ns_3_k1(i+1,2) + Ns_3_k2(i+1,2) + Ns_3_k3(i+1,2) + Ns_3_k4(i+1,2);
    
    Nns(i+1,1)   = Nns_1(i+1,1)+Nns_2(i+1,1)+Nns_3(i+1,1)+Nns_1(i+1,2)+Nns_2(i+1,2)+Nns_3(i+1,2);
    Np(i+1,1)    = Np_2(i+1,1)+Np_3(i+1,1)+Np_2(i+1,2)+Np_3(i+1,2);
    Ns(i+1,1)    = Ns_2(i+1,1)+Ns_3(i+1,1)+Ns_2(i+1,2)+Ns_3(i+1,2);
    
    density(i+1,1)=(Nns(i+1,1)+Ns(i+1,1))/Lnetwork;
    if density(i+1,1)<=kc
        speed(i+1,1)=v;
    elseif density(i+1,1)<=kj
        speed(i+1,1)=Qmax/(kc-kj)*(1-kj/density(i+1,1));
    else
        speed(i+1,1)=0;
    end
    availableparking_2(i+1,1)=A_2 - Np_2(i+1,1);   
    availableparking_3(i+1,1)=A_3*(1 - supply_proportion) - Np_3(i+1,1); %fuel vehicle parking spaces   
    availableparking_3(i+1,2)=A_3*supply_proportion - Np_3(i+1,2); %electric vehicle parking spaces 
    
    if policy == 2     
        if (findparking_3(i,2) - departparking_3(i,2)) > availableparking_3(i+1,2)
            availableparking_3(i+1,1) = availableparking_3(i+1,1) - ((findparking_3(i,2) - departparking_3(i,2)) - availableparking_3(i+1,2));
            availableparking_3(i+1,2) = 0;           
        end
    end

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
     leavethearea_1(i+1,2) = dot(enterthearea_1(1:i,2),prob); %throughtraffic
     clear prob;
%2nd garage
     prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
     starttosearch_2(i+1,1)= dot(enterthearea_2(1:i,1), prob);
     starttosearch_2(i+1,2)= dot(enterthearea_2(1:i,2), prob);
         
     clear prob
         if Ns_2(i+1,1)<availableparking_2(i+1)
             findparking_2(i+1,1)  = Ns_2(i+1,1);
         else
             findparking_2(i+1,1)  = availableparking_2(i+1);
         end
         if Ns_2(i+1,2)<availableparking_2(i+1)
             findparking_2(i+1,2)  = Ns_2(i+1,2);
         else
             findparking_2(i+1,2)  = availableparking_2(i+1);
         end
         clear prob;
         prob =fliplr(diff(gamcdf((0:i)*60*t,parameter(1),parameter(2))));
         departparking_2(i+1,1)= dot(findparking_2(1:i,1), prob(1:end));
         departparking_2(i+1,2)= dot(findparking_2(1:i,2), prob(1:end));
         clear prob;
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         leavethearea_2(i+1,1) = dot(departparking_2(1:i,1), prob);
         leavethearea_2(i+1,2) = dot(departparking_2(1:i,2), prob);
         clear prob;
%3rd on-street
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         starttosearch_3(i+1,1)= dot(enterthearea_3(1:i,1), prob);
         starttosearch_3(i+1,2)= dot(enterthearea_3(1:i,2), prob);
%        clear prob;
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,n_demand_VOT_ini]= c1_input_n_demand(i+1,demand_options);   %| n_demand
        n_demand_VOT(:,1) = n_demand_VOT_ini(:,1)'.*(1 - demand_proportion); %fuel vehicles
        n_demand_VOT(:,2) = n_demand_VOT_ini(:,1)'.*demand_proportion; %electric vehicles
        
        % demand matrix for VOT origins:
        matrix_demand_VOT(i+1,:,1) = n_demand_VOT(:,1)'.*0.77;   %| /ns for VOT, (i x V)-matrix 
        matrix_demand_VOT(i+1,:,2) = n_demand_VOT(:,2)'.*0.77;   %| /ns for VOT, (i x V)-matrix 
        
        for j= 1:4
            matrix_ns_s_VOT(i+1,j,1) = dot(matrix_demand_VOT(1:i,j,1), prob);
            matrix_ns_s_VOT(i+1,j,2) = dot(matrix_demand_VOT(1:i,j,2), prob);
        end
        clear prob;
        
%         [~,n_demand_VOT]= c1_input_n_demand(i+1,demand_options);   %| n_demand
%         
%         % demand matrix for VOT origins:
%         matrix_demand_VOT(i+1,:) = n_demand_VOT';   %| /ns for VOT, (i x V)-matrix 
%                 
%         [~,n_ns_s_VOT] = d1_transitions_n_ns_s(enterthearea(1:i,1),speed(1:i,1),t,matrix_demand_VOT);      %| ns/s
%         
%         % ns/s matrix for VOT origins:
%         matrix_ns_s_VOT(i+1,:) = n_ns_s_VOT';       %| ns/s for VOT, (i x V)-matrix

        [~, parking_probability(i+1,:), guessed_price_vector(i,1), tau(i,1), VOT_k1(i,1), VOT_k2(i,1), VOT_k3(i,1), VOT_k4(i,1), ...
            E_p_vot(i,1), travel_distance_cost(i,1), penalty_distance(i,:), ...
            total_costs(i,:), parking_pricing(i,1), delta_searching_veh_available_parking_spots(i,1),...
            delta_searching_veh_available_parking_spots_i_plus_one(i,1)] = ...
         d2_transitions_n_s_p(availableparking_3(i+1,1) + availableparking_3(i+1,2),Ns_3(i+1,1) + Ns_3(i+1,2),L,speed(i+1,1),t,matrix_ns_s_VOT(:,:,1) + matrix_ns_s_VOT(:,:,2),...
            availableparking_3(1:i,1) + availableparking_3(1:i,2),Ns_3(1:i,1) + Ns_3(1:i,2),speed(1:i,1),0,A_3, parking_pricing, delta_searching_veh_available_parking_spots,...
            Ns_3_k1(i+1,1) + Ns_3_k1(i+1,2),Ns_3_k2(i+1,1) + Ns_3_k2(i+1,2),Ns_3_k3(i+1,1) + Ns_3_k3(i+1,2),Ns_3_k4(i+1,1) + Ns_3_k4(i+1,2),decideforparking_3(1:i,1) + decideforparking_3(1:i,2),...
            starttosearch_3(1:i,1) + starttosearch_3(1:i,2));                                                %| s/p
                  
         [ACT(i+1)] = d9_avg_cruising_time(starttosearch_3(1:i,1) + starttosearch_3(1:i,2), decideforparking_3(1:i,1) + decideforparking_3(1:i,2), t);
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         if policy == 1
             findparking_3(i+1,1) = d2_transitions_n_s_p_new(availableparking_3(i+1,1),Ns_3(i+1,1),speed(i+1,1)); %fuel vehicles
             findparking_3(i+1,2) = d2_transitions_n_s_p_new(availableparking_3(i+1,2),Ns_3(i+1,2),speed(i+1,1)); %electric vehicles
         elseif policy == 2
             if availableparking_3(i+1,2) >= 1
                 findparking_3(i+1,1) = d2_transitions_n_s_p_new(availableparking_3(i+1,1),Ns_3(i+1,1),speed(i+1,1)); %fuel vehicles
                 findparking_3(i+1,2) = d2_transitions_n_s_p_new(availableparking_3(i+1,1) + availableparking_3(i+1,2) - findparking_3(i+1,1), Ns_3(i+1,2),speed(i+1,1)); %all spaces with electric vehicles
             else %electric parking is full, and electric vehicles face fuel parking spaces
                 findparking_3(i+1,1) = d2_transitions_n_s_p_new(availableparking_3(i+1,1),Ns_3(i+1,1),speed(i+1,1)); %fuel vehicles
                 findparking_3(i+1,2) = d2_transitions_n_s_p_new(availableparking_3(i+1,1) + availableparking_3(i+1,2) - findparking_3(i+1,1), Ns_3(i+1,2),speed(i+1,1)); %all spaces with electric vehicles
                 if (Ns_3(i+1,1) + Ns_3(i+1,2)) > 0
                     findparking_3(i+1,1) = (findparking_3(i+1,1) + findparking_3(i+1,2)) * Ns_3(i+1,1) / (Ns_3(i+1,1) + Ns_3(i+1,2)); %Update fuel vehicles entering a fuel parking space depending on their searching ratio
                     findparking_3(i+1,2) = (findparking_3(i+1,1) + findparking_3(i+1,2)) * Ns_3(i+1,2) / (Ns_3(i+1,1) + Ns_3(i+1,2)); %Update electric vehicles entering a fuel parking space depending on their searching ratio
                 end
             end
         end
         probability(i+1,1) = d2_transitions_n_s_p_new(availableparking_3(i+1,1) + availableparking_3(i+1,2),Ns_3(i+1,1) + Ns_3(i+1,2),speed(i+1,1))/(Ns_3(i+1,1) + Ns_3(i+1,2));
         
         if Ns_3(i+1,1) == 0
             decideforparking_3(i+1,1) = ...
                 findparking_3(i+1,1) * 1/4 * parking_probability(i+1,1) + ...
                 findparking_3(i+1,1) * 1/4 * parking_probability(i+1,2) + ...
                 findparking_3(i+1,1) * 1/4 * parking_probability(i+1,3) + ...
                 findparking_3(i+1,1) * 1/4 * parking_probability(i+1,4);
             
         else
             decideforparking_3(i+1,1) = ...
                 findparking_3(i+1,1) * Ns_3_k1(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,1) + ...
                 findparking_3(i+1,1) * Ns_3_k2(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,2) + ...
                 findparking_3(i+1,1) * Ns_3_k3(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,3) + ...
                 findparking_3(i+1,1) * Ns_3_k4(i+1,1) ./ Ns_3(i+1,1) * parking_probability(i+1,4);
             
         end
         
         if Ns_3(i+1,2) == 0         
             decideforparking_3(i+1,2) = ...
                 findparking_3(i+1,2) * 1/4 * parking_probability(i+1,1) + ...
                 findparking_3(i+1,2) * 1/4 * parking_probability(i+1,2) + ...
                 findparking_3(i+1,2) * 1/4 * parking_probability(i+1,3) + ...
                 findparking_3(i+1,2) * 1/4 * parking_probability(i+1,4);
             
         else
             decideforparking_3(i+1,2) = ...
                 findparking_3(i+1,2) * Ns_3_k1(i+1,2) ./ Ns_3(i+1,2) * parking_probability(i+1,1) + ...
                 findparking_3(i+1,2) * Ns_3_k2(i+1,2) ./ Ns_3(i+1,2) * parking_probability(i+1,2) + ...
                 findparking_3(i+1,2) * Ns_3_k3(i+1,2) ./ Ns_3(i+1,2) * parking_probability(i+1,3) + ...
                 findparking_3(i+1,2) * Ns_3_k4(i+1,2) ./ Ns_3(i+1,2) * parking_probability(i+1,4);
             
         end
         
         parking_revenue(i+1,1) = sum(decideforparking_3(i+1,1)*parking_pricing(1,1))*parkingduration_expectation_fuel/60;
         parking_revenue(i+1,2) = sum(decideforparking_3(i+1,2)*parking_pricing(1,1))*parkingduration_expectation_electric/60;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         % Determine total bias for all vehicles
         
%          bias_tot(i+1,1) = findparking_3(i+1,1) * bias(i+1,1);
         bias_tot(i+1,1) = (decideforparking_3(i+1,1) + decideforparking_3(i+1,2)) * bias(i+1,1);
         E_driven_distance(i+1,1) = driven_distance_per_time_slice(i+1,1) + bias_tot(i+1,1);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
         prob =fliplr(diff(gamcdf((0:i)*60*t,parameter(1),parameter(2))));
         departparking_3(i+1,1)= startinginarea*(1 - supply_proportion)*gampdf(i*60*t,parameter(1),parameter(2)); %fuel vehicles
         departparking_3(i+1,2)= startinginarea*supply_proportion*gampdf(i*60*t,parameter(1)*parking_duration_a_shape,parameter(2)*parking_duration_b_scale); %electric vehicles
%          departparking_3(i+1,1)= departparking_3(i+1,1)+  dot(findparking_3(1:i,1), prob(1:end));
         departparking_3(i+1,1)= departparking_3(i+1,1)+  dot(decideforparking_3(1:i,1), prob(1:end));
         departparking_3(i+1,2)= departparking_3(i+1,2)+  dot(decideforparking_3(1:i,2), prob(1:end));
         clear prob;    
                  
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         leavethearea_3(i+1,1) = dot(departparking_3(1:i,1),prob);
         leavethearea_3(i+1,2) = dot(departparking_3(1:i,2),prob);
         clear prob;
%total
    starttosearch(i+1,1)    = starttosearch_2(i+1,1)+starttosearch_3(i+1,1)+starttosearch_2(i+1,2)+starttosearch_3(i+1,2);
    findparking(i+1,1)      = findparking_2(i+1,1)+findparking_3(i+1,1)+findparking_2(i+1,2)+findparking_3(i+1,2); 
    decideforparking(i+1,1) = decideforparking_3(i+1,1)+decideforparking_3(i+1,2);
    departparking(i+1,1)    = departparking_2(i+1,1)+departparking_3(i+1,1)+departparking_2(i+1,2)+departparking_3(i+1,2); 
    leavethearea(i+1,1)     = leavethearea_1(i+1,1)+leavethearea_2(i+1,1)+leavethearea_3(i+1,1)+leavethearea_1(i+1,2)+leavethearea_2(i+1,2)+leavethearea_3(i+1,2);
%hard part of modelling each individual transition event.
%%
    clc
    disp(['== Coding line: ' num2str(i+1) ' ==']);
    i = i + 1;
    if i == 1291
       i = 1291; 
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
if speed(i,1)==0
    average_searching_time = 0;
    average_searching_distance = 0;
    final_cum_revenue = 0;
    occupancy_fuel = 0;
    occupancy_electric = 0; 
    occupancy_3_all = 0;
    ACT = 0;
    optimal_target_occupancy_rate_all_vehicles = 0;
    optimal_target_occupancy_rate_fuel_vehicles = 0;
    optimal_target_occupancy_rate_electric_vehicles = 0;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine total bias for all vehicles over all time slices
sum_bias_tot = sum(bias_tot(1:1440,1));
sum_E_driven_distance = sum(E_driven_distance(1:1440,1));
total_error = sum_bias_tot/sum_E_driven_distance;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
% occupancy_2(:,1)=1-availableparking_2(:,1)/A_2;
occupancy_3(:,1)=1-availableparking_3(:,1)/(A_3*(1 - supply_proportion)); %fuel vehicles
occupancy_fuel(:,1) = occupancy_3(:,1);
occupancy_3(:,2)=1-availableparking_3(:,2)/(A_3*supply_proportion); %electric vehicles
occupancy_electric(:,1) = occupancy_3(:,2);
occupancy_3_all(:,1)=1-(availableparking_3(:,1) + availableparking_3(:,2))/A_3; %fuel and electric vehicles
%cumulative
for j=1:i
    cum_enterthearea_1(j,1) = sum(enterthearea_1(1:j,1));
    cum_enterthearea_1(j,2) = sum(enterthearea_1(1:j,2));
    cum_leavethearea_1(j,1) = sum(leavethearea_1(1:j,1));
    cum_leavethearea_1(j,2) = sum(leavethearea_1(1:j,2));
    
    cum_enterthearea_2(j,1) = sum(enterthearea_2(1:j,1));
    cum_enterthearea_2(j,2) = sum(enterthearea_2(1:j,2));
    cum_starttosearch_2(j,1)= sum(starttosearch_2(1:j,1));
    cum_starttosearch_2(j,2)= sum(starttosearch_2(1:j,2));
    cum_findparking_2(j,1)  = sum(findparking_2(1:j,1));
    cum_findparking_2(j,2)  = sum(findparking_2(1:j,2));
    cum_departparking_2(j,1)= sum(departparking_2(1:j,1));
    cum_departparking_2(j,2)= sum(departparking_2(1:j,2));
    cum_leavethearea_2(j,1) = sum(leavethearea_2(1:j,1));
    cum_leavethearea_2(j,2) = sum(leavethearea_2(1:j,2));
    
    cum_enterthearea_3(j,1) = startinginarea*(1 - supply_proportion)+sum(enterthearea_3(1:j,1));
    cum_enterthearea_3(j,2) = startinginarea*supply_proportion+sum(enterthearea_3(1:j,2));
    cum_starttosearch_3(j,1)= startinginarea*(1 - supply_proportion)+sum(starttosearch_3(1:j,1));
    cum_starttosearch_3(j,2)= startinginarea*supply_proportion+sum(starttosearch_3(1:j,2));
    cum_findparking_3(j,1)  = startinginarea*(1 - supply_proportion)+sum(findparking_3(1:j,1));
    cum_findparking_3(j,2)  = startinginarea*supply_proportion+sum(findparking_3(1:j,2));
    cum_decideforparking_3(j,1) = startinginarea*(1 - supply_proportion)+sum(decideforparking_3(1:j,1));
    cum_decideforparking_3(j,2) = startinginarea*supply_proportion+sum(decideforparking_3(1:j,2));
    cum_departparking_3(j,1)= sum(departparking_3(1:j,1));
    cum_departparking_3(j,2)= sum(departparking_3(1:j,2));
    cum_leavethearea_3(j,1) = sum(leavethearea_3(1:j,1));
    cum_leavethearea_3(j,2) = sum(leavethearea_3(1:j,2));

    cum_enterthearea(j,1)   = startinginarea+sum(enterthearea(1:j,1));
    cum_starttosearch(j,1)  = startinginarea+sum(starttosearch(1:j,1));
    cum_findparking(j,1)    = startinginarea+sum(findparking(1:j,1));
    cum_decideforparking(j,1) = startinginarea+sum(decideforparking(1:j,1));
    cum_departparking(j,1)  = sum(departparking(1:j,1));
    cum_leavethearea(j,1)   = sum(leavethearea(1:j,1));
end
total_searchers_onstreet=sum(enterthearea_3(1:T,1)+enterthearea_3(1:T,2))
% total_searchtime_onstreet=sum(cum_starttosearch_3-cum_findparking_3)*t % unit in hours.
total_searchtime_onstreet=sum(cum_starttosearch_3(:,1) + cum_starttosearch_3(:,2) - (cum_decideforparking_3(:,1) + cum_decideforparking_3(:,2)))*t % unit in hours.
% total_searchtime_onstreet_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_findparking_3(1+11*60:16*60))*t % unit in hours.
total_searchtime_onstreet_worsthours=sum(cum_starttosearch_3(1+11*60:16*60,1) + cum_starttosearch_3(1+11*60:16*60,2))*t - sum(cum_decideforparking_3(1+11*60:16*60,1) + cum_decideforparking_3(1+11*60:16*60,2))*t % unit in hours.
total_searchtime_onstreet_deduct1min=total_searchtime_onstreet-total_searchers_onstreet*t % unit in hours.
% total_searchtime_onstreet_deduct1min_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_findparking_3(1+11*60:16*60))*t-total_searchers_onstreet*t % unit in hours.
total_searchtime_onstreet_deduct1min_worsthours=sum(cum_starttosearch_3(1+11*60:16*60,1) + cum_starttosearch_3(1+11*60:16*60,2))*t - ...
    sum(cum_decideforparking_3(1+11*60:16*60,1) + cum_decideforparking_3(1+11*60:16*60,2))*t - total_searchers_onstreet*t % unit in hours.

% average_searchtime_onstreet=total_searchtime_onstreet/sum(findparking_3)*60 % average searching time per person unit in minutes
average_searchtime_onstreet=total_searchtime_onstreet/sum(decideforparking_3(:,1) + decideforparking_3(:,2))*60 % average searching time per person unit in minutes
total_searchtime_garage=sum(cum_starttosearch_2(:,1) + cum_starttosearch_2(:,2) - (cum_findparking_2(:,1) + cum_findparking_2(:,2)))*t % unit in hours.
%%
time=1:1440;

%%
[average_searching_time, average_searching_distance, final_cum_revenue] = ...
    b_new_outputs_from_the_matrix(Nns, Ns, Np, Nns_1(:,1) + Nns_1(:,2), Nns_3(:,1) + Nns_3(:,2), speed, cum_enterthearea, cum_enterthearea_3(:,1) + cum_enterthearea_3(:,2), ...
    decideforparking_3(:,1), decideforparking_3(:,2), parking_pricing, parkingduration_expectation_fuel, parkingduration_expectation_electric);

%% -----------------------------------------------------------------------------------------------------------------------social impacts
walking_distance = 2/3 * (c18_input_network_block_length *(-1/2 + sqrt(1/4 + L/(2*c18_input_network_block_length))));

value_of_time = c7_input_value_of_time*60;
EV_infrastructure_cost= c16_input_EV_infrastructure_cost;
walking_speed = c17_input_walking_speed;

% all_n_sp = sum(sum(findparking_3(:,:),1),2);

% social_impacts = (sum(departparking_3(:,1),1)/(sum(departparking_3(:,1),1) + sum(departparking_3(:,2),1)))*parking_pricing(1,1)*parkingduration_expectation_fuel/60 + ...
%     (sum(departparking_3(:,2),1)/(sum(departparking_3(:,1),1) + sum(departparking_3(:,2),1)))*parking_pricing(1,1)*parkingduration_expectation_electric/60 + ...
%     (sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k1(:,1) + Ns_3_k1(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1))/all_n_sp * value_of_time(1)*(2*walking_distance/walking_speed + average_searching_time/60) + ...
%     (sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k2(:,1) + Ns_3_k2(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1))/all_n_sp * value_of_time(2)*(2*walking_distance/walking_speed + average_searching_time/60) + ...
%     (sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k3(:,1) + Ns_3_k3(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1))/all_n_sp * value_of_time(3)*(2*walking_distance/walking_speed + average_searching_time/60) + ...
%     (sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k4(:,1) + Ns_3_k4(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1))/all_n_sp * value_of_time(4)*(2*walking_distance/walking_speed + average_searching_time/60) + ...
%     EV_infrastructure_cost*(A_3*supply_proportion);

% all_n_sp = sum(sum(findparking_3(:,:),1),2);
%      (sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k1(:,1) + Ns_3_k1(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1))/all_n_sp * value_of_time(1) + ...
%      (sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k2(:,1) + Ns_3_k2(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1))/all_n_sp * value_of_time(2) + ...
%      (sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k3(:,1) + Ns_3_k3(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1))/all_n_sp * value_of_time(3) + ...
%      (sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k4(:,1) + Ns_3_k4(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1))/all_n_sp * value_of_time(4)

social_impacts = sum(departparking_3(:,1),1)*parking_pricing(1,1)*parkingduration_expectation_fuel/60 + ...
    sum(departparking_3(:,2),1)*parking_pricing(1,1)*parkingduration_expectation_electric/60 + ...
    sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k1(:,1) + Ns_3_k1(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1) * value_of_time(1)*(2*walking_distance/walking_speed + average_searching_time/60) + ...
    sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k2(:,1) + Ns_3_k2(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1) * value_of_time(2)*(2*walking_distance/walking_speed + average_searching_time/60) + ...
    sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k3(:,1) + Ns_3_k3(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1) * value_of_time(3)*(2*walking_distance/walking_speed + average_searching_time/60) + ...
    sum((findparking_3(:,1) + findparking_3(:,2)),1) .* sum((Ns_3_k4(:,1) + Ns_3_k4(:,2)),1) ./ sum((Ns_3(:,1) + Ns_3(:,2)),1) * value_of_time(4)*(2*walking_distance/walking_speed + average_searching_time/60) + ...
    EV_infrastructure_cost*(A_3*supply_proportion);

%% 1ST PLOT ----------------------------------------------------------------input (parking demand enterthearea)
figure
plot(time/60,cum_enterthearea_3(:,1) + cum_enterthearea_3(:,2))
xlabel('Hour of the day')
ylabel('Cumulative number of vehicles entering the area')
axis([0 24 0 2500])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:500:2500;

%% 2ND PLOT ----------------------------------------------------------------cumulative plot for first 3 curves 
figure
a=2.2;
a1=1+10*60;
a2=1+18*60;
background=(1:a:1+1439*a)';
plot(time(a1:a2)/60,cum_enterthearea_3(a1:a2,1) + cum_enterthearea_3(a1:a2,2) - background(a1:a2,1));hold on;
plot(time(a1:a2)/60,cum_starttosearch_3(a1:a2,1) + cum_starttosearch_3(a1:a2,2) - background(a1:a2,1),'--');
plot(time(a1:a2)/60,cum_findparking_3(a1:a2,1) + cum_findparking_3(a1:a2,2) - background(a1:a2,1),'k:');
% plot(time(a1:a2)/60,cum_departparking_3(a1:a2)-background(a1:a2));
% plot(time(a1:a2)/60,cum_leavethearea_3(a1:a2)-background(a1:a2));
hold off;   
axis([10 18 -380 -280])
ax = gca;
ax.XTick = 10:1:18;
% ax.YTick = -380:50:-280;
xlabel('Hour of the day')
ylabel('Transformed cumulative number of vehicles (only for parkers)')
legend('enter area','start search','found parking')

%% --------------------------------------------------------------------------occupancy
figure
plot(time/60,occupancy_3_all,'LineWidth',2)
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

%% -------------------------------------------------------------------------number of searchers
figure
plot(time/60,Ns,'--.')
hold on
plot(time/60,availableparking_3(:,1) + availableparking_3(:,2))
xlabel('Hour of the day')
ylabel('Number')
legend('Searching vehicles','Number of available parking')
axis([0 24 0 40])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:5:40;

%% PLOT --------------------------------------------------------------------share of traffic
figure 
plot(time/60,(Ns_3(:,1) + Ns_3(:,2))./(Ns + Nns));
xlabel('Hour of the day') % x-axis label
ylabel('Share of Traffic Searching for Parking') % y-axis label
axis([0 24 0 1])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis


%% -------------------------------------------------------------------------probability  
figure
plot(time/60,probability)
xlabel('Hour of the day')
ylabel('Probability of finding parking')
axis([0 24 0 1])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

%% -------------------------------------------------------------------------average searching time
for i=1:1:size(cum_findparking_3(:,1) + cum_findparking_3(:,2),2) %time 
    % the ith row corresponds to the round(cum_starttosearch_3(i),3)
for x=i:1:size(cum_findparking_3(:,1) + cum_findparking_3(:,2),2) %time
       if round(cum_findparking_3(x,1) + cum_findparking_3(x,2))==round(cum_starttosearch_3(i,1) + cum_starttosearch_3(i,2))
            searchingtime(i)=x-i;
            break
        end
end
end
figure
% plot(time/60,1./probability-1)
plot(time/60,60.*ACT)
xlabel('Hour of the day')
ylabel('Average search time before parking (minutes)')
axis([0 24 0 15])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:1:15;

%% ---------------------------------------------------------------optimal parking occupancy rate (moving average)

[occupancy_3_all_sorted,sorting_order] = sort(occupancy_3_all);
ACT_sorted = ACT(sorting_order);
ACT_movmean = movmean(ACT_sorted,11);
for r = 1:1400
%     if 60*ACT_movmean(1,r) == 60*t
    if round(60*ACT_movmean(1,r),1) == 60*t
        optimal_target_occupancy_rate_all_vehicles = occupancy_3_all_sorted(r,1);
    end
end
optimal_target_occupancy_rate_all_vehicles


figure
hold on
% ACT = ACT - 1/60;
% scatter(occupancy_3_all,60.*ACT,'bx','LineWidth',2)
[occupancy_3_all_sorted,sorting_order] = sort(occupancy_3_all);
ACT_sorted = ACT(sorting_order);
ACT_movmean = movmean(ACT_sorted,11);

% for i = 1:size(ACT_movmean,2)
%     if ACT_movmean(1,i) < 0
%         ACT_movmean(1,i) = 0;
%     end
% end

plot(occupancy_3_all_sorted,60.*ACT_movmean,'r-','LineWidth',3)
% hold on
% plot(linspace(0,1,11),ones(1,11)*average_searching_time,'b--','LineWidth',3)
xlabel('Parking occupancy rate')
axis([0 1 0 inf])
ylabel('Average searching time (minutes)')
ax = gca;
ax.XTick = 0:0.1:1;
% ax.YTick = 0:1:24;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% leg = legend('Data points for all time slices','Moving average','Average searching time for parking over one working day');
% leg = legend('Moving average','Average searching time for parking over one working day');
leg = legend('Moving average');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
hold off

%% ---------------------------------------------------------------optimal parking occupancy rate (moving average): fuel vehicles

figure
hold on
scatter(occupancy_fuel,60.*ACT,'bx','LineWidth',2)
[occupancy_fuel_sorted,sorting_order_fuel] = sort(occupancy_fuel);
ACT_sorted_fuel = ACT(sorting_order_fuel);
ACT_movmean_fuel = movmean(ACT_sorted_fuel,11);

for r = 1:1400
%     if 60*ACT_movmean_fuel(1,r) == 60*t
    if round(60*ACT_movmean_fuel(1,r),1) == 60*t
        optimal_target_occupancy_rate_fuel_vehicles = occupancy_fuel_sorted(r,1);
    end
end
optimal_target_occupancy_rate_fuel_vehicles

plot(occupancy_fuel_sorted,60.*ACT_movmean_fuel,'r-','LineWidth',3)
hold on
plot(linspace(0,1,11),ones(1,11)*average_searching_time,'g--','LineWidth',3)
xlabel('Parking occupancy rate of parking spaces for fuel vehicles')
axis([0 1 0 inf])
ylabel('Average searching time (minutes)')
ax = gca;
ax.XTick = 0:0.1:1;
% ax.YTick = 0:1:24;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
leg = legend('Data points for all time slices','Moving average','Average searching time for parking over one working day');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
hold off

%% ---------------------------------------------------------------optimal parking occupancy rate (moving average): electric vehicles

figure
hold on
scatter(occupancy_electric,60.*ACT,'bx','LineWidth',2)
[occupancy_electric_sorted,sorting_order_electric] = sort(occupancy_electric);
ACT_sorted_electric = ACT(sorting_order_electric);
ACT_movmean_electric = movmean(ACT_sorted_electric,11);

for r = 1:1400
%     if 60*ACT_movmean_electric(1,r) == 60*t
    if round(60*ACT_movmean_electric(1,r),1) == 60*t
        optimal_target_occupancy_rate_electric_vehicles = occupancy_electric_sorted(r,1);
    end
end
optimal_target_occupancy_rate_electric_vehicles

plot(occupancy_electric_sorted,60.*ACT_movmean_electric,'r-','LineWidth',3)
hold on
plot(linspace(0,1,11),ones(1,11)*average_searching_time,'g--','LineWidth',3)
xlabel('Parking occupancy rate of parking spaces for electric vehicles')
axis([0 1 0 inf])
ylabel('Average searching time (minutes)')
ax = gca;
ax.XTick = 0:0.1:1;
% ax.YTick = 0:1:24;
xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
percentsx = repmat('%', length(xticks),1);  %  equal to the size
xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
leg = legend('Data points for all time slices','Moving average','Average searching time for parking over one working day');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
hold off

%% ---------------------------------------------------------------optimal parking occupancy rate
% 
% figure
% hold on
% scatter(occupancy_3_all,60.*ACT,'bx','LineWidth',2)
% hold on
% plot(linspace(0,1,11),ones(1,11)*average_searching_time,'g--','LineWidth',3)
% xlabel('Parking occupancy rate')
% axis([0 1 0 inf])
% ylabel('Average cruising time (minutes)')
% ax = gca;
% ax.XTick = 0:0.1:1;
% % ax.YTick = 0:1:24;
% xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% percentsx = repmat('%', length(xticks),1);  %  equal to the size
% xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% legend('Data points for different parking occupancy rates', 'Average cruising time over one day')
% set(gca,'FontSize',24)
% hold off

% %% ---------------------------------------------------------------optimal parking occupancy rate incl. extrapolation
% 
% figure
% g = fittype('a*exp(b*x)+1');
% f = fit(occupancy_3_all,60.*ACT',g,'StartPoint',[0,0],'Lower',[0,20]);
% % f = fit(occupancy_3_all,60.*ACT',g,'StartPoint',[0,0],'Lower',[0,20]); % 5% Less Demand
% p = plot(f,occupancy_3_all,60.*ACT);
% p(1).Marker = 'x';
% p(2).LineWidth = 3;
% % plot(occupancy_3_all,60.*ACT,'Marker','x','LineWidth',3);
% hold on
% plot(linspace(0,1,11),ones(1,11)*average_searching_time,'g--','LineWidth',3)
% xlabel('Parking occupancy rate')
% % axis([0 1 0 inf])
% axis([0 1 0 10])
% ylabel('Average searching time (minutes)')
% ax = gca;
% ax.XTick = 0:0.1:1;
% % ax.YTick = 0:1:24;
% xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% percentsx = repmat('%', length(xticks),1);  %  equal to the size
% xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% legend('Data Points', 'Fitted Exponential Curve', 'Average searching time over one day')
% set(gca,'FontSize',24)
% hold off

% %% -----------------------------------------------------parking occupancy rate and revenue
% 
% figure
% hold on
% scatter(occupancy_3_all,parking_revenue,'bx','LineWidth',2)
% xlabel('Parking Occupancy Rate')
% axis([0 1 0 inf])
% ylabel('Parking revenue (in CHF)')
% ax = gca;
% ax.XTick = 0:0.1:1;
% % ax.YTick = 0:1:24;
% xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% percentsx = repmat('%', length(xticks),1);  %  equal to the size
% xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% set(gca,'FontSize',24)
% hold off

% %% -----------------------------------optimal parking occupancy rate (probability)
% 
% probability(1,1) = 0.9999;
% for i = 1:size(probability,1)
%     if isnan(probability(i,1)) == 1
%         probability(i,1) = 0.9999;
%     end
% end
% figure
% % scatter(occupancy_3,searchingtime)
% % g = fittype('-a./(x.^b+c)');
% % g = fittype('a*exp(b*x)');
% % f = fit(occupancy_3,1./probability-1,g,'StartPoint',[0,0],'Lower',[0,0]);
% % p = plot(f,occupancy_3,1./probability-1);
% % p(1).Marker = 'x';
% % p(2).LineWidth = 3;
% plot(occupancy_3,1./probability-1,'Marker','x','LineWidth',3);
% hold on
% plot(linspace(0,1,11),ones(1,11)*average_searching_time,'r--','LineWidth',3)
% xlabel('Parking Occupancy Rate')
% axis([0 1 0 inf])
% % axis([0 1 0 12])
% ylabel('Average searching time (minutes)')
% ax = gca;
% ax.XTick = 0:0.1:1;
% % ax.YTick = 0:1:24;
% xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% percentsx = repmat('%', length(xticks),1);  %  equal to the size
% xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% legend('Simulation Data', 'Average searching time over one day')
% set(gca,'FontSize',24)
% hold off

% %% -------------------------------------------- optimal parking occupancy rate (approach ACT)
% 
% probability(1,1) = 0.9999;
% for i = 1:size(probability,1)
%     if isnan(probability(i,1)) == 1
%         probability(i,1) = 0.9999;
%     end
% end
% figure
% % scatter(occupancy_3,searchingtime)
% % g = fittype('-a./(x.^b+c)');
% g = fittype('a*exp(b*x)');
% % f = fit(occupancy_3,60.*ACT',g,'StartPoint',[0,0],'Lower',[0,0]);
% f = fit(occupancy_3,1./probability-1,g,'StartPoint',[0,0],'Lower',[0,0]);
% % f = fit(occupancy_3,searchingtime','exp1');
% % p = plot(f,occupancy_3,60.*ACT);
% p = plot(f,occupancy_3,1./probability-1);
% p(1).Marker = 'x';
% p(2).LineWidth = 3;
% hold on
% plot(linspace(0,1,11),ones(1,11)*average_searching_time,'g--','LineWidth',3)
% xlabel('Parking Occupancy Rate')
% axis([0 1 0 inf])
% % axis([0 1 0 5])
% ylabel('Average searching time (minutes)')
% ax = gca;
% ax.XTick = 0:0.1:1;
% % ax.YTick = 0:1:24;
% xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% percentsx = repmat('%', length(xticks),1);  %  equal to the size
% xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% legend('Data Points', 'Fitted Exponential Curve', 'Average searching time over one day')
% set(gca,'FontSize',24)
% hold off
% f

% %% -------------------------------------------------optimal parking occupancy rate (old approach)
% figure
% % scatter(occupancy_3,searchingtime)
% % g = fittype('-a./(x.^b+c)');
% g = fittype('a*exp(b*x)');
% f = fit(occupancy_3,searchingtime',g,'StartPoint',[1,2],'Lower',[0,32]);
% % f = fit(occupancy_3,searchingtime','exp1');
% p = plot(f,occupancy_3,searchingtime);
% p(1).Marker = 'x';
% p(2).LineWidth = 3;
% hold on
% plot(linspace(0,1,11),ones(1,11)*average_searching_time,'g--','LineWidth',3)
% xlabel('Parking Occupancy Rate')
% axis([0 1 0 inf])
% % axis([0 1 0 5])
% ylabel('Average searching time (minutes)')
% ax = gca;
% ax.XTick = 0:0.1:1;
% % ax.YTick = 0:1:24;
% xticks = [get(gca,'xtick')]'; % There is a transpose operation here.
% percentsx = repmat('%', length(xticks),1);  %  equal to the size
% xticklabel = [num2str(100.*xticks) percentsx]; % concatenates the tick labels 
% set(gca,'xticklabel',xticklabel)% Sets tick labels back on the Axis
% legend('Data Points', 'Fitted Exponential Curve', 'Average searching time over one day')
% set(gca,'FontSize',24)
% hold off
% f

%-------------------------------------------------------------------------------------------------------
% 
% plot(time/60,availableparking_3,'LineWidth',2)
% legend();
% axis([0 24 0 539])
% xlabel('Time (hr)')
% ylabel('Number of available parking spaces')
% ax = gca;
% ax.XTick = 0:1:24;
% 
% figure
% plot(time/60,Ns_3,'LineWidth',2)
% axis([0 24 0 35])
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
% figure
% plot(occupancy_3,1./probability-1)
% xlabel('Parking occupancy') % x-axis label
% ylabel('Average searching time') % y-axis label
% xlim([0.8 1])
% 
% figure
% plot(Ns,probability)
% xlabel('# searching vehicles') % x-axis label
% ylabel('Probability') % y-axis label
% 
% figure
% plot(availableparking_3,probability)
% xlabel('# available parking spaces') % x-axis label
% ylabel('Probability of finding parking') % y-axis label
% xlim([0 30])
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
% % a=1.3;background=1:a:1+1439*a;
% figure 
% plot(time/60,cum_enterthearea_3-background,'DisplayName','cum_enterthearea');hold on;
% plot(time/60,cum_starttosearch_3-background,'DisplayName','cum_starttosearch');
% plot(time/60,cum_findparking_3-background,'DisplayName','cum_findparking');
% plot(time/60,cum_departparking_3-background,'DisplayName','cum_findparking');
% plot(time/60,cum_leavethearea_3-background,'DisplayName','cum_leavethearea');
% hold off;
% xlabel('Time (hr)')
% ylabel('Cumulative number of vehicles (only for parkers)')
% legend('cum enterthearea','cum starttosearch','cum findparking','cum departparking','cum leavethearea')
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

ACT = ACT';

% save ('Scenario_1_Optimal_target_parking_occupancy_for_all_vehicles.mat')
% save ('Scenario_2_1_Optimal_target_parking_occupancy_for_all_vehicles.mat')
% save ('Scenario_2_2_Optimal_target_parking_occupancy_for_all_vehicles.mat')
% save ('Scenario_2_3_Optimal_target_parking_occupancy_for_all_vehicles.mat')
% save ('Scenario_3_Optimal_target_parking_occupancy_for_all_vehicles.mat')
% save ('Scenario_4_1_Optimal_target_parking_occupancy_for_all_vehicles.mat')
% save ('Scenario_4_2_Optimal_target_parking_occupancy_for_all_vehicles.mat')
% save ('Scenario_4_3_Optimal_target_parking_occupancy_for_all_vehicles.mat')

end