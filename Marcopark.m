function [average_searching_time, average_searching_distance, final_cum_revenue, cum_revenue_PR, cum_revenue_congestion_toll, cum_revenue_parking_pricing, share_veh_pr_decision_mean_all_times] = ...
    Marcopark(sensitivity_parameter, sensitivity_percentage)

% clear all% SMALL CIRCLE
% clc

% Please select:
% if move_spaces_to_park_ride equals to 1, then P+R is switched on.
% if move_spaces_to_park_ride equals to 0, then basic scenario.
move_spaces_to_park_ride = 1;

global R t L A_2 A_3 PR v Lnetwork kc kj Qmax;
    R=0.1; % unit: km
    t=1/60; % STEP SLICE UNIT: hours
    L=7.7; %KM
    A_2=0; % unit: parking spaces in garages
    A_3 = 539; % unit: parking spaces on-street
    PR = 0; % unit: P+R parking spaces, !!!INPUT CONGESTION PRICING PAPER!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % supply decrease by 4%:
%     A_2 = A_2 * 0.96;
%     A_3 = A_3 * 0.96;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % supply increase by 4%:
%     A_2 = A_2 * 1.04;
%     A_3 = A_3 * 1.04;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%     v=12.5; % unit: km/hr
    v=27.933; % unit: km/hr %3D-MFD for Zurich
    Lnetwork=15.4; %km
    kc=20; % unit: veh/km
    kj=55; % unit: veh/km
    Qmax= v*kc; %unit: veh/hour
% load('ouput_enter.mat')
load('load_traffic_demand')
[demand,~] = hist(act_start_time_seconds./60, linspace(1,1440,1440));
enterarea = demand';

if move_spaces_to_park_ride == 1
   if sensitivity_parameter == 2
        A_3 = 339 + 339*sensitivity_percentage; % unit: parking spaces inside the area (scenario 1 and 2)
   elseif sensitivity_parameter == 6
        A_3 = 339 - 339*sensitivity_percentage; % unit: parking spaces inside the area (scenario 1 and 2)
   else
        A_3 = 339; % unit: parking spaces inside the area (scenario 1 and 2)       
   end
   if sensitivity_parameter == 3
        PR = 200 + 200*sensitivity_percentage; % unit: P+R parking spaces, !!!INPUT CONGESTION PRICING PAPER!!! (scenario 1 and 2)
   elseif sensitivity_parameter == 6
        PR = 200 + 339*sensitivity_percentage; % unit: P+R parking spaces, !!!INPUT CONGESTION PRICING PAPER!!! (scenario 1 and 2)
   else
        PR = 200; % unit: P+R parking spaces, !!!INPUT CONGESTION PRICING PAPER!!! (scenario 1 and 2)
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % demand decrease by 4%:
% enterarea = enterarea * 0.96;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % demand increase by 4%:
% enterarea = enterarea * 1.04;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('parameter')%parking duration parameters
load('startinginarea') % number of available parking spaces at the beginning

if move_spaces_to_park_ride == 1
   if sensitivity_parameter == 2
        startinginarea = 113 + 113*sensitivity_percentage; % (scenario 1 and 2)
   elseif sensitivity_parameter == 6
        startinginarea = 113 - 113*sensitivity_percentage; % (scenario 1 and 2)
   else
        startinginarea = 113;% (scenario 1 and 2)    
   end
end

T=size(enterarea,1);
% throughtraffic_ratio(1:size(smallcircle_in,1),1)=0; % tt is the percentage of throughtraffic corresponding to each time slice.
% garage_ratio(1:size(smallcircle_in,1),1)=0; % gg is the percentage of garage corresponding to each time slice.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CUSTOM INPUT - Congestion Pricing:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes the parking and P+R capacity real-time
% information, switch it on or off.
% parking and P+R capacity real-time information switched on = 1
% parking and P+R capacity real-time information switched off = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h1_setGlobal_parking_pr_capacity_information(0);
%%% CUSTOM INPUT - Congestion Pricing:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch_on_real_time_information_parking_pr_capacity = c19_input_switch_on_real_time_information_parking_pr_capacity;
%%

%ENTER THE AREA--------------------------------
%1st category: throughtraffic
%2nd category: vehicles to garage
%3rd category: vehicles to on-street parking
    enterthearea_1(1:T,1)=enterarea*0.23; 
%     enterthearea_2(1:T,1)=0; 
    enterthearea_3(1:T,1)=enterarea*0.77;
    enterthearea_3(1,1)=enterthearea_3(1,1)*0.5;
    enter_pr(1:T,1)=enterarea*0.77;
    enter_pr(1,1)=enter_pr(1,1)*0.5;
%total:
    enterthearea(:,1)=enterarea; 
%%
    density_car=0; %3D-MFD for Zurich
    density_PT=0; %3D-MFD for Zurich
    speed=v;  
    availableparking_2=A_2;
    availableparking_3=A_3-startinginarea;
    Npr(1,1) = 0; % unit: P+R parking spaces, !!!INPUT CONGESTION PRICING PAPER!!!
        
    if move_spaces_to_park_ride == 1
        if sensitivity_parameter == 3
            Npr(1,1) = 70 + 70*sensitivity_percentage; % unit: P+R parking spaces, !!!INPUT CONGESTION PRICING PAPER!!! (scenario 1 and 2)
        elseif sensitivity_parameter == 6
            Npr(1,1) = 70 + 113*sensitivity_percentage; % unit: P+R parking spaces, !!!INPUT CONGESTION PRICING PAPER!!! (scenario 1 and 2)   
        else
            Npr(1,1) = 70; % unit: P+R parking spaces, !!!INPUT CONGESTION PRICING PAPER!!! (scenario 1 and 2)
        end
    end

    pr_capacity = PR - Npr(1,1);
    
    average_driving_distance_PT = c18_input_network_block_length*(sqrt(c24_input_number_of_PT_stops)/2)*(-1/2 + sqrt(1/4 + L/(2*c18_input_network_block_length))) + 5*c18_input_network_block_length;
    
    PT_speed(1,1) = 0.116*speed(1,1)+9.574; %3D-MFD for Zurich
    PT_speed_pq(1,1) = 0.116*speed(1,1)+9.574; %3D-MFD for Zurich
%     PT_speed(1,1) = speed(1,1)* ((average_driving_distance_PT/speed(1,1)) / ((average_driving_distance_PT/speed(1,1)) + c25_input_average_dwell_time_PT_stops));
    
    Nns_1(1,1)=0;   Nns_2(1,1)=0;    Nns_3(1,1)=0;
    Np_2(1,1)=0;    Np_3(1,1)=startinginarea;
    Ns_2(1,1)=0;    Ns_3(1,1)=0;
    Ns_3_k1(1,1)=0;
    Ns_3_k2(1,1)=0;
    Ns_3_k3(1,1)=0;
    Ns_3_k4(1,1)=0;
    
    Nns(1,1)= Nns_1(1,1)+Nns_2(1,1)+Nns_3(1,1); %non-searching can be all, throughtraffic, garage and on-street parkers
    Np(1,1)= Np_2(1,1)+Np_3(1,1); %parking can be both garage and on-street parkers
    Ns(1,1)= Ns_2(1,1)+Ns_3(1,1); %searching can only be on-street parkers
%%
%the first row of every event--------------------------------
    leavethearea_1(1,1) =0;

    starttosearch_2(1,1)=0;
    findparking_2(1,1)  =0; 
    departparking_2(1,1)=0; 
%     leavethearea_2(1,1) =0;

    starttosearch_3(1,1)=0; 
    findparking_3(1,1)  =0; 
    decideforparking_3(1,1) =0;
    departparking_1(1,1)=0;
    enter_pr(1,1)=0;
    depart_pr(1,1)=0;
    
    starttosearch(1,1)  = starttosearch_2(1,1)+starttosearch_3(1,1);
    findparking(1,1)    = findparking_2(1,1)+findparking_3(1,1);
    decideforparking(1,1) = decideforparking_3(1,1);
    departparking(1,1)  = departparking_1(1,1) + departparking_2(1,1);
    leavethearea(1,1)   = leavethearea_1(1,1);
    
    parking_probability(1,:) = ones(1,4);
    
    [~, parkingduration_expectation,~,~] = c3_input_parkingduration(1,1); %the numbers (1,1) are irrelevant since only expecation value is needed
%     [~, gp_a2, gp_b2] = c20_input_park_ride_parkingduration(1, 1); %the numbers (1,1) are irrelevant since only expecation value is needed
%     park_ride_parkingduration_expectation = gp_a2 * gp_b2;
%     park_ride_parkingduration_expectation = parkingduration_expectation;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1_setGlobal_switch_on_damp_exp_coef(0);

% set global variable initial parking pricing to 3:
% h1_setGlobal_initial_parking_pricing(3);
% h1_setGlobal_initial_parking_pricing(0);

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

% initial bias:
bias = zeros(1440,1);
bias(1,1) = (speed*t)/2 - speed*t;
bias_tot(1,1) = findparking_3(1,1) * bias(1,1);
driven_distance_per_time_slice(1,1) = speed*t*(Nns(1,1) + Ns(1,1));
E_driven_distance(1,1) = driven_distance_per_time_slice(1,1) - bias_tot(1,1);

% VOT per user group:
if sensitivity_parameter == 1
    input_value_of_time = c7_input_value_of_time + c7_input_value_of_time.*sensitivity_percentage;
else
    input_value_of_time = c7_input_value_of_time;
end

occupancy_PT(1,1)=0;

% Initial cost values:
ACT(:) = zeros(1,1440);
% % ANST(:) = zeros(1,1440);
C_drive(:,:) = zeros(1440,4);
C_time(:,:) = zeros(1440,4);
C_veh(:,:) = zeros(1440,4);
C_veh_APR(:,:) = zeros(1440,4);
C_pr(:,:) = zeros(1440,4);
C_pr_APR(:,:) = zeros(1440,4);
share_veh_pr_decision(:,:) = zeros(1440,4);
delay_by_PT(:,1) = zeros(1440,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%
i=1;
while i<1440
    %this part is about the number of vehicles in each state
    Nns_1(i+1,1) = Nns_1(i,1) + enterthearea_1(i,1) + departparking_1(i,1) - leavethearea_1(i,1); % external destination, !!! Nnse = Nns_1 !!!
%     Nns_2(i+1,1) = Nns_2(i,1) + enterthearea_2(i,1) - starttosearch_2(i,1) + departparking_2(i,1) - leavethearea_2(i,1); % only garage parked vehicles
    Nns_3(i+1,1) = Nns_3(i,1) + enterthearea_3(i,1) - starttosearch_3(i,1); % !!! Nnsi = Nns_3 !!!
    
    Npr(i+1,1)  = Npr(i,1) + enter_pr(i,1) - depart_pr(i,1);% only P+R
    Np_3(i+1,1)  = Np_3(i,1) + decideforparking_3(i,1) - departparking_1(i,1);% only parking, !!! Np = Np_3 !!!
    
%     Ns_2(i+1,1)  = Ns_2(i,1)+starttosearch_2(i,1)-findparking_2(i,1);
    
    if Ns_3(i,1) == 0
        Ns_3_k1(i+1,1)  = Ns_3_k1(i,1) + matrix_ns_s_VOT(i,1) - findparking_3(i,1) * 1/4 * parking_probability(i,1);
        Ns_3_k2(i+1,1)  = Ns_3_k2(i,1) + matrix_ns_s_VOT(i,2) - findparking_3(i,1) * 1/4 * parking_probability(i,2);
        Ns_3_k3(i+1,1)  = Ns_3_k3(i,1) + matrix_ns_s_VOT(i,3) - findparking_3(i,1) * 1/4 * parking_probability(i,3);
        Ns_3_k4(i+1,1)  = Ns_3_k4(i,1) + matrix_ns_s_VOT(i,4) - findparking_3(i,1) * 1/4 * parking_probability(i,4);     
                
    else
        Ns_3_k1(i+1,1)  = Ns_3_k1(i,1) + matrix_ns_s_VOT(i,1) - findparking_3(i,1) * Ns_3_k1(i,1) ./ Ns_3(i,1) * parking_probability(i,1);
        Ns_3_k2(i+1,1)  = Ns_3_k2(i,1) + matrix_ns_s_VOT(i,2) - findparking_3(i,1) * Ns_3_k2(i,1) ./ Ns_3(i,1) * parking_probability(i,2);
        Ns_3_k3(i+1,1)  = Ns_3_k3(i,1) + matrix_ns_s_VOT(i,3) - findparking_3(i,1) * Ns_3_k3(i,1) ./ Ns_3(i,1) * parking_probability(i,3);
        Ns_3_k4(i+1,1)  = Ns_3_k4(i,1) + matrix_ns_s_VOT(i,4) - findparking_3(i,1) * Ns_3_k4(i,1) ./ Ns_3(i,1) * parking_probability(i,4);
                
    end
    
 
    Ns_3(i+1,1)  = Ns_3_k1(i+1,1) + Ns_3_k2(i+1,1) + Ns_3_k3(i+1,1) + Ns_3_k4(i+1,1);
%     Ns_3_alt(i+1,1)  = Ns_3(i,1)+starttosearch_3(i,1)-findparking_3(i,1);
    
    Nns(i+1,1)   = Nns_1(i+1,1)+Nns_3(i+1,1);
    Np(i+1,1)    = Np_3(i+1,1);
    Ns(i+1,1)    = Ns_3(i+1,1);
    
    density_car(i+1,1)=(Nns(i+1,1)+Ns(i+1,1))/Lnetwork; %3D-MFD for Zurich
    PT_speed_pq(i+1,1) = (0.116*v + 0.116*(-0.288)*density_car(i+1,1) + 9.574)/2 + ...
        sqrt(((0.116*v + 0.116*(-0.288)*density_car(i+1,1) + 9.574)/2)^2 + (0.116*(-5.659)*2*average_driving_distance_PT)/(c22_input_average_headway*Lnetwork)); %3D-MFD for Zurich
    density_PT(i+1,1)=2*average_driving_distance_PT/(PT_speed_pq(i+1,1)*c22_input_average_headway*Lnetwork); %3D-MFD for Zurich

    speed(i+1,1) = v + (-0.288)*density_car(i+1,1) + (-5.659)*density_PT(i+1,1); %3D-MFD for Zurich
    PT_speed(i+1,1) = 0.116*speed(i+1,1)+9.574; %3D-MFD for Zurich
    
%     if density_car(i+1,1)<=kc
%         speed(i+1,1)=v;
%     elseif density_car(i+1,1)<=kj
%         speed(i+1,1)=Qmax/(kc-kj)*(1-kj/density_car(i+1,1));
%     else
%         speed(i+1,1)=0;
%     end


%     availableparking_2(i+1,1)=A_2-Np_2(i+1,1);
    availableparking_3(i+1,1)=A_3-Np_3(i+1,1);
    pr_capacity(i+1,1) = PR - Npr(i+1,1);
           
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine bias for maximum driven distance (for each vehicle)
        bias(i+1,1) = (speed(i+1,1)*t)/2 - speed(i+1,1)*t;
        driven_distance_per_time_slice(i+1,1) = speed(i+1,1)*t*(Nns(i+1,1) + Ns(i+1,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%this part is about the number of vehicles experiencing each event----------------------------------
         dist_now = flipud(cumsum(flipud(speed(1:i,1))))*t;    
         dist_last= dist_now-(speed(i)*t);
%        prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);% UNIFORM DISTRIBUTION: 
%        in case of FIXED VALUE of z: prob = (dist_now>=z).*(dist_last<z);
%        in case of GAMMA DISTRIBUTION: prob = gamcdf(dist_now,2,2.5)-gamcdf(dist_last,2,2.5);

%2nd garage
%          prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
%      starttosearch_2(i+1,1)= dot(enterthearea_2(1:i,1), prob);
%          
%      clear prob
%          if Ns_2(i+1)<availableparking_2(i+1)
%              findparking_2(i+1,1)  = Ns_2(i+1);
%          else
%              findparking_2(i+1,1)  = availableparking_2(i+1);
%          end
%          clear prob;
%          prob =fliplr(diff(gamcdf((0:i)*60*t,parameter(1),parameter(2))));
%          departparking_2(i+1,1)= dot(findparking_2(1:i,1), prob(1:end));clear prob;
%          prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
%          leavethearea_2(i+1,1) = dot(departparking_2(1:i,1), prob);clear prob;

%3rd on-street

        [~, parking_probability(i+1,:), guessed_price_vector(i,1), tau(i,1), VOT_k1(i,1), VOT_k2(i,1), VOT_k3(i,1), VOT_k4(i,1), ...
            E_p_vot(i,1), travel_distance_cost(i,1), penalty_distance(i,:), ...
             total_costs(i,:), parking_pricing(i,1), delta_searching_veh_available_parking_spots(i,1),...
             delta_searching_veh_available_parking_spots_i_plus_one(i,1)] = ...
         d2_transitions_n_s_p(availableparking_3(i+1,1),Ns_3(i+1,1),L,speed(i+1,1),t,matrix_ns_s_VOT,availableparking_3(1:i,1),...
             Ns_3(1:i,1),speed(1:i,1),0,A_3, parking_pricing, delta_searching_veh_available_parking_spots,...
             Ns_3_k1(i+1,1),Ns_3_k2(i+1,1),Ns_3_k3(i+1,1),Ns_3_k4(i+1,1),decideforparking_3(1:i,1),starttosearch_3(1:i,1));              %| s/p
          
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Congestion pricing implementation, Start:     

         [ACT(i+1)] = d9_avg_cruising_time(starttosearch_3(1:i,1), decideforparking_3(1:i,1), t);
         average_walking_distance = d10_walking_distance(L);
         
         average_driving_distance_PT = c18_input_network_block_length*(sqrt(c24_input_number_of_PT_stops)/2)*(-1/2 + sqrt(1/4 + L/(2*c18_input_network_block_length))) + 5*c18_input_network_block_length;
         
%          PT_speed(i+1,1) = speed(i+1,1)* ((average_driving_distance_PT/speed(i+1,1)) / ...
%              ((average_driving_distance_PT/speed(i+1,1)) + c25_input_average_dwell_time_PT_stops));

         average_PT_walking_distance = ((2*c18_input_network_block_length)/(3 * sqrt(pi * c24_input_number_of_PT_stops))) * (-1/2 + sqrt(1/4 + L/(2*c18_input_network_block_length)));
         
         term_enter_pr = 0;
         if i+1-(c22_input_average_headway/t) >= 1
             for j = (i+1-(c22_input_average_headway/t)):i
                term_enter_pr = term_enter_pr + enter_pr(j,1);
             end
             occupancy_PT(i+1,1) = c26_input_occupancy_rate_passengers_per_vehicle * term_enter_pr;
         else %take first enter_pr(j,k) terms that exist
             for j = 1:i
                term_enter_pr = term_enter_pr + enter_pr(j,1);
             end
             occupancy_PT(i+1,1) = c26_input_occupancy_rate_passengers_per_vehicle * term_enter_pr;  
         end
         
         for k = 1:4
                          
            C_drive(i+1,k) = c9_input_price_per_distance*speed(i+1,1)*ACT(i+1);
            C_time(i+1,k) = input_value_of_time(k)*((0.4/speed(i+1,1)) + ACT(i+1) + 2*average_walking_distance/c17_input_walking_speed);
             
            C_veh(i+1,k) = c16_input_toll + parking_pricing(i,1)*parkingduration_expectation/60 + C_drive(i+1,k) + C_time(i+1,k);
            
            theta = c28_input_theta;
            
            if sensitivity_parameter == 4
                input_park_ride_price = c21_input_park_ride_price + c21_input_park_ride_price*sensitivity_percentage;
            else
                input_park_ride_price = c21_input_park_ride_price;
            end
            
            if sensitivity_parameter == 5
                input_PT_price = c27_input_PT_price + c27_input_PT_price*sensitivity_percentage;
            else
                input_PT_price = c27_input_PT_price;
            end
            
            C_pr(i+1,k) = input_PT_price + input_park_ride_price + ...
                input_value_of_time(k)*(c22_input_average_headway + 2*average_driving_distance_PT/PT_speed(i+1,1) + 2*average_PT_walking_distance/c17_input_walking_speed);
             
            if switch_on_real_time_information_parking_pr_capacity == 0
                 C_veh_APR(i+1,k) = C_veh(i+1,k) * PR / (A_3 + PR);
                 C_pr_APR(i+1,k) = C_pr(i+1,k) * A_3 / (A_3 + PR);                   
            elseif switch_on_real_time_information_parking_pr_capacity == 1             
                 C_veh_APR(i+1,k) = C_veh(i+1,k) * pr_capacity(i+1,1) / (availableparking_3(i+1,1) + pr_capacity(i+1,1));
                 C_pr_APR(i+1,k) = C_pr(i+1,k) * availableparking_3(i+1,1) / (availableparking_3(i+1,1) + pr_capacity(i+1,1));
            end
             
            share_veh_pr_decision(i+1,k) = exp( (C_pr_APR(i+1,k) - C_veh_APR(i+1,k)) / min(C_pr_APR(i+1,k),C_veh_APR(i+1,k)) ) /...
                 (exp( (C_pr_APR(i+1,k) - C_veh_APR(i+1,k)) / min(C_pr_APR(i+1,k),C_veh_APR(i+1,k)) ) + 1);  
         end    
         
         if enter_pr(i+1,1)*(1 - mean(share_veh_pr_decision(i+1,:),2)) <= pr_capacity(i+1,1)
             enterthearea_3(i+1,1) = enterthearea_3(i+1,1)*mean(share_veh_pr_decision(i+1,:),2);
             enter_pr(i+1,1) = enter_pr(i+1,1)*(1 - mean(share_veh_pr_decision(i+1,:),2));
         else   
             enterthearea_3(i+1,1) = enterthearea_3(i+1,1) - pr_capacity(i+1,1);
             enter_pr(i+1,1) = pr_capacity(i+1,1);
         end
         
         % Congestion pricing implementation, End:
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Start Searching:
         
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         starttosearch_3(i+1,1)= dot(enterthearea_3(1:i,1), prob);
%        clear prob;
        
%        [~,n_demand_VOT]= c1_input_n_demand(i+1);   %| n_demand
%        % demand matrix for VOT origins:
%        if enter_pr(i+1,1)*(1 - mean(share_veh_pr_decision(i+1,:),2)) <= pr_capacity(i+1,1)
%           matrix_demand_VOT(i+1,:) = n_demand_VOT'.* 0.77 * mean(share_veh_pr_decision(i+1,:),2);   %| /ns for VOT, (i x V)-matrix 
%        else   
%           matrix_demand_VOT(i+1,:) = n_demand_VOT'.* 0.77 - pr_capacity(i+1,1);
%        end
         
         for j= 1:4           
           matrix_ns_s_VOT(i+1,j) = dot(starttosearch_3(1:i,1)/4, prob);  
%          matrix_ns_s_VOT(i+1,j) = dot(matrix_demand_VOT(1:i,j), prob);
         end
         clear prob;
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Find, decide and access parking:
         
         findparking_3(i+1,1) = d2_transitions_n_s_p_new(availableparking_3(i+1,1),Ns_3(i+1,1),speed(i+1,1));               
         probability(i+1,1) = d2_transitions_n_s_p_new(availableparking_3(i+1,1),Ns_3(i+1,1),speed(i+1,1))/Ns_3(i+1,1);
         
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
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         % Determine total bias for all vehicles
         
%        bias_tot(i+1,1) = findparking_3(i+1,1) * bias(i+1,1);
         bias_tot(i+1,1) = decideforparking_3(i+1,1) * bias(i+1,1);
         E_driven_distance(i+1,1) = driven_distance_per_time_slice(i+1,1) + bias_tot(i+1,1);
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Depart Parking:
         
         prob =fliplr(diff(gamcdf((0:i)*60*t,parameter(1),parameter(2))));        
         departparking_1(i+1,1) = startinginarea*gampdf(i*60*t,parameter(1),parameter(2));
         departparking_1(i+1,1) = departparking_1(i+1,1)+  dot(decideforparking_3(1:i,1), prob(1:end));        
         clear prob;
         
%        Parking durations for P+R:
%        [~, gp_a2, gp_b2] = c20_input_park_ride_parkingduration(1,1);
         [~,~,gp_a2,gp_b2] = c3_input_parkingduration(1,1);
         delay_by_PT(i,1) = (c22_input_average_headway + 2*average_driving_distance_PT/PT_speed(i+1,1))/60;  
         park_ride_parkingduration_expectation = gp_a2 * gp_b2;
         real_park_ride_parkingduration = park_ride_parkingduration_expectation + delay_by_PT(i,1);
         gp_b2 = real_park_ride_parkingduration/gp_a2; %Update scale parameter
         prob =fliplr(diff(gamcdf((0:i)*60*t,gp_a2,gp_b2)));
         depart_pr(i+1,1) = Npr(1,1)*gampdf(i*60*t,gp_a2,gp_b2);
         depart_pr(i+1,1) = depart_pr(i+1,1)+  dot(enter_pr(1:i,1), prob(1:end));
         clear prob;    
         
%        Leave after parking
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         leavethearea_1(i+1,1) = dot(departparking_1(1:i,1),prob);
         clear prob;
%        Leave after parking + Throughtraffic
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);% UNIFORM DISTRIBUTION: 
         leavethearea_1(i+1,1) = leavethearea_1(i+1,1) + dot(enterthearea_1(1:i,1),prob); % Additionally throughtraffic
         clear prob;
         
%total
    starttosearch(i+1,1)  = starttosearch_3(i+1,1);
    findparking(i+1,1)    = findparking_3(i+1,1); 
    decideforparking(i+1,1) = decideforparking_3(i+1,1);
    departparking(i+1,1)  = departparking_1(i+1,1);
    leavethearea(i+1,1)   = leavethearea_1(i+1,1);
%hard part of modelling each individual transition event.
%%
    clc
    disp(['== Coding line: ' num2str(i+1) ' ==']);
    i = i + 1;
    if i == 8
       i = 8; 
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
% Determine total bias for all vehicles over all time slices
sum_bias_tot = sum(bias_tot(1:1440,1));
sum_E_driven_distance = sum(E_driven_distance(1:1440,1));
total_error = sum_bias_tot/sum_E_driven_distance;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
occupancy_3(:,1)=1-availableparking_3(:,1)/A_3;
occupancy_pr(:,1)=1-pr_capacity(:,1)/PR;
%cumulative
for j=1:i
    cum_enterthearea_1(j) = sum(enterthearea_1(1:j));
    cum_leavethearea_1(j) = sum(leavethearea_1(1:j));
        
    cum_enterthearea_3(j) = startinginarea+sum(enterthearea_3(1:j));
    cum_enter_pr(j) = Npr(1,1)+sum(enter_pr(1:j));
    cum_starttosearch_3(j)= startinginarea+sum(starttosearch_3(1:j));
    cum_findparking_3(j)  = startinginarea+sum(findparking_3(1:j));
    cum_decideforparking_3(j) = startinginarea+sum(decideforparking_3(1:j));

    cum_enterthearea(j)   = startinginarea+sum(enterthearea(1:j));
    cum_starttosearch(j)  = startinginarea+sum(starttosearch(1:j));
    cum_findparking(j)    = startinginarea+sum(findparking(1:j));
    cum_decideforparking(j) = startinginarea+sum(decideforparking(1:j));
    cum_departparking(j)  = sum(departparking(1:j));
    cum_leavethearea(j)   = sum(leavethearea(1:j));
end
total_searchers_onstreet=sum(enterthearea_3(1:T,1))
total_vehicles_parked=sum(enterthearea_3(1:T,1)) + sum(enter_pr(1:T,1));
% total_searchtime_onstreet=sum(cum_starttosearch_3-cum_findparking_3)*t % unit in hours.
total_searchtime_onstreet=sum(cum_starttosearch_3-cum_decideforparking_3)*t % unit in hours.
% total_searchtime_onstreet_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_findparking_3(1+11*60:16*60))*t % unit in hours.
total_searchtime_onstreet_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_decideforparking_3(1+11*60:16*60))*t % unit in hours.
total_searchtime_onstreet_deduct1min=total_searchtime_onstreet-total_searchers_onstreet*t % unit in hours.
% total_searchtime_onstreet_deduct1min_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_findparking_3(1+11*60:16*60))*t-total_searchers_onstreet*t % unit in hours.
total_searchtime_onstreet_deduct1min_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_decideforparking_3(1+11*60:16*60))*t-total_searchers_onstreet*t % unit in hours.

% average_searchtime_onstreet=total_searchtime_onstreet/sum(findparking_3)*60 % average searching time per person unit in minutes
average_searchtime_onstreet=total_searchtime_onstreet/sum(decideforparking_3)*60 % average searching time per person unit in minutes
%%
time=1:1440;

%% 1ST PLOT ----------------------------------------------------------------input (parking demand enterthearea)
figure
plot(time/60,cum_enterthearea_3)
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
background=1:a:1+1439*a;
plot(time(a1:a2)/60,cum_enterthearea_3(a1:a2)-background(a1:a2));hold on;
plot(time(a1:a2)/60,cum_starttosearch_3(a1:a2)-background(a1:a2),'--');hold on;
plot(time(a1:a2)/60,cum_decideforparking_3(a1:a2)-background(a1:a2),'k:');
hold off;   
% axis([10 18 -380 -280])
ax = gca;
% ax.XTick = 10:1:18;
% ax.YTick = -380:50:-280;
xlabel('Hour of the day')
ylabel('Transformed cumulative number of vehicles (only for parkers)')
legend('enter area','start search','decide for parking')

%% --------------------------------------------------------------------------occupancy
figure
plot(time/60,occupancy_3,'LineWidth',4)
hold on
plot(time/60,occupancy_pr,'g--','LineWidth',4)
xlabel('Hour of the day')
axis([0 24 0 1])
ylabel('Parking occupancy')
legend('Parking in area','P+R facility')
set(gca,'FontSize',24)
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

%% --------------------------------------------------------------------------occupancy with parking usage
% figure
% plot(time/60,occupancy_3_new,'k-.','LineWidth',3)
% hold on
% plot(time/60,occupancy_3_parking_usage,'k:','LineWidth',3)
% hold on
% plot(time/60,occupancy_pr,'r-','LineWidth',3)
% hold on
% plot(time/60,occupancy_pr_parking_usage,'r--','LineWidth',3)
% xlabel('Hour of the day')
% axis([0 24 0 1])
% ylabel('Parking occupancy')
% legend('Parking in area (without parking usage information)','Parking in area (with parking usage information)','P+R facility (without parking usage information)','P+R facility (with parking usage information)')
% set(gca,'FontSize',24)
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:0.1:1;
% yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
% percentsy = repmat('%', length(yticks),1);  %  equal to the size
% yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
% set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

%% -------------------------------------------------------------------------number of searchers
figure
plot(time/60,Ns,'r-','LineWidth',4)
hold on
plot(time/60,availableparking_3,'b--','LineWidth',4)
xlabel('Hour of the day')
ylabel('Number')
leg = legend('Searching cars in the area','Available parking spaces in the area');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
axis([0 24 0 30])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:5:30;

%% PLOT --------------------------------------------------------------------share of traffic
figure 
plot(time/60,Ns_3./(Ns+Nns));
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

%% ------------------------------------------------vehicles entering the area vs. P+R

figure
plot(time/60,enterthearea_1,'g:','LineWidth',4)
hold on
plot(time/60,enterthearea_3,'r-','LineWidth',4)
hold on
plot(time/60,enter_pr,'b--','LineWidth',4)
hold off
xlabel('Hour of the day')
ylabel('Number of vehicles in traffic flow')
legend('Enter the area by vehicle (external destination)', 'Enter the area by vehicle (internal destination)', 'Enter the area by P+R')
set(gca,'FontSize',24)
axis([0 24 0 inf])
ax = gca;
ax.XTick = 0:1:24;
% ax.YTick = 0:0.1:1;

figure
plot(time/60,cum_enterthearea_1,'g:','LineWidth',4)
hold on
plot(time/60,cum_enterthearea_3,'r-','LineWidth',4)
hold on
plot(time/60,cum_enter_pr,'b--','LineWidth',4)
hold off
xlabel('Hour of the day')
ylabel('Cumulated number of vehicles in traffic flow')
legend('Enter the area by vehicle (external destination)', 'Enter the area by vehicle (internal destination)', 'Enter the area by P+R')
set(gca,'FontSize',24)
axis([0 24 0 inf])
ax = gca;
ax.XTick = 0:1:24;
% ax.YTick = 0:0.1:1;

%% -------------------------------------------------------------------------enter the network, parking demand
figure
plot(time/60,cum_enterthearea,'b-','LineWidth',4)
xlabel('Hour of the day')
ylabel('Cumulated number of vehicles entering the area')
set(gca,'FontSize',24)
axis([0 24 0 inf])
ax = gca;
ax.XTick = 0:1:24;

%% -------------------------------------------------------------------------distribution of parking durations

[~, ~, a2, b2] = c3_input_parkingduration(1,1);
% a2 = 1.6;
% b2 = 142;
figure
load('agentdata')
duration=(actendtime-actstarttime)/60;%unit is minute
histogram(duration);
hold on

p_all = zeros(1200,1);
for i = 1:size(p_all,1)
    p_all(i,1) = gamcdf(i,a2,b2) - gamcdf(i-1,a2,b2); 
end
hold on
plot(1:1:1200,130000*p_all,'r-','LineWidth',3)
xlabel('Parking duration (minutes)')
ylabel('Frequency')
legend('gammaparameters =   1.6191   142.1816')
hold off
set(gca,'FontSize',24)
axis([0 1200 0 500])

%% -------------------------------------------------------------------------average searching time
for i=1:1:size(cum_findparking_3,2) %time 
    % the ith row corresponds to the round(cum_starttosearch_3(i),3)
for x=i:1:size(cum_findparking_3,2) %time
       if round(cum_findparking_3(x))==round(cum_starttosearch_3(i))
            searchingtime(i)=x-i;
            break
        end
end
end
avg_searching_time_plot = (1./probability-1)*total_searchers_onstreet/total_vehicles_parked;
figure
plot(time/60,avg_searching_time_plot,'LineWidth',4)
xlabel('Hour of the day')
ylabel('Average searching time before parking (in min)')
set(gca,'FontSize',24)
axis([0 24 0 20])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:1:20;

% % -------------------------------------------------------------------------average searching time (scenario (a) and (c))
% for i=1:1:size(cum_findparking_3,2) %time 
%     % the ith row corresponds to the round(cum_starttosearch_3(i),3)
% for x=i:1:size(cum_findparking_3,2) %time
%        if round(cum_findparking_3(x))==round(cum_starttosearch_3(i))
%             searchingtime(i)=x-i;
%             break
%         end
% end
% end
% avg_searching_time_plot = (1./probability-1)*total_searchers_onstreet/total_vehicles_parked;
% figure
% plot(time/60,avg_searching_time_plot,'r--','LineWidth',4)
% hold on
% plot(time/60,avg_searching_time_plot_scenario_c*(1-0.05),'b-','LineWidth',4)
% hold off
% xlabel('Hour of the day')
% ylabel('Average searching time before parking (in min)')
% leg = legend('Status quo (reference scenario (a))', 'Parking pricing - policy 1 (scenario (c))');
% leg.ItemTokenSize = [50,18];
% set(gca,'FontSize',24)
% axis([0 24 0 14])
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:1:14;

%% -------------------------------------------------------------------------choice of drivers for entering the area by P+R

share_veh_pr_decision_mean = mean(share_veh_pr_decision(:,:),2);
share_veh_pr_decision_mean(1,1) = share_veh_pr_decision_mean(2,1);

share_veh_pr_decision(1,:) = share_veh_pr_decision(2,:);
share_veh_pr_decision_movmean = movmean(share_veh_pr_decision,11);

figure
plot(time/60,1-share_veh_pr_decision_movmean(:,1),'r--','LineWidth',4)
hold on
plot(time/60,1-share_veh_pr_decision_movmean(:,2),'b-','LineWidth',4)
hold on
plot(time/60,1-share_veh_pr_decision_movmean(:,3),'g:','LineWidth',4)
hold on
plot(time/60,1-share_veh_pr_decision_movmean(:,4),'k-.','LineWidth',4)
hold off
xlabel('Hour of the day')
ylabel('Drivers switching to P+R')
leg = legend('VOT = 29.9 CHF/h', 'VOT = 25.4 CHF/h', 'VOT = 25.8 CHF/h', 'VOT = 17.2 CHF/h');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
axis([0 24 0 1])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(100.*yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

share_veh_pr_decision_mean_all_times = mean(share_veh_pr_decision_mean(:,:),1);

%% -------------------------------------------------------------------------choice of drivers for entering the region by vehicle

figure
plot(time/60,share_veh_pr_decision(:,1),'r--','LineWidth',4)
hold on
plot(time/60,share_veh_pr_decision(:,2),'b-','LineWidth',4)
hold on
plot(time/60,share_veh_pr_decision(:,3),'g:','LineWidth',4)
hold on
plot(time/60,share_veh_pr_decision(:,4),'k-.','LineWidth',4)
hold off
xlabel('Hour of the day')
ylabel('Choice of drivers for entering the region by vehicle')
leg = legend('VOT = 29.9 CHF/h', 'VOT = 25.4 CHF/h', 'VOT = 25.8 CHF/h', 'VOT = 17.2 CHF/h');
leg.ItemTokenSize = [50,18];
set(gca,'FontSize',24)
axis([0 24 0 1])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:0.1:1;

% ------------------------------------------------choice of drivers for entering the region by vehicle (for all scenarios)
% 
% share_veh_pr_decision_mean = mean(share_veh_pr_decision(:,:),2);
% share_veh_pr_decision_mean(1,1) = share_veh_pr_decision_mean(2,1);
% figure
% plot(time/60,share_veh_pr_decision_mean_scenario_b,'r--','LineWidth',4)
% hold on
% plot(time/60,share_veh_pr_decision_mean_scenario_c,'b-','LineWidth',4)
% hold on
% plot(time/60,share_veh_pr_decision_mean_scenario_d,'g:','LineWidth',4)
% hold on
% plot(time/60,share_veh_pr_decision_mean_scenario_e,'k-.','LineWidth',4)
% hold off
% xlabel('Hour of the day')
% ylabel('Choice of drivers for entering the region by vehicle')
% legend('Free P+R and free PT (scenario (b))', 'Parking pricing - policy 1 (scenario (c))', 'Congestion pricing - policy 2 (scenario (d))', 'Parking and congestion pricing (scenario (e))')
% set(gca,'FontSize',24)
% axis([0 24 0 1])
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:0.1:1;
% 
% share_veh_pr_decision_mean_all_times = mean(share_veh_pr_decision_mean(:,:),1);

%------------------------------------------------choice of drivers for entering the region by vehicle (extrapolation)

% share_veh_pr_decision_mean = mean(share_veh_pr_decision(:,:),2);
% share_veh_pr_decision_mean(1,1) = share_veh_pr_decision_mean(2,1);
% 
% constant = lsqcurvefit(@f, [0,0,0,0], time/60, share_veh_pr_decision_mean');
% c1 = constant(1)
% c2 = constant(2)
% c3 = constant(3)
% c4 = constant(4)
% xfit = 0:0.25:24;
% yfit = f(constant,xfit);
% 
% figure
% plot(time/60,share_veh_pr_decision_mean,'r:','LineWidth',4)
% hold on
% 
% plot(xfit,yfit,'b-','LineWidth',4)
% xlabel('Hour of the day')
% ylabel('Choice of drivers for entering the region by vehicle')
% set(gca,'FontSize',24)
% axis([0 24 0 1])
% ax = gca;
% ax.XTick = 0:1:24;
% ax.YTick = 0:0.1:1;
% 
% share_veh_pr_decision_mean_all_times = mean(share_veh_pr_decision_mean(:,:),1);

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

[average_searching_time, average_searching_distance, final_cum_revenue, cum_revenue_PR, cum_revenue_congestion_toll, cum_revenue_parking_pricing] = ...
    b_new_outputs_from_the_matrix(Nns, Ns, Np, Nns_1, Nns_3, Npr, speed, cum_enterthearea, cum_enterthearea_3, cum_enter_pr, decideforparking_3, parking_pricing, ...
    enterthearea_1, enterthearea_3, enter_pr, parkingduration_expectation);

% % average searching time during peaks (between 10th hour and 16th hour)
% average_searching_time = sum(avg_searching_time_plot(600:960,1))/361;

% save ('scenario_VOT_3DMFD.mat')
% save ('scenario_a_reference_3DMFD.mat')
% save ('scenario_b_free_PR_free_PT_3DMFD.mat')
% save ('scenario_c_parking_pricing_3DMFD.mat')
% save ('scenario_d_congestion_pricing_3DMFD.mat')
% save ('scenario_e_both_parking_and_congestion_pricing_3DMFD.mat')

end
