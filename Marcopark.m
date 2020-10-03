clear all% SMALL CIRCLE
clc
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
enterarea = demand';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % demand decrease by 4%:
% enterarea = enterarea * 0.96;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % demand increase by 4%:
% enterarea = enterarea * 1.04;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    enterthearea_1(1:T,1)=enterarea*0.23; 
    enterthearea_2(1:T,1)=0; 
    enterthearea_3(1:T,1)=enterarea*0.77; 
%total:
    enterthearea(:,1)=enterarea; 
%%
    density=0;
    speed=v;  
    availableparking_2=A_2;
    availableparking_3=A_3-startinginarea;
    
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
h1_setGlobal_initial_parking_pricing(2.5);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%
i=1;
while i<1440
    %this part is about the number of vehicles in each state
    Nns_1(i+1,1) = Nns_1(i,1)+enterthearea_1(i,1)-leavethearea_1(i,1); % only throughtraffic
    Nns_2(i+1,1) = Nns_2(i,1)+enterthearea_2(i,1)-starttosearch_2(i,1)+departparking_2(i,1)-leavethearea_2(i,1); % only garage parked vehicles
    Nns_3(i+1,1) = Nns_3(i,1)+enterthearea_3(i,1)-starttosearch_3(i,1)+departparking_3(i,1)-leavethearea_3(i,1); % only on-street parked vehicles
    
    Np_2(i+1,1)  = Np_2(i,1)+findparking_2(i,1)-departparking_2(i,1);% only garage parked vehicles
%     Np_3(i+1,1)  = Np_3(i,1)+findparking_3(i,1)-departparking_3(i,1);% only on-street parked vehicles
    Np_3(i+1,1)  = Np_3(i,1)+decideforparking_3(i,1)-departparking_3(i,1);% only on-street parked vehicles
    
    Ns_2(i+1,1)  = Ns_2(i,1)+starttosearch_2(i,1)-findparking_2(i,1);
    
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
    
    Nns(i+1,1)   = Nns_1(i+1,1)+Nns_2(i+1,1)+Nns_3(i+1,1);
    Np(i+1,1)    = Np_2(i+1,1)+Np_3(i+1,1);
    Ns(i+1,1)    = Ns_2(i+1,1)+Ns_3(i+1,1);
    
    density(i+1,1)=(Nns(i+1,1)+Ns(i+1,1))/Lnetwork;
    if density(i+1,1)<=kc
        speed(i+1,1)=v;
    elseif density(i+1,1)<=kj
        speed(i+1,1)=Qmax/(kc-kj)*(1-kj/density(i+1,1));
    else
        speed(i+1,1)=0;
    end
    availableparking_2(i+1,1)=A_2-Np_2(i+1,1);
    availableparking_3(i+1,1)=A_3-Np_3(i+1,1);
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
        
        for j= 1:4
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

        [~, parking_probability(i+1,:), guessed_price_vector(i,1), tau(i,1), VOT_k1(i,1), VOT_k2(i,1), VOT_k3(i,1), VOT_k4(i,1), ...
            E_p_vot(i,1), travel_distance_cost(i,1), penalty_distance(i,:), ...
             total_costs(i,:), parking_pricing(i,1), delta_searching_veh_available_parking_spots(i,1),...
             delta_searching_veh_available_parking_spots_i_plus_one(i,1)] = ...
         d2_transitions_n_s_p(availableparking_3(i+1,1),Ns_3(i+1,1),L,speed(i+1,1),t,matrix_ns_s_VOT,availableparking_3(1:i,1),...
             Ns_3(1:i,1),speed(1:i,1),0,A_3, parking_pricing, delta_searching_veh_available_parking_spots,...
             Ns_3_k1(i+1,1),Ns_3_k2(i+1,1),Ns_3_k3(i+1,1),Ns_3_k4(i+1,1),decideforparking_3(1:i,1),starttosearch_3(1:i,1));              %| s/p
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
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
         
%          bias_tot(i+1,1) = findparking_3(i+1,1) * bias(i+1,1);
         bias_tot(i+1,1) = decideforparking_3(i+1,1) * bias(i+1,1);
         E_driven_distance(i+1,1) = driven_distance_per_time_slice(i+1,1) + bias_tot(i+1,1);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
         prob =fliplr(diff(gamcdf((0:i)*60*t,parameter(1),parameter(2))));
         departparking_3(i+1,1)= startinginarea*gampdf(i*60*t,parameter(1),parameter(2));
%          departparking_3(i+1,1)= departparking_3(i+1,1)+  dot(findparking_3(1:i,1), prob(1:end));
         departparking_3(i+1,1)= departparking_3(i+1,1)+  dot(decideforparking_3(1:i,1), prob(1:end));
         clear prob;    
                  
         prob = unifcdf(dist_now,1*R,7*R)-unifcdf(dist_last,1*R,7*R);
         leavethearea_3(i+1,1) = dot(departparking_3(1:i,1),prob);clear prob;
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
    if i == 845
       i = 845; 
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
occupancy_2(:,1)=1-availableparking_2(:,1)/A_2;
occupancy_3(:,1)=1-availableparking_3(:,1)/A_3;
%cumulative
for j=1:i
    cum_enterthearea_1(j) = sum(enterthearea_1(1:j));
    cum_leavethearea_1(j) = sum(leavethearea_1(1:j));
    
    cum_enterthearea_2(j) = sum(enterthearea_2(1:j));
    cum_starttosearch_2(j)= sum(starttosearch_2(1:j));
    cum_findparking_2(j)  = sum(findparking_2(1:j));
    cum_departparking_2(j)= sum(departparking_2(1:j));
    cum_leavethearea_2(j) = sum(leavethearea_2(1:j));
    
    cum_enterthearea_3(j) = startinginarea+sum(enterthearea_3(1:j));
    cum_starttosearch_3(j)= startinginarea+sum(starttosearch_3(1:j));
    cum_findparking_3(j)  = startinginarea+sum(findparking_3(1:j));
    cum_decideforparking_3(j) = startinginarea+sum(decideforparking_3(1:j));
    cum_departparking_3(j)= sum(departparking_3(1:j));
    cum_leavethearea_3(j) = sum(leavethearea_3(1:j));

    cum_enterthearea(j)   = startinginarea+sum(enterthearea(1:j));
    cum_starttosearch(j)  = startinginarea+sum(starttosearch(1:j));
    cum_findparking(j)    = startinginarea+sum(findparking(1:j));
    cum_decideforparking(j) = startinginarea+sum(decideforparking(1:j));
    cum_departparking(j)  = sum(departparking(1:j));
    cum_leavethearea(j)   = sum(leavethearea(1:j));
end
total_searchers_onstreet=sum(enterthearea_3(1:T,1))
% total_searchtime_onstreet=sum(cum_starttosearch_3-cum_findparking_3)*t % unit in hours.
total_searchtime_onstreet=sum(cum_starttosearch_3-cum_decideforparking_3)*t % unit in hours.
% total_searchtime_onstreet_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_findparking_3(1+11*60:16*60))*t % unit in hours.
total_searchtime_onstreet_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_decideforparking_3(1+11*60:16*60))*t % unit in hours.
total_searchtime_onstreet_deduct1min=total_searchtime_onstreet-total_searchers_onstreet*t % unit in hours.
% total_searchtime_onstreet_deduct1min_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_findparking_3(1+11*60:16*60))*t-total_searchers_onstreet*t % unit in hours.
total_searchtime_onstreet_deduct1min_worsthours=sum(cum_starttosearch_3(1+11*60:16*60))*t -sum(cum_decideforparking_3(1+11*60:16*60))*t-total_searchers_onstreet*t % unit in hours.

% average_searchtime_onstreet=total_searchtime_onstreet/sum(findparking_3)*60 % average searching time per person unit in minutes
average_searchtime_onstreet=total_searchtime_onstreet/sum(decideforparking_3)*60 % average searching time per person unit in minutes
total_searchtime_garage=sum(cum_starttosearch_2-cum_findparking_2)*t % unit in hours.
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
plot(time(a1:a2)/60,cum_starttosearch_3(a1:a2)-background(a1:a2),'--');
plot(time(a1:a2)/60,cum_findparking_3(a1:a2)-background(a1:a2),'k:');
% plot(time(a1:a2)/60,cum_departparking_3(a1:a2)-background(a1:a2));
% plot(time(a1:a2)/60,cum_leavethearea_3(a1:a2)-background(a1:a2));
hold off;   
axis([10 18 -380 -280])
ax = gca;
ax.XTick = 10:1:18;
ax.YTick = -380:50:-280;
xlabel('Hour of the day')
ylabel('Transformed cumulative number of vehicles (only for parkers)')
legend('enter area','start search','found parking')

%% --------------------------------------------------------------------------occupancy
figure
plot(time/60,occupancy_3,'LineWidth',2)
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
plot(time/60,availableparking_3)
xlabel('Hour of the day')
ylabel('Number')
legend('Searching vehicles','Number of available parking')
axis([0 24 0 40])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:5:40;

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
figure
plot(time/60,1./probability-1)
xlabel('Hour of the day')
ylabel('Average search time before parking (minutes)')
axis([0 24 0 15])
ax = gca;
ax.XTick = 0:1:24;
ax.YTick = 0:1:15;
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