function [final_parking_fee, parking_probability, guessed_price, tau, VOT_k1, VOT_k2, VOT_k3, VOT_k4, E_p_vot, ...
            travel_distance_cost, penalty_distance, total_costs, round_delta_searching_veh_available_parking_spots_i_plus_one] = ...
                d8_parking_pricing_probability(parking_pricing, searching_veh, available_parking_spots, delta_searching_veh_available_parking_spots, ...
                matrix_ns_s_VOT, v_all_t, t, A, L, N, Ns_3_k1, Ns_3_k2, Ns_3_k3, Ns_3_k4, decideforparking_3, starttosearch_3, parking_price_oscillation, A_total, constant_pricing)

% a = last time step i
a = size(searching_veh,1);
[input_use_maximum_limit, input_parking_price_maximum_increase] = c10_input_maximum_parking_price_increase_per_time_unit;
parameter_y = c15_input_parameter_y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Delta guessed pricing with 1/parameter_y approach %%% 
% 
sum_delta_searching_veh_available_parking_spots = 0;

% Guessed pricing via delta (searching veh)/(available parking spots) and 1/parameter_y approach:
% Assumption: take last fixed number of time interval steps (c12_input_guessed_pricing_agg_time_steps) and do their increases /decreases (from i-1 to i)
if a < c12_input_guessed_pricing_agg_time_steps
    fixed_guessed_pricing_agg_time_steps = a;
else
    fixed_guessed_pricing_agg_time_steps = c12_input_guessed_pricing_agg_time_steps;
end   
for i=(round(a - fixed_guessed_pricing_agg_time_steps) + 1):a
    if i-1 ~= 0
        if round(delta_searching_veh_available_parking_spots(i-1,1),6) ~= 0
            sum_delta_searching_veh_available_parking_spots = sum_delta_searching_veh_available_parking_spots + delta_searching_veh_available_parking_spots(i,1)/delta_searching_veh_available_parking_spots(i-1,1);
        end  
    else
        sum_delta_searching_veh_available_parking_spots = 1;
    end
end
sum_delta_searching_veh_available_parking_spots = sum_delta_searching_veh_available_parking_spots / fixed_guessed_pricing_agg_time_steps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change, Paper Review, February 2017
delta_searching_veh_available_parking_spots_i_plus_one = (delta_searching_veh_available_parking_spots(a,1)) * abs(sum_delta_searching_veh_available_parking_spots);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Computation of delta guessed pricing  
% if round(delta_searching_veh_available_parking_spots_i_plus_one,2) > 0  
%     delta_guessed_parking_pricing(a,1) = c8_input_parking_price * (delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots);
%     if a-1 ~= 0       
%              if input_use_maximum_limit == 1
% %               define delta_guessed_parking_pricing increase limit:
%                 if delta_guessed_parking_pricing(a,1) > input_parking_price_maximum_increase
%                     guessed_price = parking_pricing(a,1) + input_parking_price_maximum_increase;
%                 elseif delta_guessed_parking_pricing(a,1) <= input_parking_price_maximum_increase
%                     guessed_price = parking_pricing(a,1) + delta_guessed_parking_pricing(a,1);
%                 end
%              else
% %               no increase limit:
%                 guessed_price = parking_pricing(a,1) + delta_guessed_parking_pricing(a,1); 
%              end    
%     else
%        guessed_price = c8_input_parking_price + delta_guessed_parking_pricing(a,1); 
%     end
%         
% elseif round(delta_searching_veh_available_parking_spots_i_plus_one,2) == 0
%     if delta_searching_veh_available_parking_spots(a,1) < 0
%         delta_guessed_parking_pricing(a,1) = c8_input_parking_price * (-delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots);
%     else
%         delta_guessed_parking_pricing(a,1) = c8_input_parking_price * (delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots);
%     end
%     
%     guessed_price = c8_input_parking_price + delta_guessed_parking_pricing(a,1);
%     
% elseif round(delta_searching_veh_available_parking_spots_i_plus_one,2) < 0 
%     delta_guessed_parking_pricing(a,1) = c8_input_parking_price * (-delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots);        
%     if a-1 ~= 0       
%              if input_use_maximum_limit == 1
% %               define delta_guessed_parking_pricing decrease limit:
%                 if delta_guessed_parking_pricing(a,1) > input_parking_price_maximum_increase   
%                     if parking_pricing(a,1) - input_parking_price_maximum_increase >= 0
%                         guessed_price = parking_pricing(a,1) - input_parking_price_maximum_increase;
%                     else
%                         guessed_price = 0;
%                     end 
%                 elseif delta_guessed_parking_pricing(a,1) <= input_parking_price_maximum_increase
%                     if parking_pricing(a,1) - delta_guessed_parking_pricing(a,1) >= 0
%                         guessed_price = parking_pricing(a,1) - delta_guessed_parking_pricing(a,1);
%                     else
%                         guessed_price = 0;
%                     end 
%                 end
%              else
% %               no increase limit:
%                 if parking_pricing(a,1) - delta_guessed_parking_pricing(a,1) >= 0
%                     guessed_price = parking_pricing(a,1) - delta_guessed_parking_pricing(a,1); 
%                 else
%                     guessed_price = 0;
%                 end 
%              end    
%     else
%         guessed_price = c8_input_parking_price - delta_guessed_parking_pricing(a,1); 
%     end       
% end  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expected mean value, Compute expected price for VOT origin k at time i:
V = size(matrix_ns_s_VOT,2);
E_p_vot = 0;
count_a_skip = 0;
input_value_of_time = c7_input_value_of_time;
for k=1:V
    for i=1:a
        if sum(matrix_ns_s_VOT(i,:)) ~= 0
            E_p_vot = E_p_vot + matrix_ns_s_VOT(i,k) / sum(matrix_ns_s_VOT(i,:)) * input_value_of_time(k);
        else
            count_a_skip = count_a_skip + 1;
        end
    end
end
if a ~= 0
    E_p_vot = 1/(a - (count_a_skip/V)) * E_p_vot;
    if isnan(E_p_vot)
        E_p_vot = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute expected price for VOT origin k at time i:
input_value_of_time = c7_input_value_of_time;

VOT_k1 = input_value_of_time(1);
VOT_k2 = input_value_of_time(2);
VOT_k3 = input_value_of_time(3);
VOT_k4 = input_value_of_time(4);

% VOT_k1 = Ns_3_k1/ N * input_value_of_time(1);
% VOT_k2 = Ns_3_k2/ N * input_value_of_time(2);
% VOT_k3 = Ns_3_k3/ N * input_value_of_time(3);
% VOT_k4 = Ns_3_k4/ N * input_value_of_time(4);

% VOT_k1 = 0;
% VOT_k2 = 0;
% VOT_k3 = 0;
% VOT_k4 = 0;
% count_a_skip = 0;
% for i=1:a
%     if sum(matrix_ns_s_VOT(i+1,:)) ~= 0
%         VOT_k1 = VOT_k1 + matrix_ns_s_VOT(i+1,1) / sum(matrix_ns_s_VOT(i+1,:)) * input_value_of_time(1);
%         VOT_k2 = VOT_k2 + matrix_ns_s_VOT(i+1,2) / sum(matrix_ns_s_VOT(i+1,:)) * input_value_of_time(2);
%         VOT_k3 = VOT_k3 + matrix_ns_s_VOT(i+1,3) / sum(matrix_ns_s_VOT(i+1,:)) * input_value_of_time(3);
%         VOT_k4 = VOT_k4 + matrix_ns_s_VOT(i+1,4) / sum(matrix_ns_s_VOT(i+1,:)) * input_value_of_time(4);
%     else
%         count_a_skip = count_a_skip + 1;
%     end
% end
% 
% if a ~= 0
%     VOT_k1 = 1/(a - count_a_skip) * VOT_k1;
%     VOT_k2 = 1/(a - count_a_skip) * VOT_k2;
%     VOT_k3 = 1/(a - count_a_skip) * VOT_k3;
%     VOT_k4 = 1/(a - count_a_skip) * VOT_k4;
%     if isnan(VOT_k1)
%         VOT_k1 = 0;
%     end
%     if isnan(VOT_k2)
%         VOT_k2 = 0;
%     end
%     if isnan(VOT_k3)
%         VOT_k3 = 0;
%     end
%     if isnan(VOT_k4)
%         VOT_k4 = 0;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Travel time from one parking place to next parking place
% A = parking supply at time i (=a)
% L = length of network
if A * v_all_t(a,1) == 0
    tau = inf;
else
    tau = L / (A * v_all_t(a,1)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Travel distance to next parking space

travel_distance_cost = 0;
if A ~= 0
    travel_distance_cost = c9_input_price_per_distance * L / A;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % With coding: based on static conditions and probability of finding
% % parking: !!!regarding cruising distance!!!
% ACD = 0;
% 
% % Determine average cruising distance
% if A <= A_total
%     if A >= round(searching_veh(a,1))
%        if A + 1 == 0
%            ACD = 0;
%        elseif round(searching_veh(a,1)) == 0
%            ACD = 0;
%        elseif A == 0
%            ACD = 0;  
%        else    
%            ACD = L - (L * ((1 - (1 / (A + 1))*(1 - (1 - 1/round(searching_veh(a,1)))^(A + 1))) + ...
%                     (A/round(searching_veh(a,1)))*exp(-1*A/round(searching_veh(a,1)))*...
%                     ((1/A)*(1 - round(searching_veh(a,1))) + 13/36 - 1/(4*(round(searching_veh(a,1)))^2) - ...
%                     1/(9*(round(searching_veh(a,1)))^3) + (A/round(searching_veh(a,1)))*((1/(24*(round(searching_veh(a,1)))^3)) - 1/24))));
%        end
%             
%     elseif A < round(searching_veh(a,1))
%         if searching_veh(a,1) >= 0
%             if searching_veh(a,1) <= A
%                if A + 1 == 0
%                    ACD = 0;
%                elseif round(searching_veh(a,1)) == 0
%                    ACD = 0;
%                elseif A == 0
%                    ACD = 0;  
%                else 
%                    ACD = L*(A/round(searching_veh(a,1)))^2 - (L * (((A/round(searching_veh(a,1))) - (1 / (A + 1))*(1 - (1 - 1/round(searching_veh(a,1)))^(A + 1))) + ...
%                             (A/round(searching_veh(a,1)))*exp(-1*A/round(searching_veh(a,1)))*...
%                             ((1/A) - 1 - 1/(4*(round(searching_veh(a,1)))^2) - 1/(9*(round(searching_veh(a,1)))^3) + ...
%                             (A/(24*(round(searching_veh(a,1)))^4)) + (A^2/(4*(round(searching_veh(a,1)))^2)) + (A^3/(9*(round(searching_veh(a,1)))^3)) - ...
%                             (A^4/(24*(round(searching_veh(a,1)))^4))))); 
%                end        
%             else               
% %               in case of infinite cruising distance, take old approach to determine average cruising distance:  
%                 total_searching_distance_until_now = 0;
%                 total_number_of_veh_searching = 0;
%                 
%                 for i=2:(a-1)
%                 %   total searching distance until now (in km), no consideration of
%                 %   iteration step 1 and last iteration step a
%                     total_searching_distance_until_now = total_searching_distance_until_now + searching_veh(i,1)*v_all_t(i,1)*t; % unit is km.
%                 %   total number of vehicles, that were searching (refers to transition  
%                 %   count of vehicles n_ns_s), no consideration of iteration step 1    
%                     total_number_of_veh_searching = total_number_of_veh_searching + sum(matrix_ns_s_VOT(i,:));
%                 end
% 
%                 % Compute average searching time (for vehicles still searching and already
%                 % having parked)
%                 % = average out of "average searching time for vehicles still searching" 
%                 % and "average searched time (now parked)"
%                 if total_number_of_veh_searching  ~= 0
%                     ACD = total_searching_distance_until_now / total_number_of_veh_searching;
%                 end 
%             end
%         end
%     end
% end
% 
% % Compute the penalty term for the already driven distance:    
% penalty_distance = c9_input_price_per_distance * ACD;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average cruising time

for i_minus_CT_max = 1:a
    if sum(starttosearch_3(1:i_minus_CT_max),1) > sum(decideforparking_3(1:a),1)
        break
    end
end
i_minus_CT_max = i_minus_CT_max - 1;
CT_max = a - i_minus_CT_max;
ACT = CT_max/2;

% Compute the penalty term for the already driven distance:
penalty_distance = zeros(1,4);
penalty_distance(1,1) = VOT_k1 * ACT;
penalty_distance(1,2) = VOT_k2 * ACT;
penalty_distance(1,3) = VOT_k3 * ACT;
penalty_distance(1,4) = VOT_k4 * ACT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization scheme
round_delta_searching_veh_available_parking_spots_i_plus_one = round(delta_searching_veh_available_parking_spots_i_plus_one,2);
final_parking_fee = c8_input_parking_price;
if or(constant_pricing == 0, constant_pricing == 2)
    
    if a > 1
        if constant_pricing == 2 % only occupancy-responsive pricing scheme
            if a > 900
                round_delta_searching_veh_available_parking_spots_i_plus_one = round(delta_searching_veh_available_parking_spots_i_plus_one,3);
            else
                round_delta_searching_veh_available_parking_spots_i_plus_one = round(delta_searching_veh_available_parking_spots_i_plus_one,2);
            end
        end
        
        if round_delta_searching_veh_available_parking_spots_i_plus_one > 0
            final_parking_fee = parking_pricing(a-1,1) + input_parking_price_maximum_increase;
            
        elseif round_delta_searching_veh_available_parking_spots_i_plus_one == 0
            final_parking_fee = parking_pricing(a-1,1);
            if constant_pricing == 0 % only responsive pricing scheme
                if a > 1080
%                 if a > 800 % Demand in-/decrease (+/- 4 %), Supply in-/decrease (+/- 4 %)
                    final_parking_fee = parking_pricing(a-1,1) - input_parking_price_maximum_increase;
                end
            end

        elseif round_delta_searching_veh_available_parking_spots_i_plus_one < 0
            if delta_searching_veh_available_parking_spots(a,1) < 0
                final_parking_fee = min(min(min(min((tau * VOT_k1 + travel_distance_cost + penalty_distance(1,1)) / ...
                    ((-delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots))...
                    ,(tau * VOT_k2 + travel_distance_cost + penalty_distance(1,2)) / ...
                    ((-delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)))...
                    ,(tau * VOT_k3 + travel_distance_cost + penalty_distance(1,3)) / ...
                    ((-delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)))...
                    ,(tau * VOT_k4 + travel_distance_cost + penalty_distance(1,4)) / ...
                    ((-delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)))...
                    ,(input_parking_price_maximum_increase) / ...
                    ((-delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)));
            elseif ((delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)) == 0
                final_parking_fee = parking_pricing(a-1,1) - input_parking_price_maximum_increase;
            else
                final_parking_fee = min(min(min(min((tau * VOT_k1 + travel_distance_cost + penalty_distance(1,1)) / ...
                    ((delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots))...
                    ,(tau * VOT_k2 + travel_distance_cost + penalty_distance(1,2)) / ...
                    ((delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)))...
                    ,(tau * VOT_k3 + travel_distance_cost + penalty_distance(1,3)) / ...
                    ((delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)))...
                    ,(tau * VOT_k4 + travel_distance_cost + penalty_distance(1,4)) / ...
                    ((delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)))...
                    ,(input_parking_price_maximum_increase) / ...
                    ((delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots)));
            end
            
        end
        
        if parking_pricing(a-1,1) - input_parking_price_maximum_increase > final_parking_fee
            final_parking_fee = parking_pricing(a-1,1) - input_parking_price_maximum_increase;
            
        elseif parking_pricing(a-1,1) + input_parking_price_maximum_increase < final_parking_fee
            final_parking_fee = parking_pricing(a-1,1) + input_parking_price_maximum_increase;
        end
    end
    if final_parking_fee < c8_input_parking_price
        final_parking_fee = c8_input_parking_price;
    end
elseif constant_pricing == 1
    guessed_price = final_parking_fee;
end

parking_pricing(a,1) = final_parking_fee;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of delta guessed pricing
if or(constant_pricing == 0, constant_pricing == 2)
    if round_delta_searching_veh_available_parking_spots_i_plus_one > 0
        delta_guessed_parking_pricing(a,1) = final_parking_fee * (delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots);
        if a-1 ~= 0
            if input_use_maximum_limit == 1
                % define delta_guessed_parking_pricing increase limit:
                if delta_guessed_parking_pricing(a,1) > input_parking_price_maximum_increase
                    guessed_price = parking_pricing(a,1) + input_parking_price_maximum_increase;
                elseif delta_guessed_parking_pricing(a,1) <= input_parking_price_maximum_increase
                    guessed_price = parking_pricing(a,1) + delta_guessed_parking_pricing(a,1);
                end
            else
                % no increase limit:
                guessed_price = parking_pricing(a,1) + delta_guessed_parking_pricing(a,1);
            end
        else
            guessed_price = final_parking_fee + delta_guessed_parking_pricing(a,1);
        end
        
    elseif round_delta_searching_veh_available_parking_spots_i_plus_one == 0
          guessed_price = final_parking_fee;

    elseif round_delta_searching_veh_available_parking_spots_i_plus_one < 0
        delta_guessed_parking_pricing(a,1) = final_parking_fee * (-delta_searching_veh_available_parking_spots(a,1))^(1/parameter_y) * abs(sum_delta_searching_veh_available_parking_spots);
        if a-1 ~= 0
            if input_use_maximum_limit == 1
                % define delta_guessed_parking_pricing decrease limit:
                if delta_guessed_parking_pricing(a,1) > input_parking_price_maximum_increase
                    if parking_pricing(a,1) - input_parking_price_maximum_increase >= 0
                        guessed_price = parking_pricing(a,1) - input_parking_price_maximum_increase;
                    else
                        guessed_price = 0;
                    end
                elseif delta_guessed_parking_pricing(a,1) <= input_parking_price_maximum_increase
                    if parking_pricing(a,1) - delta_guessed_parking_pricing(a,1) >= 0
                        guessed_price = parking_pricing(a,1) - delta_guessed_parking_pricing(a,1);
                    else
                        guessed_price = 0;
                    end
                end
            else
                % no increase limit:
                if parking_pricing(a,1) - delta_guessed_parking_pricing(a,1) >= 0
                    guessed_price = parking_pricing(a,1) - delta_guessed_parking_pricing(a,1);
                else
                    guessed_price = 0;
                end
            end
        else
            guessed_price = final_parking_fee - delta_guessed_parking_pricing(a,1);
        end
    end
end

if guessed_price < c8_input_parking_price
    guessed_price = c8_input_parking_price;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Average values based on n_ns_s needed for decision model:
% V = size(matrix_ns_s_VOT,2);
% user_group_percentage_k1 = 0;
% user_group_percentage_k2 = 0;
% user_group_percentage_k3 = 0;
% user_group_percentage_k4 = 0;
% count_a_skip = 0;
% 
% for i=1:a
%     if sum(matrix_ns_s_VOT(i+1,:)) ~= 0
%         user_group_percentage_k1 = user_group_percentage_k1 + matrix_ns_s_VOT(i+1,1) / sum(matrix_ns_s_VOT(i+1,:));
%         user_group_percentage_k2 = user_group_percentage_k2 + matrix_ns_s_VOT(i+1,2) / sum(matrix_ns_s_VOT(i+1,:));
%         user_group_percentage_k3 = user_group_percentage_k3 + matrix_ns_s_VOT(i+1,3) / sum(matrix_ns_s_VOT(i+1,:));
%         user_group_percentage_k4 = user_group_percentage_k4 + matrix_ns_s_VOT(i+1,4) / sum(matrix_ns_s_VOT(i+1,:));
%     else
%         count_a_skip = count_a_skip + 1;
%     end
% end
% 
% if a ~= 0
%     user_group_percentage_k1 = 1/(a - count_a_skip) * user_group_percentage_k1;
%     user_group_percentage_k2 = 1/(a - count_a_skip) * user_group_percentage_k2;
%     user_group_percentage_k3 = 1/(a - count_a_skip) * user_group_percentage_k3;
%     user_group_percentage_k4 = 1/(a - count_a_skip) * user_group_percentage_k4;
%     if isnan(user_group_percentage_k1)
%         user_group_percentage_k1 = 0;
%     end
%     if isnan(user_group_percentage_k2)
%         user_group_percentage_k2 = 0;
%     end
%     if isnan(user_group_percentage_k3)
%         user_group_percentage_k3 = 0;
%     end
%     if isnan(user_group_percentage_k4)
%         user_group_percentage_k4 = 0;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision model:

parking_probability = zeros(1,4);
total_costs = zeros(1,4);
total_costs(1,1) = guessed_price + tau * VOT_k1 + travel_distance_cost + penalty_distance(1,1);
if final_parking_fee <= guessed_price + tau * VOT_k1 + travel_distance_cost + penalty_distance(1,1)
    
    if guessed_price == inf || tau == inf
        parking_probability(1,1) = 0;
    else    
        if L ~= 0
            parking_probability(1,1) = 1;
        end
    end
else %condition not fulfilled: car keeps on driving and does not park
    parking_probability(1,1) = 0;
end

total_costs(1,2) = guessed_price + tau * VOT_k2 + travel_distance_cost + penalty_distance(1,2);
if final_parking_fee <= guessed_price + tau * VOT_k2 + travel_distance_cost + penalty_distance(1,2)
    
    if guessed_price == inf || tau == inf
        parking_probability(1,2) = 0;
    else    
        if L ~= 0
            parking_probability(1,2) = 1;
        end
    end
else %condition not fulfilled: car keeps on driving and does not park
    parking_probability(1,2) = 0;
end

total_costs(1,3) = guessed_price + tau * VOT_k3 + travel_distance_cost + penalty_distance(1,3);
if final_parking_fee <= guessed_price + tau * VOT_k3 + travel_distance_cost + penalty_distance(1,3)
    
    if guessed_price == inf || tau == inf
        parking_probability(1,3) = 0;
    else    
        if L ~= 0
            parking_probability(1,3) = 1;
        end
    end
else %condition not fulfilled: car keeps on driving and does not park
    parking_probability(1,3) = 0;
end

total_costs(1,4) = guessed_price + tau * VOT_k4 + travel_distance_cost + penalty_distance(1,4);
if final_parking_fee <= guessed_price + tau * VOT_k4 + travel_distance_cost + penalty_distance(1,4)
    
    if guessed_price == inf || tau == inf
        parking_probability(1,4) = 0;
    else    
        if L ~= 0
            parking_probability(1,4) = 1;
        end
    end
else %condition not fulfilled: car keeps on driving and does not park
    parking_probability(1,4) = 0;
end

% if parking_probability > 1
%     parking_probability = 1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute final probability value:
% Assumption: all parking spots have the same price

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change Start, Paper Review, February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% final_parking_fee = parking_pricing(a,1);
% if constant_pricing == 0
%     if guessed_price < c8_input_parking_price
%         guessed_price = c8_input_parking_price;
%     end
%     
% %     final_parking_fee = guessed_price + tau * E_p_vot + travel_distance_cost + penalty_distance;    
% %     final_parking_fee = guessed_price;
%     final_parking_fee = guessed_price + (tau * E_p_vot + travel_distance_cost + penalty_distance)/60;
%     total_costs = final_parking_fee;
% 
%     if a > 1
%         if parking_pricing(a-1,1) - input_parking_price_maximum_increase > final_parking_fee
%             final_parking_fee = parking_pricing(a-1,1) - input_parking_price_maximum_increase;
%             
% %           if final_parking_fee > guessed_price + tau * E_p_vot + travel_distance_cost + penalty_distance
% %               final_parking_fee = guessed_price + tau * E_p_vot + travel_distance_cost + penalty_distance;
% %           end
% 
%         elseif parking_pricing(a-1,1) + input_parking_price_maximum_increase < final_parking_fee
%             final_parking_fee = parking_pricing(a-1,1) + input_parking_price_maximum_increase;
%         end
%     end
%     if final_parking_fee < c8_input_parking_price
%         final_parking_fee = c8_input_parking_price;
%     end
% elseif constant_pricing == 1
%     guessed_price = final_parking_fee;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change End, Paper Review, February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if final_parking_fee <= guessed_price
% if final_parking_fee <= guessed_price + tau * E_p_vot + travel_distance_cost + penalty_distance
% 
%     if guessed_price == inf || tau == inf
%         parking_probability = 0;
%     else    
%         if L ~= 0
%             parking_probability = 1; % probability = A * 1 / A (all parking spots behave equally)
% %           parking_probability = A/L * (guessed_price + tau * E_p_vot + penalty_distance); 
% %           parking_probability = (1/L * (guessed_price + tau * E_p_vot + penalty_distance))^A; 
%         end
%     end
% else %condition not fulfilled: car keeps on driving and does not park
%     parking_probability = 0;
% end
% 
% if parking_probability > 1
%     parking_probability = 1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute probabilty, that there are no other parking spots between
% two parking spots in the used parking_probability computation
% if A == 0
%     parking_probability_no_other_parking_spots = 1;
% else    
%     parking_probability_no_other_parking_spots = (1 - 1/A)^(A-1);
% end
% parking_probability_no_other_parking_spots = 1;
% 
% if parking_probability_no_other_parking_spots > 1
%     parking_probability_no_other_parking_spots = 1;
% end

% parking_pricing(a,1)
% searching_veh
% available_parking_spots
% guessed_price
% tau
% E_p_vot
% t_average
% travel_distance_cost
% penalty_distance
% parking_probability
% parking_probability_no_other_parking_spots

end