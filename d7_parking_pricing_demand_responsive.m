function [parking_pricing, delta_searching_veh_available_parking_spots, final_delta_searching_veh_available_parking_spots] = ...
    d7_parking_pricing_demand_responsive(searching_veh, available_parking_spots, parking_price_oscillation, old_parking_pricing, old_delta_searching_veh_available_parking_spots, only_occupancy)
% Test:
% searching_veh=[6;6;7;8;8;9;9;9;9;9;8;8;8;7;7;7;6;6;5;5;4;4;4;3;3;3;2;2;2;2;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;];
% available_parking_spots=[10;10;9;9;9;8;8;8;8;8;7;7;7;6;6;6;5;5;4;4;3;3;3;2;2;2;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;];

[input_use_maximum_limit, input_parking_price_maximum_increase] = c10_input_maximum_parking_price_increase_per_time_unit;
parameter_y = c15_input_parameter_y;
a = size(searching_veh,1);
parking_pricing = zeros(a,1);
delta_parking_pricing = zeros(a,1);
delta_searching_veh_available_parking_spots = zeros(a,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change, Paper Review, February 2017 
for i = 1:size(old_parking_pricing,1)
    parking_pricing(i,1) = old_parking_pricing(i,1);
    delta_searching_veh_available_parking_spots(i,1) = old_delta_searching_veh_available_parking_spots(i,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial price
parking_pricing(1,1) = c8_input_parking_price;

if only_occupancy == 0
    delta_searching_veh_available_parking_spots(1,1) = searching_veh(1,1) / available_parking_spots(1,1);
elseif only_occupancy == 1
    delta_searching_veh_available_parking_spots(1,1) = 1 / available_parking_spots(1,1);
end
delta_parking_pricing(1,1) = parking_pricing(1,1) * (delta_searching_veh_available_parking_spots(1,1))^(1/parameter_y);

if a > 1
    for i = a
    % for i=2:a    
        if available_parking_spots(i,1) == 0
            %no free parking spaces available -> infinity fee
            parking_pricing(i,1) = inf;
        else
            if only_occupancy == 0
                delta_searching_veh_available_parking_spots(i,1) = searching_veh(i,1) / available_parking_spots(i,1) - searching_veh(i-1,1) / available_parking_spots(i-1,1);
            elseif only_occupancy == 1
                delta_searching_veh_available_parking_spots(i,1) = 1 / available_parking_spots(i,1) - 1 / available_parking_spots(i-1,1);
            end

            if delta_searching_veh_available_parking_spots(i,1) >= 0  
%%%%%%%%%%%% pricing related to delta_searching_veh_available_parking_spots - increasing %%%%%%%%%%%% 
                if h2_getGlobal_switch_on_damp_exp_coef == 1
                    delta_parking_pricing(i,1) = exp(-i * c11_input_parking_price_oscillation(parking_price_oscillation)) * c8_input_parking_price * (delta_searching_veh_available_parking_spots(i,1))^(1/parameter_y);  
                elseif h2_getGlobal_switch_on_damp_exp_coef == 0
                    delta_parking_pricing(i,1) = c8_input_parking_price * (delta_searching_veh_available_parking_spots(i,1))^(1/parameter_y);    
                end

                 if input_use_maximum_limit == 1
%               define delta_parking_pricing increase limit:
                    if delta_parking_pricing(i,1) > input_parking_price_maximum_increase     
                        parking_pricing(i,1) = parking_pricing(i-1,1) + input_parking_price_maximum_increase;
                    elseif delta_parking_pricing(i,1) <= input_parking_price_maximum_increase
                        parking_pricing(i,1) = parking_pricing(i-1,1) + delta_parking_pricing(i,1);
                    end
                 else
%               no increase limit:
                    parking_pricing(i,1) = parking_pricing(i-1,1) + delta_parking_pricing(i,1);
                 end    

            elseif delta_searching_veh_available_parking_spots(i,1) < 0  
%%%%%%%%%%%% pricing related to delta_searching_veh_available_parking_spots - decreasing %%%%%%%%%%%% 
                if h2_getGlobal_switch_on_damp_exp_coef == 1
                    delta_parking_pricing(i,1) = exp(-i * c11_input_parking_price_oscillation(parking_price_oscillation)) * c8_input_parking_price * (-delta_searching_veh_available_parking_spots(i,1))^(1/parameter_y);  
                elseif h2_getGlobal_switch_on_damp_exp_coef == 0
                    delta_parking_pricing(i,1) = c8_input_parking_price * (-delta_searching_veh_available_parking_spots(i,1))^(1/parameter_y);
                end

                 if input_use_maximum_limit == 1
%               define delta_parking_pricing decrease limit:
                    if delta_parking_pricing(i,1) > input_parking_price_maximum_increase
                        if parking_pricing(i-1,1) - input_parking_price_maximum_increase >= c8_input_parking_price
                            parking_pricing(i,1) = parking_pricing(i-1,1) - input_parking_price_maximum_increase;
                        else
                            parking_pricing(i,1) = c8_input_parking_price;
                        end   
                    elseif delta_parking_pricing(i,1) <= input_parking_price_maximum_increase
                        if parking_pricing(i-1,1) - delta_parking_pricing(i,1) >= c8_input_parking_price
                            parking_pricing(i,1) = parking_pricing(i-1,1) - delta_parking_pricing(i,1); 
                        else
                            parking_pricing(i,1) = c8_input_parking_price;
                        end   
                    end    
                 else
%               no increase limit:
                    if parking_pricing(i-1,1) - delta_parking_pricing(i,1) >= c8_input_parking_price
                        parking_pricing(i,1) = parking_pricing(i-1,1) - delta_parking_pricing(i,1); 
                    else
                        parking_pricing(i,1) = c8_input_parking_price;
                    end 
                 end    
            end    

    % %     Option: parking pricing with piece-wise function    
    %       if searching_veh(i,1) / available_parking_spots(i,1) <= 1
    %         parking_pricing(i,1) = c8_input_parking_price;  
    %       elseif searching_veh(i,1) / available_parking_spots(i,1) <= c10_input_parking_price_maximum
    % %         parking_pricing(i,1) = c8_input_parking_price * searching_veh(i,1) / available_parking_spots(i,1);  
    %         parking_pricing(i,1) = exp(- searching_veh(i,1) / available_parking_spots(i,1)) * c8_input_parking_price * searching_veh(i,1) / available_parking_spots(i,1);
    %       else
    %         parking_pricing(i,1) = c8_input_parking_price * c10_input_parking_price_maximum;      
    %       end

    %     Options without piece-wise function
    %     parking_pricing(i,1) = c8_input_parking_price * searching_veh(i,1) / available_parking_spots(i,1);
    %     parking_pricing(i,1) = parking_pricing(i-1,1) * searching_veh(i,1) / available_parking_spots(i,1);
        end
    end
end

% final delta searching vehicles / available parking spaces:
final_delta_searching_veh_available_parking_spots = delta_searching_veh_available_parking_spots(a,1);

end