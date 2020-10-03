function c_outputs_plot_delta_guessed_pricing(delta_searching_veh_available_parking_spots, delta_searching_veh_available_parking_spots_i_plus_one)

figure
hold on

final_delta_searching_veh_available_parking_spots_i_plus_one = zeros(size(delta_searching_veh_available_parking_spots_i_plus_one,1) - 1,1);
for i = 1:(size(delta_searching_veh_available_parking_spots_i_plus_one,1) - 1)
   if abs(delta_searching_veh_available_parking_spots_i_plus_one(i,1)) > 0.7
       final_delta_searching_veh_available_parking_spots_i_plus_one(i,1) = delta_searching_veh_available_parking_spots(i+1,1); 
   else
       final_delta_searching_veh_available_parking_spots_i_plus_one(i,1) = delta_searching_veh_available_parking_spots_i_plus_one(i,1);
   end
    
end

% x_value = (1:size(delta_searching_veh_available_parking_spots_i_plus_one,1))'./60;
% y_value = delta_searching_veh_available_parking_spots_i_plus_one(1:size(delta_searching_veh_available_parking_spots_i_plus_one,1));
x_value = (1:size(final_delta_searching_veh_available_parking_spots_i_plus_one,1))'./60;
y_value = final_delta_searching_veh_available_parking_spots_i_plus_one(1:size(final_delta_searching_veh_available_parking_spots_i_plus_one,1));

red = plot(x_value, y_value, 'r--', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Predicted');

x_value2 = (1:(size(delta_searching_veh_available_parking_spots,1) - 1))'./60;
y_value2 = delta_searching_veh_available_parking_spots(2:size(delta_searching_veh_available_parking_spots,1));

blue = plot(x_value2, y_value2, 'b-', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Real');

xlabel('Time (hours)')
ylabel('Ratio between searchers and available parking spaces')
set(gca,'FontSize',24)
% axis([0,24,-0.9,0.9])
axis([9.5,16.5,-0.9,0.9])
legend([red blue],'Location','northeast')

hold off

percentage_delta_searching_veh_available_parking_spots = zeros(size(final_delta_searching_veh_available_parking_spots_i_plus_one,1) - 1,1);
for i = 1:(size(final_delta_searching_veh_available_parking_spots_i_plus_one,1) - 1)
   percentage_delta_searching_veh_available_parking_spots(i,1) = ...
       (final_delta_searching_veh_available_parking_spots_i_plus_one(i,1) - delta_searching_veh_available_parking_spots(i+1,1));
%    percentage_delta_searching_veh_available_parking_spots(i,1) = ...
%        (delta_searching_veh_available_parking_spots_i_plus_one(i,1) - delta_searching_veh_available_parking_spots(i+1,1)) / delta_searching_veh_available_parking_spots(i+1,1);   
end
mean(percentage_delta_searching_veh_available_parking_spots,1)
    
end