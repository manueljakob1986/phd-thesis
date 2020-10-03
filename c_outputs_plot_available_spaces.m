function c_outputs_plot_available_spaces(availableparking_3, total_garage_capacity)

% Total revenue from on-street parking:
% Idea: Number of vehicles n_s_p * on-street parking fee

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on

red = plot((1:1:size(availableparking_3,1))'./60, availableparking_3(1:1:size(availableparking_3,1)), ...
    'r-', 'MarkerSize', 10, 'LineWidth', 3, 'DisplayName', 'On-street');
blue = plot((1:1:size(total_garage_capacity,1))'./60, total_garage_capacity(1:1:size(total_garage_capacity,1)), ...
    'b-', 'MarkerSize', 10, 'LineWidth', 4, 'DisplayName', 'Garage');

xlabel('Time (hours)')
ylabel('Available parking spaces')
set(gca,'FontSize',24)
% axis([0,24,2,150])
legend([blue red],'Location','northeast')

hold off

end