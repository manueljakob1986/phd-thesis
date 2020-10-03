function c_outputs_plot_parking_pricing_comparison(parking_pricing, guessed_price_vector)

figure
hold on


x_value = (1:5:size(parking_pricing,1))'./60;
y_value = round(2*parking_pricing(1:5:size(parking_pricing,1)))/2;

blue = plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4, 'DisplayName', 'Parking pricing');

x_value2 = (1:5:size(guessed_price_vector,1))'./60;
y_value2 = round(2*guessed_price_vector(1:5:size(guessed_price_vector,1)))/2;

red = plot(x_value2, y_value2, 'r-', 'MarkerSize', 10, 'LineWidth', 3, 'DisplayName', 'Predicted parking pricing');

xlabel('Time (hours)')
ylabel('Parking price (in CHF)')
set(gca,'FontSize',24)
axis([0,24,2,7.5])
legend([blue red],'Location','northeast')

hold off

end