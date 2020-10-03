function c_outputs_plot_parking_pricing(parking_pricing)

figure
hold on

% plot(linspace(1,size(parking_pricing,1),size(parking_pricing,1)),parking_pricing,'b', 'LineWidth', 4, 'MarkerSize', 10)
% 
% % plot grouped over every 4 minutes:
% plot(4*linspace(1,size(parking_pricing(1:4:size(parking_pricing,1)),1),size(parking_pricing(1:4:size(parking_pricing,1)),1)),parking_pricing(1:4:size(parking_pricing,1)),'b', 'LineWidth', 4, 'MarkerSize', 10)

% % plot grouped over every 4 minutes and rounded to next 0.5 values (round(value*2)/2):
% plot(4*linspace(1,size(parking_pricing(1:4:size(parking_pricing,1)),1),size(parking_pricing(1:4:size(parking_pricing,1)),1)),...
%     round(2*parking_pricing(1:4:size(parking_pricing,1)))/2,'b', 'LineWidth', 4, 'MarkerSize', 10)
% axis([0,150,2,12])

% % plot grouped over every 5 minutes and rounded to next 0.5 values (round(value*2)/2):
% plot(1:5:size(parking_pricing,1),...
%     round(2*parking_pricing(1:5:size(parking_pricing,1)))/2, 'b', 'LineWidth', 4, 'MarkerSize', 10)

% % plot grouped over every 10 minutes and rounded to next 0.5 values (round(value*2)/2):
% plot(1:10:size(parking_pricing,1),...
%     round(2*parking_pricing(1:10:size(parking_pricing,1)))/2, 'b', 'LineWidth', 4, 'MarkerSize', 10)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_value = (linspace(1,size(parking_pricing,1),size(parking_pricing,1)))';
x_value = (1:5:size(parking_pricing,1))'./60;
% y_value = parking_pricing;
y_value = round(2*parking_pricing(1:5:size(parking_pricing,1)))/2;
% y_value = round(parking_pricing(1:5:size(parking_pricing,1)),1);

plot(x_value, y_value, 'b-', 'MarkerSize', 10, 'LineWidth', 4)

% coeffs = polyfit(x_value, y_value, 14);
% xfit = linspace(x_value(1), max(x_value), size(x_value,1));
% yfit = polyval(coeffs, xfit);
% hold on;
% plot(xfit, yfit, 'r-', 'LineWidth', 4);
% grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if h2_getGlobal_switch_on_damp_exp_coef == 0
%     plot(linspace(1,size(parking_pricing,1),size(parking_pricing,1)),parking_pricing,'b', 'LineWidth', 4, 'MarkerSize', 10)
% elseif h2_getGlobal_switch_on_damp_exp_coef == 1 
%     plot(linspace(1,size(parking_pricing,1),size(parking_pricing,1)),parking_pricing,'b--', 'LineWidth', 4, 'MarkerSize', 10)
% end    
 
% title('Parking pricing')
xlabel('Time (hours)')
ylabel('Parking price (in CHF)')
set(gca,'FontSize',24)
% axis([0,24,2,10])
% axis([0,24,2.4,4])

hold off

end