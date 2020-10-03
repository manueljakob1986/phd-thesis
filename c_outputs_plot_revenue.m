function c_outputs_plot_revenue(matrix,parking_pricing)
matrix = matrix; %input from a_matrix.m

% Total revenue from on-street parking:
% Idea: Number of vehicles n_s_p * on-street parking fee

figure
revenue(:,1) = zeros(size(matrix(:,10),1) - 1,1);
cum_revenue(:,1) = zeros(size(matrix(:,10),1) - 1,1);
correct_parking_pricing(:,1) = zeros(size(matrix(:,10),1) - 1,1);

hold on
% Get the correct parking pricing value (update only every 5 minutes and
% rounded to next 0.5 CHF):
for j = 1:5:size(parking_pricing(:,1),1)
    
    correct_parking_pricing(j,1) = round(2*parking_pricing(j,1))/2;
    if j ~= size(parking_pricing(:,1),1)
        correct_parking_pricing(j + 1,1) = correct_parking_pricing(j,1);
        correct_parking_pricing(j + 2,1) = correct_parking_pricing(j,1);
        correct_parking_pricing(j + 3,1) = correct_parking_pricing(j,1);
        correct_parking_pricing(j + 4,1) = correct_parking_pricing(j,1);
    end
end
    
% Get revenue value from corrected parking pricing value for every 5 minutes:   
for i = 2:size(matrix(:,10),1)
    revenue(i-1,1) = matrix(i,10)*correct_parking_pricing(i-1,1);

    if i ~= 2
        cum_revenue(i-1,1) = cum_revenue(i-2,1) + revenue(i-1,1);
    elseif i == 2
        cum_revenue(i-1,1) = revenue(i-1,1);
    end
end

plot(linspace(1,size(cum_revenue(:,1),1),size(cum_revenue(:,1),1)),cum_revenue(:,1),'b', 'LineWidth', 4, 'MarkerSize', 10)

% title('Cumulative revenue resulting from parking pricing')
xlabel('Time (min)')
ylabel('Cumulative revenue (in CHF)')
set(gca,'FontSize',24)
axis([0,1440,0,13000])
% axis([0,1440,0,7000])

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % use of polyfit:
% % Total revenue from on-street parking:
% % Idea: Number of vehicles n_s_p * on-street parking fee
% 
% figure
% 
% revenue_fit(:,1) = zeros(size(matrix(:,10),1) - 1,1);
% cum_revenue_fit(:,1) = zeros(size(matrix(:,10),1) - 1,1);
% correct_parking_pricing_fit(:,1) = zeros(size(matrix(:,10),1) - 1,1);
% 
% hold on
% 
% x_value = (1:5:size(parking_pricing,1))';
% y_value = round(2*parking_pricing(1:5:size(parking_pricing,1)))/2;
% coeffs = polyfit(x_value, y_value, 14);
% xfit = linspace(x_value(1), max(x_value), size(x_value,1));
% yfit = polyval(coeffs, xfit);
% 
% % Get the correct parking pricing value (update only every 5 minutes and
% % rounded to next 0.5 CHF):
% for j = 1:5:size(parking_pricing(:,1),1)
%     
%     correct_parking_pricing_fit(j,1) = round(2*yfit(1,(j-1)/5 + 1))/2;
%     if j ~= size(parking_pricing(:,1),1)
%         correct_parking_pricing_fit(j + 1,1) = correct_parking_pricing_fit(j,1);
%         correct_parking_pricing_fit(j + 2,1) = correct_parking_pricing_fit(j,1);
%         correct_parking_pricing_fit(j + 3,1) = correct_parking_pricing_fit(j,1);
%         correct_parking_pricing_fit(j + 4,1) = correct_parking_pricing_fit(j,1);
%     end
% end
%     
% % Get revenue value from corrected parking pricing value for every 5 minutes:   
% for i = 2:size(matrix(:,10),1)
%     revenue_fit(i-1,1) = matrix(i,10)*correct_parking_pricing_fit(i-1,1);
% 
%     if i ~= 2
%         cum_revenue_fit(i-1,1) = cum_revenue_fit(i-2,1) + revenue_fit(i-1,1);
%     elseif i == 2
%         cum_revenue_fit(i-1,1) = revenue_fit(i-1,1);
%     end
% end
% 
% plot(linspace(1,size(cum_revenue_fit(:,1),1),size(cum_revenue_fit(:,1),1)),cum_revenue_fit(:,1),'b', 'LineWidth', 4, 'MarkerSize', 10)
% 
% % title('Cumulative revenue resulting from parking pricing')
% xlabel('Time (min)')
% ylabel('Cumulative revenue (in CHF)')
% set(gca,'FontSize',24)
% axis([0,1440,0,13000])
% % axis([0,1440,0,7000])
% 
% hold off

end