function c_outputs_plots(matrix)

% matrix = matrix;    %input from a_matrix.m
a = size(matrix,1);

parking_pricing_switched_on = c13_input_switch_on_parking_pricing;

% figure
hold on

if or(parking_pricing_switched_on == 1 , or(parking_pricing_switched_on == 2, parking_pricing_switched_on == 3))
    plot(linspace(1,a,a),matrix(:,14),'b--',...       %| cumulative /ns
         linspace(1,a,a),matrix(:,15),'r--',...       %| cumulative ns/s
         linspace(1,a,a),matrix(:,16),'g--',...       %| cumulative s/p
         linspace(1,a,a),matrix(:,17),'m--',...       %| cumulative p/ns
         linspace(1,a,a),matrix(:,18),'k--')          %| cumulative ns/    

elseif parking_pricing_switched_on == 0
    plot(linspace(1,a,a),matrix(:,14),'b',...       %| cumulative /ns
         linspace(1,a,a),matrix(:,15),'r',...       %| cumulative ns/s
         linspace(1,a,a),matrix(:,16),'g',...       %| cumulative s/p
         linspace(1,a,a),matrix(:,17),'m',...       %| cumulative p/ns
         linspace(1,a,a),matrix(:,18),'k')          %| cumulative ns/
end 

hold off
 
title('Queuing diagram of vehicles (with Parking Pricing)')
xlabel('Time (min)')
ylabel('Cumulative number of vehicles')

% legend('Enter the area','Start to search','Access parking','Depart parking','Leave the area')
legend('Enter the area (with Parking Pricing)','Start to search (with Parking Pricing)','Access parking (with Parking Pricing)','Depart parking (with Parking Pricing)','Leave the area (with Parking Pricing)')

end