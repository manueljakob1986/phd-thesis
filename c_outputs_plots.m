function c_outputs_plots(matrix)

% matrix = matrix;    %input from a_matrix.m
a = size(matrix,1);

parking_pricing_switched_on = c13_input_switch_on_parking_pricing;

% figure
hold on

if or(parking_pricing_switched_on == 1 , parking_pricing_switched_on == 2)
    
%     plot(2*linspace(1,a/2,a/2),matrix(1:2:a,20),'r',...       %| cumulative ns/dgp
%          2*linspace(1,a/2,a/2),matrix(1:2:a,22),'k',...       %| cumulative s/dgp
%          2*linspace(1,a/2,a/2),matrix(1:2:a,33),'c',...       %| cumulative dgp/gp
%          2*linspace(1,a/2,a/2),matrix(1:2:a,34),'g',...       %| cumulative dgp/s
%          2*linspace(1,a/2,a/2),matrix(1:2:a,24),'b')          %| cumulative gp/ns

%     plot(2*linspace(1,a/2,a/2),matrix(1:2:a,18),'c',...       %| cumulative /ns
%          2*linspace(1,a/2,a/2),matrix(1:2:a,19),'r',...       %| cumulative ns/s
%          2*linspace(1,a/2,a/2),matrix(1:2:a,21),'k',...       %| cumulative s/p
%          2*linspace(1,a/2,a/2),matrix(1:2:a,23),'b',...       %| cumulative p/ns
%          2*linspace(1,a/2,a/2),matrix(1:2:a,25),'g')          %| cumulative ns/ 
     
%     plot(2*linspace(1,a/2,a/2),matrix(1:2:a,18),'b',...       %| cumulative /ns
%          2*linspace(1,a/2,a/2),matrix(1:2:a,19),'r',...       %| cumulative ns/s
%          2*linspace(1,a/2,a/2),matrix(1:2:a,20),'c',...       %| cumulative ns/dgp
%          2*linspace(1,a/2,a/2),matrix(1:2:a,21),'g',...       %| cumulative s/p
%          2*linspace(1,a/2,a/2),matrix(1:2:a,22),'m',...       %| cumulative s/dgp
%          2*linspace(1,a/2,a/2),matrix(1:2:a,33),'k',...       %| cumulative dgp/gp
%          2*linspace(1,a/2,a/2),matrix(1:2:a,34),'y',...       %| cumulative dgp/s
%          2*linspace(1,a/2,a/2),matrix(1:2:a,23),'b--',...     %| cumulative p/ns
%          2*linspace(1,a/2,a/2),matrix(1:2:a,24),'g--',...     %| cumulative gp/ns
%          2*linspace(1,a/2,a/2),matrix(1:2:a,25),'k--')        %| cumulative ns/

%     plot(linspace(1,a,a),matrix(:,20),'r',...       %| cumulative ns/dgp
%          linspace(1,a,a),matrix(:,22),'k',...       %| cumulative s/dgp
%          linspace(1,a,a),matrix(:,33),'c',...       %| cumulative dgp/gp
%          linspace(1,a,a),matrix(:,34),'g',...       %| cumulative dgp/s
%          linspace(1,a,a),matrix(:,24),'b')          %| cumulative gp/ns
%     
%     plot(linspace(1,a,a),matrix(:,18),'c',...       %| cumulative /ns
%          linspace(1,a,a),matrix(:,19),'r',...       %| cumulative ns/s
%          linspace(1,a,a),matrix(:,21),'k',...       %| cumulative s/p
%          linspace(1,a,a),matrix(:,23),'b',...       %| cumulative p/ns
%          linspace(1,a,a),matrix(:,25),'g')          %| cumulative ns/ 
    
    plot(linspace(1,a,a),matrix(:,18),'b',...       %| cumulative /ns
         linspace(1,a,a),matrix(:,19),'r',...       %| cumulative ns/s
         linspace(1,a,a),matrix(:,20),'c',...       %| cumulative ns/dgp
         linspace(1,a,a),matrix(:,21),'g',...       %| cumulative s/p
         linspace(1,a,a),matrix(:,22),'m',...       %| cumulative s/dgp
         linspace(1,a,a),matrix(:,33),'k',...       %| cumulative dgp/gp
         linspace(1,a,a),matrix(:,34),'y',...       %| cumulative dgp/s
         linspace(1,a,a),matrix(:,23),'b--',...     %| cumulative p/ns
         linspace(1,a,a),matrix(:,24),'g--',...     %| cumulative gp/ns
         linspace(1,a,a),matrix(:,25),'k--')        %| cumulative ns/    

elseif parking_pricing_switched_on == 0
    plot(linspace(1,a,a),matrix(:,18),'b--',...     %| cumulative /ns
         linspace(1,a,a),matrix(:,19),'r--',...     %| cumulative ns/s
         linspace(1,a,a),matrix(:,20),'c--',...     %| cumulative ns/dgp
         linspace(1,a,a),matrix(:,21),'g--',...     %| cumulative s/p
         linspace(1,a,a),matrix(:,22),'m--',...     %| cumulative s/dgp
         linspace(1,a,a),matrix(:,33),'k--',...     %| cumulative dgp/gp
         linspace(1,a,a),matrix(:,34),'y--',...     %| cumulative dgp/s
         linspace(1,a,a),matrix(:,23),'b-.',...     %| cumulative p/ns
         linspace(1,a,a),matrix(:,24),'g-.',...     %| cumulative gp/ns
         linspace(1,a,a),matrix(:,25),'k-.')        %| cumulative ns/
end 

hold off
 
title('Queuing diagram of vehicles')
xlabel('Time (min)')
ylabel('Cumulative number of vehicles')
set(gca,'FontSize',16) 

% legend('Go to off-street parking',...
% 'Switch to off-street parking','Access off-street parking',...
% 'Not access off-street parking','Depart off-street parking')

% legend('Enter the area','Start to search for on-street parking',...
% 'Access on-street parking',...
% 'Depart on-street parking',...
% 'Leave the area')

legend('Enter the area','Start to search for on-street parking','Go to off-street parking',...
'Access on-street parking','Switch to off-street parking','Access off-street parking',...
'Not access off-street parking','Depart on-street parking','Depart off-street parking',...
'Leave the area')

% legend('Enter the area (with Parking Pricing)','Start to search for on-street parking (with Parking Pricing)','Go to off-street parking (with Parking Pricing)',...
% 'Access on-street parking (with Parking Pricing)','Switch to off-street parking (with Parking Pricing)','Access off-street parking (with Parking Pricing)',...
% 'Not access off-street parking (with Parking Pricing)','Depart on-street parking (with Parking Pricing)','Depart off-street parking (with Parking Pricing)',...
% 'Leave the area (with Parking Pricing)')

end