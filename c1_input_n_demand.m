function [n_demand_total,n_demand]=c1_input_n_demand(i,demand_options)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATSim data import:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('load_traffic_demand')

[demand,~] = hist(act_start_time_seconds./60, linspace(1,1440,1440));
n_demand_total = demand(1,i);


% This is the demand per VOT origin:
value_of_time = c7_input_value_of_time;
V = size(value_of_time,1);
n_demand = zeros(V,1);

n_demand(1,1) = cumulative_demand_VOT1(i,1);
n_demand(2,1) = cumulative_demand_VOT2(i,1);
n_demand(3,1) = cumulative_demand_VOT3(i,1);
n_demand(4,1) = cumulative_demand_VOT4(i,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if demand_options == 0
    % initial demand
else
    n_demand_total = n_demand_total * demand_options;
    n_demand(1,1) = n_demand(1,1) * demand_options;
    n_demand(2,1) = n_demand(2,1) * demand_options;
    n_demand(3,1) = n_demand(3,1) * demand_options;
    n_demand(4,1) = n_demand(4,1) * demand_options;
end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Poisson distribution:
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% lambda = 20;     % average arrival rate per hour = 20
% N_total = 160;   % this is the total demand.
% 
% parking_pricing_switched_on = c13_input_switch_on_parking_pricing;
%  
% % This is the demand per VOT origin:
% % Please check that sum over entries is equal to N_total !!!!!!
% 
% if or(parking_pricing_switched_on == 1 , parking_pricing_switched_on == 2)
%      N = [    
%          30
%          50
%          80];
% elseif parking_pricing_switched_on == 0
%      N = 160;
% end 
%  
% n_demand_total = N_total * poisspdf(i,lambda);
% 
%  
% value_of_time = c7_input_value_of_time;
% V = size(value_of_time,1);
% n_demand = zeros(V,1);
% for j=1:V
%   n_demand(j,1) = N(j) * poisspdf(i,lambda); 
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Gamma distribution:
% % A gamma distribution is used to define the entry time of the vehicles to 
% % the area, where the average arrival time is 20 min after the observation 
% % period starts (i.e., the shape parameter is 4 and the scale parameter is 5).  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  a1=4;          % SHAPE OF THE FUNCTION (similar to gauss)
%  b1=5;          % SCALE OF THE FUNCTION (mean=20, st.dev.=10), E(X) = a1*b1, VAR(X) = a1 * b1^2
%  N_total=200;   % this is the total demand.
%  
%  parking_pricing_switched_on = c13_input_switch_on_parking_pricing;
%  
% % This is the demand per VOT origin:
% % Please check that sum over entries is equal to N_total !!!!!!
% 
% if or(parking_pricing_switched_on == 1 , parking_pricing_switched_on == 2)
%      N = [    
%          40
%          60
%          100];
% elseif parking_pricing_switched_on == 0
%      N = 200;
% end
% 
%  n_demand_total=N_total*gampdf(i,a1,b1);
% 
%  
%  value_of_time = c7_input_value_of_time;
%  V = size(value_of_time,1);
%  n_demand = zeros(V,1);
%  for j=1:V
%    n_demand(j,1) = N(j)*gampdf(i,a1,b1); 
%  end    
 
end

