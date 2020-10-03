function [matrix, parking_pricing, guessed_price_vector, E_p_vot, tau, cost_travelling_next_space, penalty_distance, total_costs, bias] = a_matrix(parking_price_oscillation)
%%---------MAIN INFO OF THE MATRIX:
% This file provides the main body of the transition matrix,
% each row represents a time slice where the 1st represents the initial conditions;
% and each colomn represents an variables, there are 18 colomns.
% The colomns represent: 
% |  states  |       4      |   5   |  6  |    7    |     five  transitions    |         13        |  CUMULATIVE transitions  |      19       |          20               |    21        |        22          |
% | Nns;Ns;Np|parking supply|density|speed| n_demand| /ns; ns/s; s/p; p/ns; ns/|cumulative n_demand| /ns; ns/s; s/p; p/ns; ns/|nleave_parkers | cumulative_nleave_parkers | Nns_parkers  | Nns_throughtraffic | 

%%---------ASSUMPTIONS (INPUTS):
% We categorize the inputs to this script into two:
% 1st part is fixed values defined in this script. They include, 
% 1st, size/length of the network, L;
% 2nd, Parking capacity, A; 
% 3rd, time length of each time slice, t; 
% 4th, maximum speed (free flow speed), v; 
% 5th, maximum flow, Qmax; 
% 6th, critical/optimal traffic density, kc; 
% 7th, jam density, kj. 

        L = 15.4; % unit:km
        A = 539; % unit: parking spaces.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % supply decrease by 4%:
%         A = A * 0.96; % unit: parking spaces.        
%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % supply increase by 4%:
%         A = A * 1.04; % unit: parking spaces.     
%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        t = 1/60; % unit: hours/60.
        v = 12.5; % unit: km/hr.
        kc = 20; % unit: veh/km.
        kj = 55; % unit: veh/km.
        Qmax = v * kc; %unit: veh/hour.

% 2nd part is the values that are provided by other scipts so that they can be
% easily adjusted (e.g., to another distribution). They include, 
% d0, n_demand.
% c2, beta^i; 
% c3, distribution of parking durations;
% d1, n'enter the area' (transition event and travel demand); 
% d2, l_ns/s; 
% d3, l_p/.
% d4, l_/; 

%---------PROCESS TO OBTAIN THE MATRIX:
% Individual scripts are written to generate colomns 8-12 (for the same time slice), based on colomns 1-7.
% Colomns 1-7 (for the next time slice) are generated within this scripts, which replies on colomns 8-12.
        
        %__________________________________________________________________
        % Build Matrix(the first row )                                    |
        matrix(1,1:22)=0;  %                                              |
        matrix(1,4)=A;%/                                                  |
        matrix(1,6)=v;%                                                   |
        [matrix(1,7),n_demand_VOT]= c1_input_n_demand(1);%/               |
        matrix(1,8)=matrix(1,7);%                                         |
        matrix(1,13)= matrix(1,7);%/                                      |
        matrix(1,14)=matrix(1,13);
        
        % demand matrix for VOT origins:
        matrix_demand_VOT(1,:) = n_demand_VOT'; 
        
        % initial price:
        parking_pricing(1,1) = c8_input_parking_price;
        % initial delta searching:
        delta_searching_veh_available_parking_spots = matrix(1,2)/matrix(1,4);
              
        % ns/s matrix for VOT origins:
        V = size(n_demand_VOT,1);
        matrix_ns_s_VOT(1,:) = zeros(1,V); 
        
        % initial values:
        matrix(1,3) = 183;
        
        % initial bias:
        bias = zeros(1440,1);
        bias(1,1) = (matrix(1,6)*t)/2 - matrix(1,6)*t;
        %_________________________________________________________________|
        
        
                %_________________________________________________________________________________
                % Build Matrix(From the second row on)
                i=1;                                                                            %|
                while i < 1440
%                while matrix(i,25) < matrix(i,17)-0.0000005                                     %|  
                % the loop keeps on,                                                            %|
                % when the cumulative no. of vehicles left the area < the cumulative no. of vehicles entered the area.
                
                    matrix(i+1,1)=matrix(i,1)+matrix(i,8)+matrix(i,11)-matrix(i,9)-matrix(i,12);%| Nns
                    matrix(i+1,2)=matrix(i,2)+matrix(i,9)-matrix(i,10);                         %| Ns
                    matrix(i+1,3)=matrix(i,3)+matrix(i,10)-matrix(i,11);                        %| Np
                    matrix(i+1,4)=A-matrix(i+1,3);                                              %| parking supply
                    matrix(i+1,5)=(matrix(i+1,1)+matrix(i+1,2))/L;                              %| density
                    if matrix(i+1,5)<=kc                                                        %| 
                        matrix(i+1,6)=v;                                                        %| speed
                    elseif matrix(i+1,5)<=kj                                                    %|
                        matrix(i+1,6)=Qmax/(kc-kj)*(1-kj/matrix(i+1,5));                        %| speed
                    else 
                        matrix(i+1,6)=0;                                                        %| speed
                    end                                                                         
                [matrix(i+1,7),n_demand_VOT]= c1_input_n_demand(i+1);                           %| n_demand
%               matrix(i+1,8)= matrix(i+1,7);                                                   %| /ns (old coding)   
                matrix(i+1,8)= d0_transitions_n_enter(matrix(i+1,7));                           %| /ns  
                
%               demand matrix for VOT origins:
                matrix_demand_VOT(i+1,:) = n_demand_VOT';                                       %| /ns for VOT, (i x V)-matrix 
                
                [matrix(i+1,9),n_ns_s_VOT] = d1_transitions_n_ns_s(matrix(1:i,8),matrix(1:i,6),t,matrix_demand_VOT);      %| ns/s

%               ns/s matrix for VOT origins:
                matrix_ns_s_VOT(i+1,:) = n_ns_s_VOT';                                           %| ns/s for VOT, (i x V)-matrix
                [matrix(i+1,10), guessed_price_vector(i,1), tau(i,1), E_p_vot(i,1), travel_distance_cost(i,1), penalty_distance(i,1), ...
                    total_costs(i,1), parking_pricing(i,1), delta_searching_veh_available_parking_spots(i,1)] = ...
                                d2_transitions_n_s_p(matrix(i+1,4),matrix(i+1,2),L,matrix(i+1,6),t,matrix_ns_s_VOT,matrix(1:i,4),...
                                matrix(1:i,2),matrix(1:i,6),parking_price_oscillation,A, parking_pricing, ...
                                delta_searching_veh_available_parking_spots);                                               %| s/p
                            
                            cost_travelling_next_space(i,1) = tau(i,1) * E_p_vot(i,1) + travel_distance_cost(i,1);                     
                            
                matrix(i+1,11)= d3_transitions_n_p_ns(matrix(1:i,10),t);                                                    %| p/ns
                matrix(i+1,12)= d4_transitions_n_ns_leave(matrix(1:i,8),matrix(1:i,11),matrix(1:i,6),t);                    %| ns/
                matrix(i+1,13)= matrix(i,13)+matrix(i+1,7);                                     %| cumulative n_demand
                matrix(i+1,14)= matrix(i,14)+matrix(i+1,8);                                     %| cumulative /ns
                matrix(i+1,15)= matrix(i,15)+matrix(i+1,9);                                     %| cumulative ns/s
                matrix(i+1,16)= matrix(i,16)+matrix(i+1,10);                                    %| cumulative s/p
                matrix(i+1,17)= matrix(i,17)+matrix(i+1,11);                                    %| cumulative p/ns
                matrix(i+1,18)= matrix(i,18)+matrix(i+1,12);                                    %| cumulative ns/
                matrix(i+1,19)= d5_transitions_n_ns_leave_parkers(matrix(1:i,8),matrix(1:i,11),matrix(1:i,6),t);    %| the number of parkers that leave the area in this slice
                matrix(i+1,20)= matrix(i,20)+matrix(i+1,19);                                                        %| cumulative_nleave_parkers
                matrix(i+1,21)= matrix(i+1,14)*(1-c2_input_beta)-matrix(i+1,15)+matrix(i+1,17)-matrix(i+1,20);      %| Nns_parkers
               %Nns_parkers   =  N enter   * park rate    -Nstarttosearch  + Ndepartparking-(N leave- Nleave_throughtraffic)           
                matrix(i+1,22)= matrix(i+1,14)*c2_input_beta-(matrix(i+1,18)-matrix(i+1,20));                       %| Nns_throughtraffic

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Determine bias for maximum driven distance
                bias(i+1,1) = (matrix(i+1,6)*t)/2 - matrix(i+1,6)*t;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                i = i + 1                                                                      %|
                if i == 1212
                    i
                end
                
                end                                                                            %|
                %________________________________________________________________________________
                
%                 matrix
%                 guessed_price_vector
%                 tau 
%                 E_p_vot
%                 penalty_distance
%                 total_costs
%                 parking_pricing

end

