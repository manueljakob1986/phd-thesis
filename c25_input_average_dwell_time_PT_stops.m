function average_dwell_time_PT_stops = c25_input_average_dwell_time_PT_stops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This input parameter describes the average dwell time , i.e., the average
% time each PT vehicle spends at PT stops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
average_dwell_time_PT_stops = c24_input_number_of_PT_stops * (30/60)/60; %30 seconds per PT stop

end