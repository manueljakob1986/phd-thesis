function [p, gp_a2, gp_b2] = c20_input_park_ride_parkingduration(td, t)

% Garage parking durations:
gp_a2 = 2.1323;
gp_b2 = 137.4145;

% No time control
timecontrol_gp = 1440;
           
% % when there is no time control, use this:
% p = gamcdf(td,a2,b2) - gamcdf(td-60*t,a2,b2);

% % Consider whole day: 1440 minutes
    if td <= timecontrol_gp
       p = gamcdf(td,gp_a2,gp_b2) - gamcdf(td-60*t,gp_a2,gp_b2);   
    else
       p = 0;
    end

end