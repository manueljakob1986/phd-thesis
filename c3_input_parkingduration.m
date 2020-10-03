function p=c3_input_parkingduration(td,t)
a2 = 1.6;
b2 = 142;
%the average parking duration=a2*b2

% % when there is no time control, use this:
% p = gamcdf(td,a2,b2) - gamcdf(td-60*t,a2,b2);

% when there is a duration time control, use below:

% timecontrol = 60; % the longerst parking duration is 60 minutes.
timecontrol = c14_input_on_street_parking_duration;

    if td < timecontrol
       p = gamcdf(td,a2,b2) - gamcdf(td-60*t,a2,b2);
    elseif td == timecontrol
       p = 1-gamcdf(td-60*t,a2,b2);
    else
       p = 0;
    end

end