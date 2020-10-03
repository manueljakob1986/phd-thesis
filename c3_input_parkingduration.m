function p=c3_input_parkingduration(td, t, on_or_off_street_parking)

[a2, b2, on_a2, on_b2, gp_a2, gp_b2] = c18_input_parking_duration_distribution_parameter;

pdf_gammamixture = @(duration_garage,p,shape1,shape2,scale1,scale2) ...
                         p*gamcdf(duration_garage,shape1,scale1) + (1-p)*gamcdf(duration_garage,shape2,scale2);              

% % when there is no time control, use this:
% p = gamcdf(td,a2,b2) - gamcdf(td-60*t,a2,b2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add parking duration control (MJA, June 2017)
% when there is a duration time control, use below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% on-street parking:
if on_or_off_street_parking == 1

% timecontrol_op = c14_input_on_street_parking_duration;
timecontrol_gp = c17_input_off_street_parking_duration;

% % Consider whole day: 1440 minutes
    if td <= timecontrol_gp
       p = gamcdf(td,on_a2,on_b2) - gamcdf(td-60*t,on_a2,on_b2);
%        p = gamcdf(td,a2,b2) - gamcdf(td-60*t,a2,b2);
    else
       p = 0;
    end

% garage parking:
elseif on_or_off_street_parking == 2

timecontrol_gp = c17_input_off_street_parking_duration;

% % Consider whole day: 1440 minutes
    if td <= timecontrol_gp
       p = gamcdf(td,gp_a2,gp_b2) - gamcdf(td-60*t,gp_a2,gp_b2);
%        p = gamcdf(td,a2,b2) - gamcdf(td-60*t,a2,b2);  
      
%        p = pdf_gammamixture(td, 0.4185, 5.1378, 2.0146, 79.4054, 104.3817) - ...
%            pdf_gammamixture(td-60*t, 0.4185, 5.1378, 2.0146, 79.4054, 104.3817);
       
    else
       p = 0;
    end

% % garage parking:
% elseif on_or_off_street_parking == 2
% 
% timecontrol_gp = c17_input_off_street_parking_duration;
% 
% % % Consider whole day: 1440 minutes
%     if td <= timecontrol_gp
%        p = gamcdf(td,gp_a2,gp_b2) - gamcdf(td-60*t,gp_a2,gp_b2);
%     else
%        p = 0;
%     end
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % on-street parking:
% if on_or_off_street_parking == 1
% 
% timecontrol_op = c14_input_on_street_parking_duration;
% 
%     if td < timecontrol_op
%        p = gamcdf(td,a2,b2) - gamcdf(td-60*t,a2,b2);
%     elseif td == timecontrol_op
%        p = 1-gamcdf(td-60*t,a2,b2);
%     else
%        p = 0;
%     end
% 
% % garage parking:
% elseif on_or_off_street_parking == 2
% 
% timecontrol_gp = c17_input_off_street_parking_duration;
% 
%     if td < timecontrol_gp
%        p = gamcdf(td,a2,b2) - gamcdf(td-60*t,a2,b2); 
%     elseif td == timecontrol_gp
%        p = 1 - gamcdf(td-60*t,a2,b2); 
%     else    
%        p = 0;
%     end
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end