function off_street_parking_pricing_switched_on = c16_input_switch_on_off_street_parking_pricing

% This input parameter switches on the demand-responsive off-street 
% parking pricing.
% In case it is switched off, the off-street parking pricing will be
% considered as constant.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demand-responsive off-street parking pricing = 
% switched on to demand-responsive 
% (optimize garage and on-street parking capacity),
% then: off_street_parking_pricing_switched_on = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demand-responsive off-street parking pricing = 
% switched to constant, then:
% off_street_parking_pricing_switched_on = 2
% And the off-street parking pricing is considered as constant.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demand-responsive off-street parking pricing = 
% switched off (not existent), then:
% off_street_parking_pricing_switched_on = 0
% And the off-street parking pricing is considered as zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
off_street_parking_pricing_switched_on = 2;

end