function [total_revenue, total_on_street_revenue, total_garage_revenue] = d18_total_revenue(depart_parking_on_street, ndepart_garage)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in January 2018 to compute total revenue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of total revenue.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

on_street_parking_fee = h2_getGlobal_initial_parking_pricing;

garage_parking_fee = h2_getGlobal_initial_garage_parking_pricing;

duration_on_street_parking = ...
    (25*129 + 75*299 + 125*269 + 175*196 + 225*147 + 275*49)/...
    (129 + 299 + 269 + 196 + 147 + 49);
% average_parking_durations_on_street = 128.6731
% average_parking_durations_on_street = 2.1446 (in hours)

duration_garage_parking = ...
    (25*77 + 75*179 + 125*162 + 175*118 + 225*88 + 275*145 + 325*130 + ...
    375*91 + 425*83 + 475*105 + 525*76 + 575*71 + 625*41 + 675*29 + ...
    725*21 + 775*8 + 825*6 + 875*5 + 925*5 + 975*1 + 1025*2 + 1075*1 + 1125*1)/...
    (77 + 179 + 162 + 118 + 88 + 145 + 130 + ...
    91 + 83 + 105 + 76 + 71 + 41 + 29 + ...
    21 + 8 + 6 + 5 + 5 + 1 + 2 + 1 + 1);
% average_parking_durations_garage = 307.2491
% average_parking_durations_garage = 5.1208 (in hours)

total_revenue = sum(depart_parking_on_street(:,:),1)*on_street_parking_fee*(duration_on_street_parking/60) + ...
    sum(ndepart_garage(:,:),1)*garage_parking_fee*(duration_garage_parking/60);

total_on_street_revenue = sum(depart_parking_on_street(:,:),1)*on_street_parking_fee*(duration_on_street_parking/60);

total_garage_revenue = sum(ndepart_garage(:,:),1)*garage_parking_fee*(duration_garage_parking/60);

end