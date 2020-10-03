function [ACT] = d15_avg_cruising_time(starttosearch_3, decideforparking_3, n_vehicles_s_dgp_VOT, n_dgp_searching_VOT, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJA: Changed in June 2017 to include parking garages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see coding in:
% Average cruising time (see: [...] = d8_parking_pricing_probability(...))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = size(starttosearch_3,1);

n_dgp_searching = sum(n_dgp_searching_VOT,2);
n_vehicles_s_dgp = sum(n_vehicles_s_dgp_VOT,2);

% Average cruising time
for i_minus_CT_max = 1:a
    if sum(starttosearch_3(1:i_minus_CT_max) + n_dgp_searching(1:i_minus_CT_max),1) > sum(decideforparking_3(1:a) + n_vehicles_s_dgp(1:a),1)
        break
    end
end
i_minus_CT_max = i_minus_CT_max - 1;
CT_max = a - i_minus_CT_max;
ACT = CT_max/2;

ACT = ACT * t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end