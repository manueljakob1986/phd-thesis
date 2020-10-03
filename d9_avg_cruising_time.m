function [ACT] = d9_avg_cruising_time(starttosearch_3, decideforparking_3, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please see coding in:
% Average cruising time (see: [...] = d8_parking_pricing_probability(...))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = size(starttosearch_3,1);

% Average cruising time
for i_minus_CT_max = 1:a
    if sum(starttosearch_3(1:i_minus_CT_max),1) > sum(decideforparking_3(1:a),1)
        break
    end
end
i_minus_CT_max = i_minus_CT_max - 1;
CT_max = a - i_minus_CT_max;
ACT = CT_max/2;

ACT = ACT * t;

end