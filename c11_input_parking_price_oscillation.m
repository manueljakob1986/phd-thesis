function dumping_coefficient = c11_input_parking_price_oscillation(new_dumping_coefficient)

% This factors is the oscillation reduction parameter for the parking
% pricing computation.
% parking_price_oscillation denominator smaller -> higher reduction of
% oscillations for parking pricing
% parking_price_oscillation denominator higher -> smaller reduction of
% oscillations for parking pricing

% Assumption: Maximum parking price is equal for all on-street parking spots:
 
dumping_coefficient = new_dumping_coefficient;

end