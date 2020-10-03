function guessed_pricing_agg_time_steps = c12_input_guessed_pricing_agg_time_steps

% This factors is essential for the appoximation of the guessed parking
% pricing in the future time step i+1. In case an aggregation method
% approach is used to compute the guessed parking pricing vector, then we
% use a fixed number of last time interval steps, that should be used to
% aggregate the guessed parking price.
% This is more reasonable than having a percentage value, because for a
% longer simulation time, too many past time interval steps could be
% considered even if they should not influence the result.

% This input parameter returns the fixed number of time interval steps,
% that should be considered in the aggregation of the guessed parking
% pricing vector:
 
guessed_pricing_agg_time_steps = 10;

end