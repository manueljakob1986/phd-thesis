% function ndepart=d3_transitions_n_p_ns(naccess,t,parking_duration_a_shape,parking_duration_b_scale)
% % naccess=[0;0;0;0;1;2;2;3;5;6;];
% % t=1/60;
% [a b]=size(naccess);
% ndepart=0;
% for i=1:a
%    td=(a+1-i)*60*t;% unit: minutes;
%    n(i,1)=naccess(i,1)*c3_input_parkingduration(td,t,parking_duration_a_shape,parking_duration_b_scale);
% % n(i,1) is the number of vehicles who arrived in time slice i and leaving
% % in time slice a, in other words, the cars who arrived in time slice i and
% % stayed for a period of td.
%    ndepart=ndepart+n(i,1);
% end
% end