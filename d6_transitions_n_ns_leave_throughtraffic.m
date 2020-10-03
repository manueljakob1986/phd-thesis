function nleave_throughtraffic=d6_transitions_n_ns_leave_throughtraffic(nenter, ndepart, speed, t)

[a b]=size(nenter);
nleave_throughtraffic=0;
for i=1:a
    distancetraveled(i,1)=sum(speed(i:a,1))*t;
    distancetraveled(i,2)=sum(speed(i:a-1,1))*t;
    if and(distancetraveled(i,1)>=c5_input__l_throughtraffic_leave,distancetraveled(i,2)<c5_input__l_throughtraffic_leave)
        nleave_throughtraffic=nleave_throughtraffic+nenter(i,1)*(c2_input_beta);
    end
end

end