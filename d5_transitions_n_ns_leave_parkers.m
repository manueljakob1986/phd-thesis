function nleave_parkers=d5_transitions_n_ns_leave_parkers(nenter, ndepart, speed, t)

[a b]=size(nenter);
nleave_parkers=0;
for i=1:a
    distancetraveled(i,1)=sum(speed(i:a,1))*t;
    distancetraveled(i,2)=sum(speed(i:a-1,1))*t;
    if and(distancetraveled(i,1)>=c6_input__l_p_leave,distancetraveled(i,2)<c6_input__l_p_leave)
        nleave_parkers=nleave_parkers+ndepart(i,1);
    end
end

end