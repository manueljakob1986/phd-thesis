function findparking_3=d2_transitions_n_s_p_new(availableparking,Ns,v)
global t L;

if Ns<=availableparking
    if v*t<=L/Ns
        findparking_3=Ns*(1-(1-v*t/L)^availableparking);
    elseif v*t<=L
        findparking_3=Ns*(1+(1-1/Ns)^availableparking*log(v*t/L)/log(Ns));
    else
        findparking_3=Ns;
    end
else
    if v*t<=L/Ns
        findparking_3=Ns*(1-(1-v*t/L)^availableparking);
    elseif v*t<=L*availableparking/Ns
        findparking_3=Ns*(availableparking/Ns+(availableparking/Ns-1+(1-1/Ns)^availableparking)*log(Ns/availableparking*v*t/L)/log(availableparking));
    else
        findparking_3=availableparking;
    end
end
end