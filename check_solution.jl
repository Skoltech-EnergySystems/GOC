
println("################################################")
println("Check if the solution violates the accuracy threshold")
accur_tresh = 1e-6;
println("In Phase 0 allowed violatation of all constraints: ", accur_tresh)
println("Number of violated values in REACTIVE generation technical limits")
println(sum(.!(Qmin - accur_tresh.<=getvalue(q)[:,1].<= Qmax + accur_tresh)))
println(sum(.!(Qmin - accur_tresh.<=getvalue(q)[:,2].<= Qmax + accur_tresh)))
println("Number of violated values in ACTIVE generation technical limits")
println(sum(.!(Pmin - accur_tresh.<=getvalue(p)[:,1].<= Pmax + accur_tresh)))
println(sum(.!(Pmin - accur_tresh.<=getvalue(p)[:,2].<= Pmax + accur_tresh)))

##
println("Number of violated values related with VoltageLimits")
println(sum(.!(-pi/2 - accur_tresh.<= getvalue(δ) .<= pi/2 + accur_tresh)))
println(sum(.!(Vmin - accur_tresh.<= getvalue(V) .<= Vmax + accur_tresh)))

##
println("Number of violated values related with FlowLimits")
PLine = zeros(NBus,NBus,NK);
QLine = PLine;
S_ij = zeros(NBr,NK);
S_ji = zeros(NBr,NK);
Losses = zeros(NBr,NK);

for k=1:NK
    for l=1:NBr
        i = links[l][1];
        j = links[l][2];

        PLine[i,j,k] = (getvalue(V)[i,k]*getvalue(V)[j,k]*(GL[i,j]*cos(getvalue(δ)[i,k]-getvalue(δ)[j,k]) + BL[i,j]*sin(getvalue(δ)[i,k]-getvalue(δ)[j,k])) - getvalue(V)[i,k]^2*GL[i,j])*aL[k,l]
        PLine[j,i,k] = (getvalue(V)[j,k]*getvalue(V)[i,k]*(GL[j,i]*cos(getvalue(δ)[j,k]-getvalue(δ)[i,k]) + BL[j,i]*sin(getvalue(δ)[j,k]-getvalue(δ)[i,k])) - getvalue(V)[j,k]^2*GL[j,i])*aL[k,l]

        QLine[i,j,k] = (getvalue(V)[i,k]*getvalue(V)[j,k]*(GL[i,j]*sin(getvalue(δ)[i,k]-getvalue(δ)[j,k]) - BL[i,j]*cos(getvalue(δ)[i,k]-getvalue(δ)[j,k])) + getvalue(V)[i,k]^2*BL[i,j])*aL[k,l];
        QLine[j,i,k] = (getvalue(V)[j,k]*getvalue(V)[i,k]*(GL[j,i]*sin(getvalue(δ)[j,k]-getvalue(δ)[i,k]) - BL[j,i]*cos(getvalue(δ)[j,k]-getvalue(δ)[i,k])) + getvalue(V)[j,k]^2*BL[j,i])*aL[k,l];

        S_ij[l,k] = sqrt(PLine[i,j,k].^2 + QLine[i,j,k].^2);
        S_ji[l,k] = sqrt(PLine[j,i,k].^2 + QLine[j,i,k].^2);
        Losses[l,k] = abs(S_ij[l,k] + S_ji[l,k]);

    end
end

println(sum(.!(S_ij[:,1] .<= RATE_A + accur_tresh)))
println(sum(.!(S_ij[:,2] .<= RATE_A + accur_tresh)))
println(sum(.!(S_ji[:,1] .<= RATE_A + accur_tresh)))
println(sum(.!(S_ji[:,2] .<= RATE_A + accur_tresh)))

##
println("Number of violated values related with NodalBalance P")
UnBalanceP_list = [];
UnBalanceP = zeros(NBus,NK)
for i=1:NBus
    for k=1:NK
        if isdefined(:NGen)
            for g=1:NGen
                UnBalanceP[i,k] = sum(Ω[i,g]*getvalue(p)[g,k] for g=1:NGen) - Pd[i] - getvalue(V)[i,k]*sum(getvalue(V)[j,k]*(GL[i,j]*cos(getvalue(δ)[i,k]-getvalue(δ)[j,k]) + BL[i,j]*sin(getvalue(δ)[i,k]-getvalue(δ)[j,k])) for j=1:NBus);
            end
        else
            UnBalanceP[i,k] =  getvalue(p)[i,k] - Pd[i] - getvalue(V)[i,k]*sum(getvalue(V)[j,k]*(GL[i,j]*cos(getvalue(δ)[i,k]-getvalue(δ)[j,k]) + BL[i,j]*sin(getvalue(δ)[i,k]-getvalue(δ)[j,k])) for j=1:NBus);
        end
        UnBalanceP_val = (.!(-accur_tresh.<=UnBalanceP[i,k].<= accur_tresh))
        UnBalanceP_list = [UnBalanceP_list; UnBalanceP_val]
    end
end
println(sum(UnBalanceP_list))
#println(UnBalanceP)
if sum(abs.(UnBalanceP_list)) > 0
    println("MAX UnBalanceP: ", maximum(abs.(UnBalanceP)))
end

println("Number of violated values related with NodalBalance Q")
UnBalanceQ_list = [];
UnBalanceQ = zeros(NBus,NK)
for i=1:NBus
    for k=1:NK
        if isdefined(:NGen)
            for g=1:NGen
                UnBalanceQ[i,k] = sum(Ω[i,g]*getvalue(q)[g,k] for g=1:NGen) - Qd[i] - getvalue(V)[i,k]*sum(getvalue(V)[j,k]*(GL[i,j]*sin(getvalue(δ)[i,k]-getvalue(δ)[j,k]) - BL[i,j]*cos(getvalue(δ)[i,k]-getvalue(δ)[j,k])) for j=1:NBus);
            end
        else
            UnBalanceQ[i,k] =  getvalue(q)[i,k] - Qd[i] - getvalue(V)[i,k]*sum(getvalue(V)[j,k]*(GL[i,j]*sin(getvalue(δ)[i,k]-getvalue(δ)[j,k]) - BL[i,j]*cos(getvalue(δ)[i,k]-getvalue(δ)[j,k])) for j=1:NBus);
        end
        UnBalanceQ_val = (.!(-accur_tresh.<=UnBalanceQ[i,k].<= accur_tresh))
        UnBalanceQ_list = [UnBalanceQ_list; UnBalanceQ_val]
    end
end
println(sum(UnBalanceQ_list))
#println(UnBalanceQ)
if sum(abs.(UnBalanceQ_list)) > 0
    println("MAX UnBalanceQ: ", maximum(abs.(UnBalanceQ)))
end

println("Number of violated PV-PQ: qmax")
UnBalanceGmax_list = [];
UnBalanceGmax = zeros(NGen);
# Paired external and internal location of geberators
Bus_Gen = hcat(collect(gData[k].ID_ext for k in sort(collect(keys(gData)))),collect(gData[k].ID_int for k in sort(collect(keys(gData)))));

for k = 1:size(Bus_Gen,1)
    UnBalanceGmax[k] = min(max( 0,getvalue(V[Bus_Gen[k,1],1]) - getvalue(V[Bus_Gen[k,1],2]) ), max( 0,Qmax[Bus_Gen[k,2]] - getvalue(q[Bus_Gen[k,2],2]) ) )
    UnBalanceGmax_val = (.!(-accur_tresh.<=UnBalanceGmax[k].<= accur_tresh))
    UnBalanceGmax_list = [UnBalanceGmax_list; UnBalanceGmax_val]
end
println(sum(UnBalanceGmax_list))
if sum(abs.(UnBalanceGmax_list)) > 0
    println("MAX UnBalanceGmax: ", maximum(abs.(UnBalanceGmax)))
end

println("Number of violated PV-PQ: qmin")
UnBalanceGmin_list = [];
UnBalanceGmin = zeros(NGen);

for k = 1:size(Bus_Gen,1)
    UnBalanceGmin[k] = min(max( 0,getvalue(V[Bus_Gen[k,1],2]) - getvalue(V[Bus_Gen[k,1],1]) ), max( 0,getvalue(q[Bus_Gen[k,2],2])  - Qmin[Bus_Gen[k,2]]) )
    UnBalanceGmin_val = (.!(-accur_tresh.<=UnBalanceGmin[k].<= accur_tresh))
    UnBalanceGmin_list = [UnBalanceGmin_list; UnBalanceGmin_val]
end
println(sum(UnBalanceGmin_list))
if sum(abs.(UnBalanceGmin_list)) > 0
    println("MAX UnBalanceGmin: ", maximum(abs.(UnBalanceGmin)))
end
