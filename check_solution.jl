
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
        for g=1:NGen
            UnBalanceP[i,k] =  Ω[i,g]*getvalue(p)[g,k] - Pd[i] - getvalue(V)[i,k]*sum(getvalue(V)[j,k]*(GL[i,j]*cos(getvalue(δ)[i,k]-getvalue(δ)[j,k]) + BL[i,j]*sin(getvalue(δ)[i,k]-getvalue(δ)[j,k])) for j=1:NBus);
        end
        UnBalanceP_val = (.!(-accur_tresh.<=UnBalanceP[i,k].<= accur_tresh))
        UnBalanceP_list = [UnBalanceP_list; UnBalanceP_val]
    end
end
println(sum(UnBalanceP_list))
#println(UnBalanceP)

println("Number of violated values related with NodalBalance Q")
UnBalanceQ_list = [];
UnBalanceQ = zeros(NBus,NK)
for i=1:NBus
    for k=1:NK
        for g=1:NGen
            UnBalanceQ[i,k] =  Ω[i,g]*getvalue(q)[g,k] - Qd[i] - getvalue(V)[i,k]*sum(getvalue(V)[j,k]*(GL[i,j]*sin(getvalue(δ)[i,k]-getvalue(δ)[j,k]) - BL[i,j]*cos(getvalue(δ)[i,k]-getvalue(δ)[j,k])) for j=1:NBus);
        end
        UnBalanceQ_val = (.!(-accur_tresh.<=UnBalanceQ[i,k].<= accur_tresh))
        UnBalanceQ_list = [UnBalanceQ_list; UnBalanceQ_val]
    end
end
println(sum(UnBalanceQ_list))
#println(UnBalanceQ)
