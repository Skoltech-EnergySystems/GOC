function dPrint(PSCOPF, genSeg, BLGS, Br, contingency, baseMVA)
    # Extracting calculated generation values for base case
    Pgen = [];
    for y in getvalue(PSCOPF[:p])[:,1]
        if abs(y) > 0.0001
            push!(Pgen,y*100)
        end
    end

    Qgen = [];
    for z in getvalue(PSCOPF[:q])[:,1]
        if abs(z) > 0.0001
            push!(Qgen,z*100)
        end
    end

    #------------- Defining structure of the output file "solution1.txt" ----------
    sizePg = size(genSeg)[1];
    pr = Array{Any}(sizePg,4);
    # Bus ID
    pr[1:sizePg,1] = genSeg[:,1];
    # Unit ID
    pr[1:sizePg,2] = repeat(["'1 '"],outer=[sizePg]);
    # real power in megawatts
    pr[1:sizePg,3] = Pgen;
    # reactive power in megaVar
    pr[1:sizePg,4] = Qgen;

    open("solution1.txt", "w") do f1
        write(f1, "--generation dispatch \n")
        write(f1, "bus id,unit id,pg(MW),qg(MVar) \n")
        for i in 1:sizePg
            write(f1, join(pr[i,:], ","))
            write(f1, " \n")
        end
        write(f1, "--end of generation dispatch")
    end

    # Extracting calculated generation values for contingency case
    Pgenk = [];
    for y in getvalue(PSCOPF[:p])[:,2]
        if abs(y) > 0.0001
            push!(Pgenk,y*100)
        end
    end

    Qgenk = [];
    for z in getvalue(PSCOPF[:q])[:,2]
        if abs(z) > 0.0001
            push!(Qgenk,z*100)
        end
    end

    #------- Defining structure of the output file "solution2.txt" ---------------
    pr2 = Array{Any}(sizePg,5);

    ## 1. contingency generator dispatch
    # contingency ID
    pr2[1:sizePg,1] = repeat(["1"],outer=[sizePg]);
    # generator ID (not used in the evaluation process)
    for i in 1:sizePg
        pr2[i,2] = "l_$i"
    end
    # bus ID
    pr2[1:sizePg,3] = genSeg[:,1];
    # unit ID
    pr2[1:sizePg,4] = repeat(["'1 '"],outer=[sizePg]);
    # Reactive power in megaVar
    pr2[1:sizePg,5] = Qgenk;

    # Number of buses
    NBus = size(BLGS,1);
    # Number of branches (power lines and transformers)
    NBr = size(Br,1);

    # All lines and transformers are online
    aL = transpose(ones(size(Br,1),1));
    # Contingency FROM buses
    Cont_FR = float(contingency[:,3])[1,1];
    # Contingency TO buses
    Cont_TO = float(contingency[:,4])[1,1];
    # Number of contingencies
    Cont_size = size(contingency,1);
    # Contingency matrix
    Br_ind=hcat(collect(1:size(Br,1)),Br[:,1:2]);
    aL_loop = transpose(ones(size(Br,1),1))
    for i=1:size(Br,1)
        if (Br_ind[i,2]==Cont_FR && Br_ind[i,3]==Cont_TO)
            aL_loop[1,i]=0;
            aL=vcat(aL,aL_loop)
        end
    end

    # Number of contingencies
    NK = size(aL,1);

    ## 2. contingency bus information
    pr3 = Array{Any}(2*NBus,4);
    # contingency ID
    pr3[1:NBus,1] = repeat(["0"],outer=[NBus]); # 0. base case
    pr3[NBus+1:2*NBus,1] = repeat(["1"],outer=[NBus]); # 1. contingency case
    # bus ID
    pr3[1:NBus,2] = round.(Int,BLGS[:,1]); # 0. base case
    pr3[NBus+1:2*NBus,2] = round.(Int,BLGS[:,1]); # 1. contingency case
    # Voltage in per unit
    pr3[1:NBus,3] = getvalue(PSCOPF[:V])[:,1]; # 0. base case
    pr3[NBus+1:2*NBus,3] = getvalue(PSCOPF[:V])[:,2]; # 1. contingency case
    # Voltage angle in degree
    pr3[1:NBus,4] = rad2deg.(getvalue(PSCOPF[:δ]))[:,1]; # 0. base case
    pr3[NBus+1:2*NBus,4] = rad2deg.(getvalue(PSCOPF[:δ]))[:,2]; # 1. contingency case

    ## 3. contingency delta
    pr4 = Array{Any}(NK-1,2);
    pr4[1,1] = '1';
    pr4[1,2] = getvalue(PSCOPF[:Delta])*100;

    ## 4. contingency line flow information
    pr5 = Array{Any}(2*NBr,9)
    # contingency ID
    pr5[1:NBr,1] = repeat(["0"],outer=[NBr]); # 0. base case
    pr5[NBr+1:2*NBr,1] = repeat(["1"],outer=[NBr]); # 1. contingency case
    # line ID (not used in the evaluation process)
    for i in 1:NBr
        pr5[i,2] = "i_$i" # 0. base case
    end

    for i in 1:NBr
        pr5[NBr+i,2] = "i_$i" # 1. contingency case
    end
    # origin bus ID
    pr5[1:NBr,3] = round.(Int,Br[:,1]); # 0. base case
    pr5[NBr+1:2*NBr,3] = round.(Int,Br[:,1]); # 1. contingency case
    # destination bus ID
    pr5[1:NBr,4] = round.(Int,Br[:,2]); # 0. base case
    pr5[NBr+1:2*NBr,4] = round.(Int,Br[:,2]); # 1. contingency case
    # circuit ID
    pr5[1:2*NBr,5] = repeat(["'BL'"],outer=[2*NBr]);
    # real power in megawatts at origin
    originP = [];
    for (x,y) in zip(round.(Int,Br[:,1]), round.(Int,Br[:,2]))
        push!(originP,getvalue(PSCOPF[:pl])[x,y,1]*baseMVA) # 0. base case
    end
    pr5[1:NBr,6] = originP;

    originP = [];
    for (x,y) in zip(round.(Int,Br[:,1]), round.(Int,Br[:,2]))
        push!(originP,getvalue(PSCOPF[:pl])[x,y,2]*baseMVA) # 1. contingency case
    end
    pr5[NBr+1:2*NBr,6] = originP;
    # reactive power in MVar at origin
    originQ = [];
    for (x,y) in zip(round.(Int,Br[:,1]), round.(Int,Br[:,2]))
        push!(originQ,getvalue(PSCOPF[:ql])[x,y,1]*baseMVA) # 0. base case
    end
    pr5[1:NBr,7] = originQ;

    originQ = [];
    for (x,y) in zip(round.(Int,Br[:,1]), round.(Int,Br[:,2]))
        push!(originQ,getvalue(PSCOPF[:ql])[x,y,2]*baseMVA) # 1. contingency case
    end
    pr5[NBr+1:2*NBr,7] = originQ;

    # real power in megawatts at destination
    destP = [];
    for (x,y) in zip(round.(Int,Br[:,2]), round.(Int,Br[:,1]))
        push!(destP,getvalue(PSCOPF[:pl])[x,y,1]*baseMVA) # 0. base case
    end
    pr5[1:NBr,8] = destP;

    destP = [];
    for (x,y) in zip(round.(Int,Br[:,2]), round.(Int,Br[:,1]))
        push!(destP,getvalue(PSCOPF[:pl])[x,y,2]*baseMVA) # 1. contingency case
    end
    pr5[NBr+1:2*NBr,8] = destP;

    # reactive power in MVar at destination
    destQ = [];
    for (x,y) in zip(round.(Int,Br[:,2]), round.(Int,Br[:,1]))
        push!(destQ,getvalue(PSCOPF[:ql])[x,y,1]*baseMVA) # 0. base case
    end
    pr5[1:NBr,9] = destQ;

    destQ = [];
    for (x,y) in zip(round.(Int,Br[:,2]), round.(Int,Br[:,1]))
        push!(destQ,getvalue(PSCOPF[:ql])[x,y,2]*baseMVA) # 1. contingency case
    end
    pr5[NBr+1:2*NBr,9] = destQ;
    # Deleting row of contingency case
    for i in size(pr5)[1]
        if pr5[i,1] == "1"
            for j in 1:NBr
                if aL[2,j] == 0
                    pr5 = pr5[setdiff(1:end,j+NBr),:]
                end
            end
        end
    end

    # Writing into the file solution2.txt
    open("solution2.txt", "w") do f2
        write(f2, "--contingency generator \n")
        write(f2, "contingency id,genID,bus id,unit id,q(MW) \n")
        for i in 1:size(pr2)[1]
            write(f2, join(pr2[i,:], ","))
            write(f2, " \n")
        end
        write(f2, "--end of contingency generator \n")
        write(f2, "--bus \n")
        write(f2, "contingency id,bus id,v(pu),theta(deg) \n")
        for i in 1:size(pr3)[1]
            write(f2, join(pr3[i,:], ","))
            write(f2, " \n")
        end
        write(f2, "--end of bus \n")
        write(f2, "--Delta \n")
        write(f2, "contingency id,Delta(MW) \n")
        for i in 1:size(pr4)[1]
            write(f2, join(pr4[i,:], ","))
            write(f2, " \n")
        end
        write(f2, "--end of Delta \n")
        write(f2, "--line flow \n")
        write(f2, "contingency id,line id,origin bus id,destination bus id,circuit id, p_origin(MW) q_origin(MVar), p_destination(MW), q_destination(MVar) \n")
        for i in 1:size(pr5)[1]
            write(f2, join(pr5[i,:], ","))
            write(f2, " \n")
        end
        write(f2, "--end of line flow")
    end
#    return f1, f2
end
