function dPrint(PSCOPF,NGen,NBus,NBr,NBrc,NK,aL,Ct,Ctc,Cf,Cfc,Yf,Yfc,Yt,Ytc,br_int,br_intc,genSeg,busSeg,brList,brDList,baseMVA)
    println("Deriving FROM and TO injections")
    tic()
    # Computed V in complex form for pre-contingency case
    Vk0 = getvalue(PSCOPF[:V])[:,1];
    δk0 = getvalue(PSCOPF[:δ])[:,1];
    Vk0c = Vk0.*exp.(im*δk0);

    # Computed V in complex form for post-contingency case
    Vk1 = getvalue(PSCOPF[:V])[:,2];
    δk1 = getvalue(PSCOPF[:δ])[:,2];
    Vk1c = Vk1.*exp.(im*δk1);

    # Apparent power FROM for pre-contingency case
    Sf_k0 = Diagonal(Cf*Vk0c)*conj(Yf)*conj(Vk0c);
    # Apparent power FROM for post-contingency case
    Sf_k1 = Diagonal(Cfc*Vk1c)*conj(Yfc)*conj(Vk1c);

    # Apparent power TO for pre-contingency case
    St_k0 = Diagonal(Ct*Vk0c)*conj(Yt)*conj(Vk0c);
    # Apparent power TO for post-contingency case
    St_k1 = Diagonal(Ctc*Vk1c)*conj(Ytc)*conj(Vk1c);

    # Active power FROM for pre-contingency case
    Pf_k0 = real(Sf_k0);
    # Active power FROM for post-contingency case
    Pf_k1 = real(Sf_k1);

    # Reactive power FROM for pre-contingency case
    Qf_k0 = imag(Sf_k0);
    # Real power FROM for post-contingency case
    Qf_k1 = imag(Sf_k1);

    # Active power TO for pre-contingency case
    Pt_k0 = real(St_k0);
    # Active power TO for post-contingency case
    Pt_k1 = real(St_k1);

    # Reactive power TO for pre-contingency case
    Qt_k0 = imag(St_k0);
    # Real power TO for post-contingency case
    Qt_k1 = imag(St_k1);

    toc()
    ############################ OUTPUT FILES #####################################
    ################################## SOLUTION1.TXT ##############################
    println("Writing results")
    tic()
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
    pr = Array{Any}(NGen,4);
    # Bus ID
    pr[1:NGen,1] = genSeg[:,1];
    # Unit
    pr[1:NGen,2] = repeat(["'1 '"],outer=[NGen]);
    # real power in megawatts
    pr[1:NGen,3] = Pgen;
    # reactive power in megaVar
    pr[1:NGen,4] = Qgen;

    # Writing into the file "solution1.txt"
    open("solution1.txt", "w") do f1
        write(f1, "--generation dispatch\r\n")
        write(f1, "bus id,unit id,pg(MW),qg(MVar)\r\n")
        for i in 1:NGen
            write(f1, join(pr[i,:], ","))
            write(f1, "\r\n")
        end
        write(f1, "--end of generation dispatch")
    end
    #
    ################################## SOLUTION2.TXT ##############################
    # Extracting calculated generation values for contingency case
    Qgenk = [];
    for z in getvalue(PSCOPF[:q])[:,2]
        if abs(z) > 0.0001
            push!(Qgenk,z*100)
        end
    end

    #------- Defining structure of the output file "solution2.txt" ---------------
    pr2 = Array{Any}(NGen,5);

    ## 1. contingency generator dispatch
    # contingency ID
    pr2[1:NGen,1] = repeat(["1"],outer=[NGen]);
    # generator ID (not used in the evaluation process)
    for i in 1:NGen
        pr2[i,2] = "l_$i"
    end
    # bus ID
    pr2[1:NGen,3] = genSeg[:,1];
    # unit ID
    pr2[1:NGen,4] = repeat(["'1 '"],outer=[NGen]);
    # Reactive power in megaVar
    pr2[1:NGen,5] = Qgenk;

    ## 2. contingency bus information
    pr3 = Array{Any}(2*NBus,4);
    # contingency ID
    pr3[1:NBus,1] = repeat(["0"],outer=[NBus]); # 0. base case
    pr3[NBus+1:2*NBus,1] = repeat(["1"],outer=[NBus]); # 1. contingency case
    # bus ID
    pr3[1:NBus,2] = busSeg[:,1]; # 0. base case
    pr3[NBus+1:2*NBus,2] = busSeg[:,1]; # 1. contingency case
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
    pr5 = Array{Any}(2*NBr-1,9)
    # contingency ID
    pr5[1:NBr,1] = repeat(["0"],outer=[NBr]); # 0. base case
    pr5[NBr+1:2*NBr-1,1] = repeat(["1"],outer=[NBrc]); # 1. contingency case
    # line ID (not used in the evaluation process)
    for i in 1:NBr
        pr5[i,2] = "i_$i" # 0. base case
    end

    for i in 1:NBrc
        pr5[NBr+i,2] = "i_$i" # 1. contingency case
    end

    # List of origin buses
    #orig_bus = collect(brData[k].From for k in sort(collect(keys(brDList))));

    # List of destination buses
    #dest_bus = collect(brData[k].To for k in sort(collect(keys(brDList))));

    # origin bus ID
    pr5[1:NBr,3] = br_int[:,1]; # 0. base case
    pr5[NBr+1:2*NBr-1,3] = br_intc[:,1]; # 1. contingency case

    # destination bus ID
    pr5[1:NBr,4] = br_int[:,2]; # 0. base case
    pr5[NBr+1:2*NBr-1,4] = br_intc[:,2]; # 1. contingency case

    # circuit ID
    pr5[1:2*NBr-1,5] = repeat(["'BL'"],outer=[2*NBr-1]);

    # real power in megawatts at origin
    #=
    originP0 = [];
    for (x,y) in zip(orig_bus, dest_bus)
        push!(originP0,getvalue(pl)[x,y,1]*baseMVA) # 0. base case
    end
    =#
    pr5[1:NBr,6] = Pf_k0*100;
    #maximum(originP0 - Pf_k0*100)
    #=
    originP1 = [];|
    for (x,y) in zip(orig_bus, dest_bus)
        push!(originP1,getvalue(pl)[x,y,2]*baseMVA) # 1. contingency case
    end
    =#
    pr5[NBr+1:2*NBr-1,6] = Pf_k1*100
    #maximum(originP1 - Pf_k1*100)
    # reactive power in MVar at origin
    #=
    originQ0 = [];
    for (x,y) in zip(orig_bus, dest_bus)
        push!(originQ0,getvalue(ql)[x,y,1]*baseMVA) # 0. base case
    end
    =#
    pr5[1:NBr,7] = Qf_k0*100;
    #maximum(originQ0 - Qf_k0*100)
    #=
    originQ1 = [];
    for (x,y) in zip(orig_bus, dest_bus)
        push!(originQ1,getvalue(ql)[x,y,2]*baseMVA) # 1. contingency case
    end
    =#
    pr5[NBr+1:2*NBr-1,7] = Qf_k1*100
    #maximum(originQ1 - Qf_k1*100)

    # real power in megawatts at destination
    #=
    destP0 = [];
    for (x,y) in zip(dest_bus, orig_bus)
        push!(destP0,getvalue(pl)[x,y,1]*baseMVA) # 0. base case
    end
    =#
    pr5[1:NBr,8] = Pt_k0*100;
    #maximum(destP0 - Pt_k0*100)
    #=
    destP1 = [];
    for (x,y) in zip(dest_bus, orig_bus)
        push!(destP1,getvalue(pl)[x,y,2]*baseMVA) # 1. contingency case
    end
    =#
    pr5[NBr+1:2*NBr-1,8] = Pt_k1*100
    #maximum(destP1 - Pt_k1*100)

    # reactive power in MVar at destination
    #=
    destQ0 = [];
    for (x,y) in zip(dest_bus, orig_bus)
        push!(destQ0,getvalue(ql)[x,y,1]*baseMVA) # 0. base case
    end
    =#
    pr5[1:NBr,9] = Qt_k0*100;
    #maximum(destQ0 - Qt_k0*100)
    #=
    destQ1 = [];
    for (x,y) in zip(dest_bus, orig_bus)
        push!(destQ1,getvalue(ql)[x,y,2]*baseMVA) # 1. contingency case
    end
    =#
    pr5[NBr+1:2*NBr-1,9] = Qt_k1*100
    #maximum(destQ1 - Qt_k1*100)
    # Deleting row of contingency case
    #=
    for i in size(pr5)[1]
        if pr5[i,1] == "1"
            for j in 1:NBr
                if aL[2,j] == 0
                    pr5 = pr5[1:size(pr5,1) .!= j+NBr,: ]
                end
            end
        end
    end
    #
    =#
    pr5
    # Writing into the file solution2.txt
    open("solution2.txt", "w") do f2
        write(f2, "--contingency generator\r\n")
        write(f2, "conID,genID,busID,unitID,q(MVar)\r\n")
        for i in 1:size(pr2)[1]
            write(f2, join(pr2[i,:], ","))
            write(f2, "\r\n")
        end
        write(f2, "--end of contingency generator\r\n")
        write(f2, "--bus\r\n")
        write(f2, "contingency id,bus id,v(pu),theta(deg)\r\n")
        for i in 1:size(pr3)[1]
            write(f2, join(pr3[i,:], ","))
            write(f2, "\r\n")
        end
        write(f2, "--end of bus\r\n")
        write(f2, "--Delta\r\n")
        write(f2, "contingency id,Delta(MW)\r\n")
        for i in 1:size(pr4)[1]
            write(f2, join(pr4[i,:], ","))
            write(f2, "\r\n")
        end
        write(f2, "--end of Delta\r\n")
        write(f2, "--line flow\r\n")
        write(f2, "-contingency id,line id,origin bus id,destination bus id,circuit id,p_origin(MW),q_origin(MVar),p_destination(MW),q_destination(MVar)\r\n")
        for i in 1:size(pr5)[1]
            write(f2, join(pr5[i,:], ","))
            write(f2, "\r\n")
        end
        write(f2, "--end of line flow\r\n")
    end
end
