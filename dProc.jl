# preprocessing: read in the data and divide it into segments
using CSV;
function dProc(rawFile, genFile, contFile)

    rawData = readdlm(rawFile,',', skipblanks=false);
    n,m = size(rawData);
    busStartL = 4; # first three lines are headers
    busEndL = 0;
    loadStartL = 0;
    loadEndL = 0;
    shuntStartL = 0;
    shuntEndL = 0;
    genStartL = 0;
    genEndL = 0;
    branchStartL = 0;
    branchEndL = 0;
    transformerStartL = 0;
    transformerEndL = 0;
    for i in 1:n
      if contains(string(rawData[i,1]),"END OF BUS DATA")
        busEndL = i - 1;
        loadStartL = i + 1;
      end
      if contains(string(rawData[i,1]),"END OF LOAD DATA")
        loadEndL = i - 1;
        shuntStartL = i + 1;
      end
      if contains(string(rawData[i,1]),"END OF FIXED SHUNT DATA")
        shuntEndL = i - 1;
        genStartL = i + 1;
      end
      if contains(string(rawData[i,1]),"END OF GENERATOR DATA")
        genEndL = i - 1;
        branchStartL = i + 1;
      end
      if contains(string(rawData[i,1]),"END OF BRANCH DATA")
        branchEndL = i - 1;
        transformerStartL = i + 1;
      end
      if contains(string(rawData[i,1]),"END OF TRANSFORMER DATA")
        transformerEndL = i - 1;
      end
    end

    baseMVA = rawData[1,2];
    busSeg = rawData[busStartL:busEndL,:];
    loadSeg = rawData[loadStartL:loadEndL,:];
    shuntSeg = rawData[shuntStartL:shuntEndL,:];
    genSeg = rawData[genStartL:genEndL,:];
    branchSeg = rawData[branchStartL:branchEndL,:];
    transformerSeg = rawData[transformerStartL:transformerEndL,:];

    # Generator cost function
    gen_cost = CSV.read(genFile);
    gen_cost = convert(Array,gen_cost);

    # Contingency
    contingency = CSV.read(contFile);

    #---------------------------------- 2. DEFINING VARIABLES -----------------------------------#
    #################### Bus, load, Generator, Shunt table ########################
    # 1: Bus ID
    # 2: base KV
    # 3: bus type
    # 4: Vm, pu
    # 5: Va, rad
    # 6: normal Vmax, pu
    # 7: normal Vmin, pu
    # 8: emergency Vmax, pu
    # 9: emergency Vmin, pu
    # 10: Pd, pu
    # 11: Qd, pu
    # 12: genSeg, pu
    # 13: Qg, pu
    # 14: Qmax, pu
    # 15: Qmin, pu
    # 16: Vg
    # 17: baseMVA;
    # 18: Pmax, pu;
    # 19: Pmin, pu;
    # 20: Gs, pu
    # 21: Bs, pu
    # 22. Constant term c in cost function
    # 23. Linear term b in cost function
    # 24. Quadratic term a in cost function
    # 25. Participation factor of generator

    # Own created BLGS for data extraction with 25 columns
    BLGS = zeros(size(busSeg[:,1],1),25);
    # Bus ID
    BLGS[:,1] = busSeg[:,1];
    # base kVA
    BLGS[:,2] = busSeg[:,3];
    # bus Type
    BLGS[:,3] = zeros(size(busSeg[:,1],1),1);
    # Initial voltage magnitude
    BLGS[:,4] = ones(size(busSeg[:,1],1),1);
    # Initial voltage Angle
    BLGS[:,5] = zeros(size(busSeg[:,1],1),1);
    # normal Vmax
    BLGS[:,6] = busSeg[:,10];
    # normal Vmin
    BLGS[:,7] = busSeg[:,11];
    # emergency Vmax
    BLGS[:,8] = busSeg[:,10];
    # emergency Vmin
    BLGS[:,9] = busSeg[:,11];
    # P demand
    Pd = hcat(loadSeg[:,1],loadSeg[:,6]/baseMVA);
    for i = 1:size(Pd,1)
        for j = 1:size(BLGS,1)
            if BLGS[j,1] == Pd[i,1]
                BLGS[j,10] = Pd[i,2]
            end
        end
    end

    # Q demand
    Qd = hcat(loadSeg[:,1],loadSeg[:,7]/baseMVA);
    for i = 1:size(Qd,1)
        for j = 1:size(BLGS,1)
            if BLGS[j,1] == Qd[i,1]
                BLGS[j,11] = Qd[i,2]
            end
        end
    end
    # P generation
    BLGS[:,12] = ones(size(busSeg[:,1],1),1);
    # Q generation
    BLGS[:,13] = ones(size(busSeg[:,1],1),1);
    # Q max
    Qmax = hcat(genSeg[:,1],genSeg[:,5]/baseMVA);
    for i = 1:size(Qmax,1)
        for j = 1:size(BLGS,1)
            if BLGS[j,1] == Qmax[i,1]
                BLGS[j,14] = Qmax[i,2]
            end
        end
    end
    # Q min
    Qmin = hcat(genSeg[:,1],genSeg[:,6]/baseMVA);
    for i = 1:size(Qmin,1)
        for j = 1:size(BLGS,1)
            if BLGS[j,1] == Qmin[i,1]
                BLGS[j,15] = Qmin[i,2]
            end
        end
    end
    # Voltage magnitude on generator
    BLGS[:,16] = ones(size(busSeg[:,1],1),1);
    # P max
    Pmax = hcat(genSeg[:,1],genSeg[:,17]/baseMVA);
    for i = 1:size(Pmax,1)
        for j = 1:size(BLGS,1)
            if BLGS[j,1] == Pmax[i,1]
                BLGS[j,18] = Pmax[i,2]
            end
        end
    end
    # P Min
    Pmin = hcat(genSeg[:,1],genSeg[:,18]/baseMVA);
    for i = 1:size(Pmin,1)
        for j = 1:size(BLGS,1)
            if BLGS[j,1] == Pmin[i,1]
                BLGS[j,19] = Pmin[i,2]
            end
        end
    end
    # Conductance Gs
    Gs = hcat(shuntSeg[:,1],shuntSeg[:,4]/baseMVA);
    for i = 1:size(Gs,1)
        for j = 1:size(BLGS,1)
            if BLGS[j,1] == Gs[i,1]
                BLGS[j,20] = Gs[i,2]
            end
        end
    end
    # Susceptance Bs
    Bs = hcat(shuntSeg[:,1],shuntSeg[:,5]/baseMVA);
    for i = 1:size(Bs,1)
        for j = 1:size(BLGS,1)
            if BLGS[j,1] == Bs[i,1]
                BLGS[j,21] = Bs[i,2]
            end
        end
    end

    ## Generation cost and participation factor
    for i in unique(gen_cost[:,1])
        index = findfirst(gen_cost[:,1],i)
        BLGS[i,22:25] = gen_cost[index:index+3,4]
    end

    ################### BRANCHES TABLE (power lines and transformers )#############
    # 1 : FROM
    # 2 : TO
    # 3 : r
    # 4 : x
    # 5 : b
    # 6 : RATE A
    # 7 : Kt - transformation ratio

    # size of branch matrix
    # Number of lines
    NL = size(branchSeg,1);
    # NUmber of transformers
    NT = Int(size(transformerSeg,1)/4);
    # size of branch matrix
    Br = zeros(NL + NT,7);
    # From buses
    Br[1:NL,1] = branchSeg[:,1];
    Br[NL + 1 : NL + NT,1] = transformerSeg[range(1,4,NT),1];
    # To buses
    Br[1:size(branchSeg,1),2] = branchSeg[:,2];
    Br[size(branchSeg,1)+1:size(branchSeg,1)+3,2] = transformerSeg[[1,5,9],2];
    # r
    Br[1:NL,3] = branchSeg[:,4];
    Br[NL + 1 : NL + NT,3] = transformerSeg[range(2,4,NT),1];
    # x
    Br[1:NL,4] = branchSeg[:,5];
    Br[NL + 1 : NL + NT,4] = transformerSeg[range(2,4,NT),2];
    # b
    Br[1:NL,5] = branchSeg[:,6];
    Br[NL + 1 : NL + NT,5] = repeat([0.0],outer=[NT]);
    # RATE A
    Br[1:NL,6] = branchSeg[:,7]/baseMVA;
    Br[NL+1:NL+NT,6] = transformerSeg[range(2,4,NT),3]/baseMVA;
    # Kt
    Br[:,7] = ones(NL + NT);
    Br[NL + 1 : NL + NT,7] = transformerSeg[range(3,4,NT),1];
    return BLGS, Br, contingency, genSeg, baseMVA
end
