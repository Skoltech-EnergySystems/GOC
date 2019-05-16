using CSV, Ipopt, JuMP;

function buildMod(fData,uData, contDList, contingency)
    PSCOPF = Model(solver=IpoptSolver(print_level=0))

    baseMVA = fData.baseMVA;
    bList = fData.busList;
    bData = fData.busDList;
    gList = fData.genList;
    gData = fData.genDList;
    brList = fData.brList;
    brData = fData.brDList;
    S = uData.contList;
    contData = uData.contDList;
    # Constant term c in cost function
    # Constant term c in cost function
    Dictcost = collect(gData[k].cParams for k in sort(collect(keys(gData))));
    c = collect(Dictcost[k][:0] for k in sort(collect(keys(Dictcost))));
    # Linear term b in cost function
    b = collect(Dictcost[k][:1] for k in sort(collect(keys(Dictcost))));
    # Quadratic term a in cost function
    a = collect(Dictcost[k][:2] for k in sort(collect(keys(Dictcost))));
    # Participation factor of generator
    part_factor = collect(gData[k].alpha for k in sort(collect(keys(gData))));

    #_______________________________#3. Solution variables#________________________
    # Number of buses
    NBus = size(bList)[1];
    # Number of branches (power lines and transformers)
    NBr = size(brList)[1];
    # Number of generators
    NGen = size(gList)[1];
    ### Connection matrices
    # list of "from" buses

    br_int = round.(Int,zeros(size(brList,1),2))
    for i = 1:size(bList,1)
        for j = 1:size(brList,1)
            if collect(bData[k].ID_ext for k in sort(collect(keys(bData))))[i] == collect(brData[k].ID_ext for k in sort(collect(keys(brData))))[j][1]
                br_int[j,1] = collect(bData[k].ID_int for k in sort(collect(keys(bData))))[i]
            end
            if collect(bData[k].ID_ext for k in sort(collect(keys(bData))))[i] == collect(brData[k].ID_ext for k in sort(collect(keys(brData))))[j][2]
                br_int[j,2] = collect(bData[k].ID_int for k in sort(collect(keys(bData))))[i]
            end
        end
    end
    #br_int
    f = br_int[:,1];

    # list of "to" buses
    t = br_int[:,2];
    #connection matrix for line & from buses
    Cf = sparse(collect(1:NBr), f, ones(NBr), NBr, NBus);
    # connection matrix for line & to buses
    Ct = sparse(collect(1:NBr), t, ones(NBr), NBr, NBus);
    # series admittance
    Ys = 1 ./ (collect(brData[k].r for k in sort(collect(keys(brData)))) + im * collect(brData[k].x for k in sort(collect(keys(brData)))))
    # line charging susceptance
    Bc = collect(brData[k].bc for k in sort(collect(keys(brData))))
    # transformer phase shirters
    shift = zeros(size(brList,1),1)
    # tap ratios
    tap = collect(brData[k].tau for k in sort(collect(keys(brData)))).*exp.(im*pi/180*shift)

    # admittance co-matrices
    Ytt = Ys + im*Bc/2;
    Yff = Ytt ./ (tap .* conj(tap));
    Yft = - Ys ./ conj(tap);
    Ytf = - Ys ./ tap;
    # Shunt admittance
    Ysh = (collect(bData[k].gsh for k in sort(collect(keys(bData)))) + im * collect(bData[k].bsh for k in sort(collect(keys(bData)))))
    # From and To admittance co-matrices
    i = [1:NBr; 1:NBr]; # In Matpower this line is transposed
    Yf = sparse(i, [f; t], [Yff; Yft][:,1], NBr, NBus);
    Yt = sparse(i, [f; t], [Ytf; Ytt][:,1], NBr, NBus);
    # Ybus
    Ybus = transpose(Cf) * Yf + transpose(Ct) * Yt + sparse(1:NBus, 1:NBus, Ysh, NBus, NBus);
    # Conductance and susceptance of Ybus
    GL = real.(Ybus);
    BL = imag.(Ybus);
    # Generator limits
    Pmax = collect(gData[k].Pmax for k in sort(collect(keys(gData))))
    Pmin = collect(gData[k].Pmin for k in sort(collect(keys(gData))))
    Qmax = collect(gData[k].Qmax for k in sort(collect(keys(gData))))
    Qmin = collect(gData[k].Qmin for k in sort(collect(keys(gData))))
    # Branch limits (MVA)
    RATE_A = collect(brData[k].t for k in sort(collect(keys(brData))))
    # Voltage limits
    Vmax = collect(bData[k].Vmax for k in sort(collect(keys(bData))))
    Vmin = collect(bData[k].Vmin for k in sort(collect(keys(bData))))
    # Demand
    Pd = collect(bData[k].Pd for k in sort(collect(keys(bData))))
    Qd = collect(bData[k].Qd for k in sort(collect(keys(bData))))

    ##################### Contingencies ##########################################
    # All lines and transformers are online
    aL = transpose(ones(NBr,1));
    # Contingency FROM buses
    Cont_FR = float(contingency[:,3])[1,1]
    # Contingency TO buses
    Cont_TO = float(contingency[:,4])[1,1]
    # Number of contingencies
    Cont_size = size(contingency,1);
    # Contingency matrix
    Br_ind = hcat(collect(1:NBr),collect(brData[k].From for k in sort(collect(keys(brData)))),collect(brData[k].To for k in sort(collect(keys(brData)))));
    aL_loop = transpose(ones(NBr,1));
    for i = 1:NBr
        if (Br_ind[i,2]==Cont_FR && Br_ind[i,3]==Cont_TO)
            aL_loop[1,i]=0;
            aL=vcat(aL,aL_loop)
        end
    end

    ############################# SET SIZES ######################################
    # Number of contingencies
    NK = size(aL,1);

    ############################# SET LINES ######################################
    links = Array{Tuple{Int16, Int16}}(NBr);
    for k=1:NBr
      links[k] = (collect(brData[k].From for k in sort(collect(keys(brData))))[k,1], collect(brData[k].To for k in sort(collect(keys(brData))))[k,1])   # line from-to (directed)
    end

    println("Building the model")
