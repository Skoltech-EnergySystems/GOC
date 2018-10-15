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

    ######################### MODEL DEFINITION ####################################
    #----------------------------- VARIABLES -------------------------------------
    # Nodal active power generation
    @variable(PSCOPF, p[g=1:NGen, k=1:NK])
    # Nodal reactive power generation
    @variable(PSCOPF, q[g=1:NGen, k=1:NK])
    # Branch active power flows
    @variable(PSCOPF, pl[i=1:NBus, j=1:NBus, k=1:NK])
    # Branch reactive power flows
    @variable(PSCOPF, ql[i=1:NBus, j=1:NBus, k=1:NK])
    # Voltage magnitude
    @variable(PSCOPF, V[i=1:NBus, k=1:NK], start = 1.0)
    #,Voltage phase angle
    @variable(PSCOPF, δ[i=1:NBus, k=1:NK], start = 0.0)
    # Active power generation cost
    @variable(PSCOPF, Cost, lowerbound=0);
    # Delta - powerfall
    @variable(PSCOPF, Delta)

    # "generators-buses" incidence matrix
    Ω = zeros(NBus,NGen);

    for g = 1:NGen
        for i = 1:NBus
            if collect(gData[k].ID_ext for k in sort(collect(keys(gData))))[g] == collect(bData[k].ID_ext for k in sort(collect(keys(bData))))[i] &&
                 g in collect(gData[k].ID_int for k in sort(collect(keys(gData)))) &&
                 i in collect(bData[k].ID_int for k in sort(collect(keys(bData))))
                Ω[i,g] = 1
            end
        end
    end

    # Check of correctness: in total there are 99 generators
    #sum(Ω)
    #--------------------------------- CONSTRAINTS -------------------------------
    # Inner expression for active and reactive power definitions (in nodal and branches)
    @NLexpression(PSCOPF, InnerP[i=1:NBus,j=1:NBus,k=1:NK], GL[i,j]*cos(δ[i,k]-δ[j,k]) + BL[i,j]*sin(δ[i,k]-δ[j,k]))
    @NLexpression(PSCOPF, InnerQ[i=1:NBus,j=1:NBus,k=1:NK], GL[i,j]*sin(δ[i,k]-δ[j,k]) - BL[i,j]*cos(δ[i,k]-δ[j,k]))
    # Nodal active power balance
    @NLconstraint(PSCOPF, EqConstrP[i=1:NBus,k=1:NK], sum(Ω[i,g]*p[g,k] for g=1:NGen) - Pd[i] == V[i,k]*sum(V[j,k]*InnerP[i,j,k] for j=1:NBus));
    # Nodal reactive power balance
    @NLconstraint(PSCOPF, EqConstrQ[i=1:NBus,k=1:NK], sum(Ω[i,g]*q[g,k] for g=1:NGen) - Qd[i] == V[i,k]*sum(V[j,k]*InnerQ[i,j,k] for j=1:NBus));

    # Active power flow
    # br_int has the same functions as "links", but is done for consequtive bus numbers
    @NLconstraint(PSCOPF,DefFlowP_ij[l=1:NBr,k=1:NK], pl[br_int[l,1],br_int[l,2],k] ==
    (V[br_int[l,1],k]*V[br_int[l,2],k]*InnerP[br_int[l,1],br_int[l,2],k] - V[br_int[l,1],k]^2*GL[br_int[l,1],br_int[l,2]])*aL[k,l])

    @NLconstraint(PSCOPF,DefFlowP_ji[l=1:NBr,k=1:NK], pl[br_int[l,2],br_int[l,1],k] ==
    (V[br_int[l,2],k]*V[br_int[l,1],k]*InnerP[br_int[l,2],br_int[l,1],k] - V[br_int[l,2],k]^2*GL[br_int[l,2],br_int[l,1]])*aL[k,l])

    # Reactive power flow
    @NLconstraint(PSCOPF,DefFlowQ_ij[l=1:NBr,k=1:NK], ql[br_int[l,1],br_int[l,2],k] ==
    (V[br_int[l,1],k]*V[br_int[l,2],k]*InnerQ[br_int[l,1],br_int[l,2],k] + V[br_int[l,1],k]^2*BL[br_int[l,1],br_int[l,2]])*aL[k,l]);

    @NLconstraint(PSCOPF,DefFlowQ_ji[l=1:NBr,k=1:NK], ql[br_int[l,2],br_int[l,1],k] ==
    (V[br_int[l,2],k]*V[br_int[l,1],k]*InnerQ[br_int[l,2],br_int[l,1],k] + V[br_int[l,2],k]^2*BL[br_int[l,2],br_int[l,1]])*aL[k,l]);

    # Line limits on Apparent Power
    @NLconstraint(PSCOPF, FlowLimits_ij[l=1:NBr,k=1:NK], (pl[br_int[l,1],br_int[l,2],k]^2 + ql[br_int[l,1],br_int[l,2],k]^2) <= RATE_A[l]^2);
    @NLconstraint(PSCOPF, FlowLimits_ji[l=1:NBr,k=1:NK], (pl[br_int[l,2],br_int[l,1],k]^2 + ql[br_int[l,2],br_int[l,1],k]^2) <= RATE_A[l]^2);

    # Generator active power output bound
    # slack bus in node 1
    @constraint(PSCOPF, GenLimitsP[g=1:NGen,k=1:NK], Pmin[g] <= p[g,k] <= Pmax[g]);

    # Generator reactive power output bound
    # slack bus in node 1
    @constraint(PSCOPF, GenLimitsQ[g=1:NGen,k=1:NK], Qmin[g] <= q[g,k] <= Qmax[g]);

    # Voltage magnitude limits
    @constraint(PSCOPF, VoltageLimits[i=1:NBus,k=1:NK], Vmin[i] <= V[i,k] <= Vmax[i]);

    # Angle limits
    @constraint(PSCOPF, slackBus[k=1:NK], δ[1,k] == 0);
    @constraint(PSCOPF, limitAngle[i=2:NBus,k=1:NK], -pi/2 <= δ[i,k] <= pi/2);

    # Delta - powerfall constraint
    @constraint(PSCOPF, Delta == sum(p[g,2] for g = 1:NGen) - sum(p[g,1] for g = 1:NGen));

    # Post-contingency real power shortfall
    total_part_factor = sum(part_factor);
    @constraint(PSCOPF, PostCont[g=1:NGen], p[g,2] == p[g,1] + (part_factor[g]/total_part_factor).*Delta);
    #@constraint(PSCOPF, PostCont[g=1:NGen], p[g,2] == p[g,1] + part_factor[g]*Delta);

    ## Objective function
    @NLconstraint(PSCOPF, ObjectiveFunction, (Cost == sum(a[g]*(p[g,1])^2 + b[g]*p[g,1] + c[g]  for g=1:NGen)));
    @objective(PSCOPF, Min, Cost)

    ################################# AGC #########################################
    # Voltage control only on generator buses
    @constraint(PSCOPF, AGCLower[g=1:NGen, i=1:NBus], (q[g,2] - Qmin[g])*(V[i,2] - V[i,1]) <= 0);
    @constraint(PSCOPF, AGCUpper[g=1:NGen, i=1:NBus], (q[g,2] - Qmax[g])*(V[i,2] - V[i,1]) <= 0);

    ## resolution
    return (PSCOPF,NGen,NBus,NBr,NK,aL)
end
