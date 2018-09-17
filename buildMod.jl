using CSV, Ipopt, JuMP;

function buildMod(BLGS, Br, contingency, genSeg)
    PSCOPF = Model(solver=IpoptSolver(print_level=1))
    ####################################################################################### ALL UP GOES TO dProc function
    # Constant term c in cost function
    c = BLGS[:,22];
    # Linear term b in cost function
    b = BLGS[:,23];
    # Quadratic term a in cost function
    a = BLGS[:,24];
    # Participation factor of generator
    part_factor = BLGS[:,25];

    #_______________________________#3. Solution variables#________________________
    # Number of buses
    NBus = size(BLGS,1);
    # Number of branches (power lines and transformers)
    NBr = size(Br,1);

    ### Connection matrices
    # list of "from" buses
    f = round.(Int,Br[:,1]);
    # list of "to" buses
    t = round.(Int,Br[:,2]);
    #connection matrix for line & from buses
    Cf = sparse(collect(1:NBr), f, ones(NBr), NBr, NBus);
    # connection matrix for line & to buses
    Ct = sparse(collect(1:NBr), t, ones(NBr), NBr, NBus);
    # series admittance
    Ys = 1 ./ (Br[:, 3] + im * Br[:, 4]);
    # line charging susceptance
    Bc = Br[:, 5];
    # transformer phase shirters
    shift = zeros(size(Br,1),1);
    # tap ratios
    tap = Br[:,7].*exp.(im*pi/180*shift);
    # 4 admittance co-matrices
    Ytt = Ys + im*Bc/2;
    Yff = Ytt ./ (tap .* conj(tap));
    Yft = - Ys ./ conj(tap);
    Ytf = - Ys ./ tap;
    # Shunt admittance
    Ysh = (BLGS[:, 20] + im * BLGS[:, 21]);
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
    Pmax = BLGS[:,18];
    Pmin = BLGS[:,19];
    Qmax = BLGS[:,14];
    Qmin = BLGS[:,15];

    # Branch limits (MVA)
    RATE_A = Br[:,6]

    # Voltage limits
    Vmax = BLGS[1,6]; # All values are the same
    Vmin = BLGS[1,7]; # All values are the same

    # Demand
    Pd = BLGS[:,10];
    Qd = BLGS[:,11];

    ####################### Incidence Matrix [I x L] #############################
    A = zeros(size(BLGS,1),size(Br,1));
    for i = 1:size(BLGS,1)
        for j = 1:size(Br,1)
            if Br[j,1] == BLGS[i,1]
                A[i,j] = 1
            elseif Br[j,2] == BLGS[i,1]
                A[i,j] = -1
            else A[i,j] = 0
            end
        end
    end

    ##################### Contingencies ##########################################
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
    ############################# SET SIZES ######################################
    # Number of contingencies
    NK = size(aL,1);

    ############################# SET LINES ######################################
    links = Array{Tuple{Int16, Int16}}(NBr);
    for l=1:NBr
      links[l] = (Br[l,1], Br[l,2])   # line from-to (directed)
    end
    ######################### MODEL DEFINITION ####################################
    #----------------------------- VARIABLES -------------------------------------
    # Nodal active power generation
    @variable(PSCOPF, p[i=1:NBus, k=1:NK])
    # Nodal reactive power generation
    @variable(PSCOPF, q[i=1:NBus, k=1:NK])
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
    @variable(PSCOPF, Delta);

    #--------------------------------- CONSTRAINTS -------------------------------
    # Inner expression for active and reactive power definitions (in nodal and branches)
    @NLexpression(PSCOPF, InnerP[i=1:NBus,j=1:NBus,k=1:NK], GL[i,j]*cos(δ[i,k]-δ[j,k]) + BL[i,j]*sin(δ[i,k]-δ[j,k]))
    @NLexpression(PSCOPF, InnerQ[i=1:NBus,j=1:NBus,k=1:NK], GL[i,j]*sin(δ[i,k]-δ[j,k]) - BL[i,j]*cos(δ[i,k]-δ[j,k]))
    # Nodal active power balance
    @NLconstraint(PSCOPF, EqConstrP[i=1:NBus,k=1:NK], p[i,k] - Pd[i] == V[i,k]*sum(V[j,k]*InnerP[i,j,k] for j=1:NBus));
    # Nodal reactive power balance
    @NLconstraint(PSCOPF, EqConstrQ[i=1:NBus,k=1:NK], q[i,k] - Qd[i] == V[i,k]*sum(V[j,k]*InnerQ[i,j,k] for j=1:NBus));

    # Active power flow
    @NLconstraint(PSCOPF,DefFlowP_ij[l=1:NBr,k=1:NK], pl[links[l][1],links[l][2],k] ==
    (V[links[l][1],k]*V[links[l][2],k]*InnerP[links[l][1],links[l][2],k] - V[links[l][1],k]^2*GL[links[l][1],links[l][2]])*aL[k,l])

    @NLconstraint(PSCOPF,DefFlowP_ji[l=1:NBr,k=1:NK], pl[links[l][2],links[l][1],k] ==
    (V[links[l][2],k]*V[links[l][1],k]*InnerP[links[l][2],links[l][1],k] - V[links[l][2],k]^2*GL[links[l][2],links[l][1]])*aL[k,l])

    # Reactive power flow
    @NLconstraint(PSCOPF,DefFlowQ_ij[l=1:NBr,k=1:NK], ql[links[l][1],links[l][2],k] ==
    (V[links[l][1],k]*V[links[l][2],k]*InnerQ[links[l][1],links[l][2],k] + V[links[l][1],k]^2*BL[links[l][1],links[l][2]])*aL[k,l]);

    @NLconstraint(PSCOPF,DefFlowQ_ji[l=1:NBr,k=1:NK], ql[links[l][2],links[l][1],k] ==
    (V[links[l][2],k]*V[links[l][1],k]*InnerQ[links[l][2],links[l][1],k] + V[links[l][2],k]^2*BL[links[l][2],links[l][1]])*aL[k,l]);

    # Line limits on Apparent Power
    @NLconstraint(PSCOPF, FlowLimits_ij[l=1:NBr,k=1:NK], (pl[links[l][1],links[l][2],k]^2 + ql[links[l][1],links[l][2],k]^2) <= RATE_A[l]^2);
    @NLconstraint(PSCOPF, FlowLimits_ji[l=1:NBr,k=1:NK], (pl[links[l][2],links[l][1],k]^2 + ql[links[l][2],links[l][1],k]^2) <= RATE_A[l]^2);

    # Generator active power output bound
    # slack bus in node 1
    @constraint(PSCOPF, GenLimitsP[i=1:NBus,k=1:NK], Pmin[i] <= p[i,k] <= Pmax[i]);

    # Generator reactive power output bound
    # slack bus in node 1
    @constraint(PSCOPF, GenLimitsQ[i=1:NBus,k=1:NK], Qmin[i] <= q[i,k] <= Qmax[i]);

    # Voltage magnitude limits
    @constraint(PSCOPF, VoltageLimits[i=1:NBus,k=1:NK], Vmin <= V[i,k] <= Vmax);

    # Angle limits
    @constraint(PSCOPF, slackBus[k=1:NK], δ[1,k] == 0);
    @constraint(PSCOPF, limitAngle[i=2:NBus,k=1:NK], -pi/2 <= δ[i,k] <= pi/2);

    ## Objective function
    @NLconstraint(PSCOPF, ObjectiveFunction, (Cost == sum(a[i]*(p[i])^2 + b[i]*p[i] + c[i]  for i=1:NBus)));
    @objective(PSCOPF, Min, Cost)

    ################################# AGC #########################################
    # Voltage control only on generator buses
    @constraint(PSCOPF, AGCLower[i in genSeg[:,1]], (q[i,2] - Qmin[i])*(V[i,2] - V[i,1]) <= 0);
    @constraint(PSCOPF, AGCUpper[i in genSeg[:,1]], (q[i,2] - Qmax[i])*(V[i,2] - V[i,1]) <= 0);
    total_part_factor = sum(part_factor);
    # Delta - powerfall constraint
    @constraint(PSCOPF, Delta == sum(p[i,2] for i in genSeg[:,1]) - sum(p[i,1] for i in genSeg[:,1]));

    @constraint(PSCOPF, GCLower[i in genSeg[:,1]], (p[i,2] - Pmin[i])*(p[i,2] - p[i,1] - (part_factor[i]/total_part_factor).*Delta) <= 0);
    @constraint(PSCOPF, GCUpper[i in genSeg[:,1]], (p[i,2] - Pmax[i])*(p[i,2] - p[i,1] - (part_factor[i]/total_part_factor).*Delta) <= 0);

    #@constraint(PSCOPF, GCLower[i in genSeg[:,1]], (p[i,2] - Pmin[i])*(p[i,2] - p[i,1] - part_factor[i] * (sum(p[i,2] for i in genSeg[:,1]) - sum(p[i,1] for i in genSeg[:,1]) )) <= 0);
    #@constraint(PSCOPF, GCUpper[i in genSeg[:,1]], (p[i,2] - Pmax[i])*(p[i,2] - p[i,1] - part_factor[i] * (sum(p[i,2] for i in genSeg[:,1]) - sum(p[i,1] for i in genSeg[:,1]) )) <= 0);

    # Print the model
    #print(PSCOPF)

    ## resolution
    return PSCOPF
end
