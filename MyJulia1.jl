#using GLPKMathProgInterface; # Chnage to Gurobi
#Grid Optimization Competition
using CSV, Ipopt, JuMP;

function MyJulia1(rawFile, genFile, contFile)
# Preventive Security Constrained Optimal Power Flow
PSCOPF = Model(solver=IpoptSolver(print_level=0))
#PSCOPF = Model(solver=IpoptSolver())
#---------------------------------- 1. IMPORTING DATA -----------------------------------#
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
# 17: baseMVA
# 18: Pmax, pu
# 19: Pmin, pu
# 20: Gs, pu
# 21: Bs, pu
# 22. Constant term c in cost function
# 23. Linear term b in cost function
# 24. Quadratic term a in cost function
# 25. Participation factor of generator

BLGS = zeros(14,25);
# Bus ID
BLGS[:,1] = busSeg[:,1];
# base kVA
BLGS[:,2] = busSeg[:,3];
# bus Type
BLGS[:,3] = zeros(14,1);
# Initial voltage magnitude
BLGS[:,4] = ones(14,1);
# Initial voltage Angle
BLGS[:,5] = zeros(14,1);
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
BLGS[:,12] = ones(14,1);
# Q generation
BLGS[:,13] = ones(14,1);
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
BLGS[:,16] = ones(14,1);
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

# Constant term c in cost function
c = BLGS[:,22];
# Linear term b in cost function
b = BLGS[:,23];
# Quadratic term a in cost function
a = BLGS[:,24];
# Participation factor of generator
part_factor = BLGS[:,25];
#convert(DataFrame,BLGS)
#rename(BLGS, :x1 => :Bus_ID)
################### BRANCHES TABLE (power lines and transformers )#############
# 1 : FROM
# 2 : TO
# 3 : r
# 4 : x
# 5 : b
# 6 : RATE A
# 7 : Kt
# size of branch matrix
Br = zeros(size(branchSeg,1)+3,7)
# From buses
Br[1:size(branchSeg,1),1] = branchSeg[:,1];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,1] = transformerSeg[[1,5,9],1];
# To buses
Br[1:size(branchSeg,1),2] = branchSeg[:,2];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,2] = transformerSeg[[1,5,9],2];
# r
Br[1:size(branchSeg,1),3] = branchSeg[:,4];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,3] = transformerSeg[[2,6,10],1]
# x
Br[1:size(branchSeg,1),4] = branchSeg[:,5];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,4] = transformerSeg[[2,6,10],2]
# b
Br[1:size(branchSeg,1),5] = branchSeg[:,6];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,5] = [0.0,0.0,0.0]
# RATE A
Br[1:size(branchSeg,1),6] = branchSeg[:,7]/baseMVA
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,6] = transformerSeg[[2,6,10],3]/baseMVA
# Kt
Br[:,7] = ones(size(branchSeg,1)+3)
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,7] = (transformerSeg[[3,7,11],1])

#_______________________________#3. Solution variables#________________________
# Number of buses
NBus = size(BLGS,1);
# Number of branches (power lines and transformers)
NBr = size(Br,1);
# Reactance, pu
X = Br[:,4]
# Resistance, pu
R = Br[:,3];
# Impedance, pu
Z = (X.^2+R.^2).^0.5;

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
# Post-contingency real power shortfall
@expression(PSCOPF, Delta, sum(p[i,2] for i in genSeg[:,1]) - sum(p[i,1] for i in genSeg[:,1]))

#@constraint(PSCOPF, GCLower[i in genSeg[:,1]], (p[i,2] - Pmin[i])*(p[i,2] - p[i,1] - (part_factor[i]/total_part_factor).*Delta) <= 0);
#@constraint(PSCOPF, GCUpper[i in genSeg[:,1]], (p[i,2] - Pmax[i])*(p[i,2] - p[i,1] - (part_factor[i]/total_part_factor).*Delta) <= 0);

@constraint(PSCOPF, GCLower[i in genSeg[:,1]], (p[i,2] - Pmin[i])*(p[i,2] - p[i,1] - part_factor[i]*Delta) <= 0);
@constraint(PSCOPF, GCUpper[i in genSeg[:,1]], (p[i,2] - Pmax[i])*(p[i,2] - p[i,1] - part_factor[i]*Delta) <= 0);

# Print the model
#print(PSCOPF)

## resolution
solve(PSCOPF)
############################ OUTPUT FILES ##################################
#cd("C:/Users/Ильгиз/Documents/Документы_Ильгиз/Skoltech_2018/Grid Competition/Code")
# Extracting calculated generation values for base case
Pgen = [];
for y in getvalue(p)[:,1]
    if abs(y) > 0.0001
        push!(Pgen,y*100)
    end
end

Qgen = [];
for z in getvalue(q)[:,1]
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

# Writing into the file "solution1.txt"
open("solution1.txt", "w") do f1
    write(f1, "--generation dispatch\r\n")
    write(f1, "bus id,unit id,genSeg(MW),qg(MVar)\r\n")
    for i in 1:sizePg
        write(f1, join(pr[i,:], ","))
        write(f1, "\r\n")
    end
    write(f1, "--end of generation dispatch")
end

# Extracting calculated generation values for contingency case
Pgenk = [];
for y in getvalue(p)[:,2]
    if abs(y) > 0.0001
        push!(Pgenk,y*100)
    end
end
Pgenk

Qgenk = [];
for z in getvalue(q)[:,2]
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

## 2. contingency bus information
pr3 = Array{Any}(2*NBus,4);
# contingency ID
pr3[1:NBus,1] = repeat(["0"],outer=[NBus]); # 0. base case
pr3[NBus+1:2*NBus,1] = repeat(["1"],outer=[NBus]); # 1. contingency case
# bus ID
pr3[1:NBus,2] = round.(Int,BLGS[:,1]); # 0. base case
pr3[NBus+1:2*NBus,2] = round.(Int,BLGS[:,1]); # 1. contingency case
# Voltage in per unit
pr3[1:NBus,3] = getvalue(V)[:,1]; # 0. base case
pr3[NBus+1:2*NBus,3] = getvalue(V)[:,2]; # 1. contingency case
# Voltage angle in degree
pr3[1:NBus,4] = rad2deg.(getvalue(δ))[:,1]; # 0. base case
pr3[NBus+1:2*NBus,4] = rad2deg.(getvalue(δ))[:,2]; # 1. contingency case

## 3. contingency delta
pr4 = Array{Any}(NK-1,2);
pr4[1,1] = '1';
pr4[1,2] = getvalue(Delta)*100;

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
    push!(originP,getvalue(pl)[x,y,1]*baseMVA) # 0. base case
end
pr5[1:NBr,6] = originP;

originP = [];
for (x,y) in zip(round.(Int,Br[:,1]), round.(Int,Br[:,2]))
    push!(originP,getvalue(pl)[x,y,2]*baseMVA) # 1. contingency case
end
pr5[NBr+1:2*NBr,6] = originP;
# reactive power in MVar at origin
originQ = [];
for (x,y) in zip(round.(Int,Br[:,1]), round.(Int,Br[:,2]))
    push!(originQ,getvalue(ql)[x,y,1]*baseMVA) # 0. base case
end
pr5[1:NBr,7] = originQ;

originQ = [];
for (x,y) in zip(round.(Int,Br[:,1]), round.(Int,Br[:,2]))
    push!(originQ,getvalue(ql)[x,y,2]*baseMVA) # 1. contingency case
end
pr5[NBr+1:2*NBr,7] = originQ;

# real power in megawatts at destination
destP = [];
for (x,y) in zip(round.(Int,Br[:,2]), round.(Int,Br[:,1]))
    push!(destP,getvalue(pl)[x,y,1]*baseMVA) # 0. base case
end
pr5[1:NBr,8] = destP;

destP = [];
for (x,y) in zip(round.(Int,Br[:,2]), round.(Int,Br[:,1]))
    push!(destP,getvalue(pl)[x,y,2]*baseMVA) # 1. contingency case
end
pr5[NBr+1:2*NBr,8] = destP;

# reactive power in MVar at destination
destQ = [];
for (x,y) in zip(round.(Int,Br[:,2]), round.(Int,Br[:,1]))
    push!(destQ,getvalue(ql)[x,y,1]*baseMVA) # 0. base case
end
pr5[1:NBr,9] = destQ;

destQ = [];
for (x,y) in zip(round.(Int,Br[:,2]), round.(Int,Br[:,1]))
    push!(destQ,getvalue(ql)[x,y,2]*baseMVA) # 1. contingency case
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
