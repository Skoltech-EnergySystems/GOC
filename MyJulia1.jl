using GLPKMathProgInterface;
#Grid Optimization Competition
using DataFrames, DataArrays, CSV, Ipopt, JuMP, ConditionalJuMP;

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
Br = zeros(size(branchSeg,1)+3,7);
# From buses
Br[1:size(branchSeg,1),1] = branchSeg[:,1];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,1] = transformerSeg[[1,5,9],1];
# To buses
Br[1:size(branchSeg,1),2] = branchSeg[:,2];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,2] = transformerSeg[[1,5,9],2];
# r
Br[1:size(branchSeg,1),3] = branchSeg[:,4];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,3] = transformerSeg[[2,6,10],1];
# x
Br[1:size(branchSeg,1),4] = branchSeg[:,5];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,4] = transformerSeg[[2,6,10],2];
# b
Br[1:size(branchSeg,1),5] = branchSeg[:,6];
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,5] = [0.0,0.0,0.0];
# RATE A
Br[1:size(branchSeg,1),6] = branchSeg[:,7]/baseMVA;
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,6] = transformerSeg[[2,6,10],3]/baseMVA;
# Kt
Br[:,7] = ones(size(branchSeg,1)+3);
Br[size(branchSeg,1)+1:size(branchSeg,1)+3,7] = (transformerSeg[[3,7,11],1]);

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
RATE_A = Br[:,6];

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

############################# SET LINES ######################################
links = Array{Tuple{Int16, Int16}}(NBr);
for l=1:NBr
  links[l] = (Br[l,1], Br[l,2])   # line from-to (directed)
end
######################### MODEL DEFINITION ####################################
#----------------------------- VARIABLES -------------------------------------
# Nodal active power generation
@variable(PSCOPF, p[i=1:NBus]);
# Nodal reactive power generation;
@variable(PSCOPF, q[i=1:NBus]);
# Branch active power flows
@variable(PSCOPF, pl[i=1:NBus, j=1:NBus]);
# Branch reactive power flows
@variable(PSCOPF, ql[i=1:NBus, j=1:NBus]);
# Voltage magnitude
@variable(PSCOPF, V[i=1:NBus], start = 1.0);
#,Voltage phase angle
@variable(PSCOPF, δ[i=1:NBus], start = 0.0);
# Active power generation cost
@variable(PSCOPF, Cost, lowerbound=0);

#--------------------------------- CONSTRAINTS -------------------------------
# Inner expression for active and reactive power definitions (in nodal and branches)
@NLexpression(PSCOPF, InnerP[i=1:NBus,j=1:NBus], GL[i,j]*cos(δ[i]-δ[j]) + BL[i,j]*sin(δ[i]-δ[j]));
@NLexpression(PSCOPF, InnerQ[i=1:NBus,j=1:NBus], GL[i,j]*sin(δ[i]-δ[j]) - BL[i,j]*cos(δ[i]-δ[j]));
# Nodal active power balance
@NLconstraint(PSCOPF, EqConstrP[i=1:NBus], p[i] - Pd[i] == V[i]*sum(V[j]*InnerP[i,j] for j=1:NBus));
# Nodal reactive power balance
@NLconstraint(PSCOPF, EqConstrQ[i=1:NBus], q[i] - Qd[i] == V[i]*sum(V[j]*InnerQ[i,j] for j=1:NBus));

# Active power flow
@NLconstraint(PSCOPF,DefFlowP_ij[l=1:NBr], pl[links[l][1],links[l][2]] ==
(V[links[l][1]]*V[links[l][2]]*InnerP[links[l][1],links[l][2]] - V[links[l][1]]^2*GL[links[l][1],links[l][2]]));

@NLconstraint(PSCOPF,DefFlowP_ji[l=1:NBr], pl[links[l][2],links[l][1]] ==
(V[links[l][2]]*V[links[l][1]]*InnerP[links[l][2],links[l][1]] - V[links[l][2]]^2*GL[links[l][2],links[l][1]]));

# Reactive power flow
@NLconstraint(PSCOPF,DefFlowQ_ij[l=1:NBr], ql[links[l][1],links[l][2]] ==
(V[links[l][1]]*V[links[l][2]]*InnerQ[links[l][1],links[l][2]] + V[links[l][1]]^2*BL[links[l][1],links[l][2]]));

@NLconstraint(PSCOPF,DefFlowQ_ji[l=1:NBr], ql[links[l][2],links[l][1]] ==
(V[links[l][2]]*V[links[l][1]]*InnerQ[links[l][2],links[l][1]] + V[links[l][2]]^2*BL[links[l][2],links[l][1]]));

# Line limits on Apparent Power
@NLconstraint(PSCOPF, FlowLimits_ij[l=1:NBr], (pl[links[l][1],links[l][2]]^2 + ql[links[l][1],links[l][2]]^2) <= RATE_A[l]^2);
@NLconstraint(PSCOPF, FlowLimits_ji[l=1:NBr], (pl[links[l][2],links[l][1]]^2 + ql[links[l][2],links[l][1]]^2) <= RATE_A[l]^2);

# Generator active power output bound
# slack bus in node 1
@constraint(PSCOPF, GenLimitsP[i=1:NBus], Pmin[i] <= p[i] <= Pmax[i]);

# Generator reactive power output bound
# slack bus in node 1
@constraint(PSCOPF, GenLimitsQ[i=1:NBus], Qmin[i] <= q[i] <= Qmax[i]);

# Voltage magnitude limits
@constraint(PSCOPF, VoltageLimits[i=1:NBus], Vmin <= V[i] <= Vmax);

# Angle limits
@constraint(PSCOPF, δ[1] == 0);
@constraint(PSCOPF, limitAngle[i=2:NBus], -pi/2 <= δ[i] <= pi/2);

## Objective function
@NLconstraint(PSCOPF, ObjectiveFunction, (Cost == sum(a[i]*(p[i])^2 + b[i]*p[i] + c[i]  for i=1:NBus)));
@objective(PSCOPF, Min, Cost)

# Print the model
#print(PSCOPF)

## resolution
solve(PSCOPF)
############################ OUTPUT FILES ##################################
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
    writedlm(f1, "--generation dispatch\r\n")
    writedlm(f1, "bus id,unit id,genSeg(MW),qg(MVar)\r\n")
    for i in 1:sizePg
        writedlm(f1, join(pr[i,:], ","))
        writedlm(f1, "\r\n")
    end
    writedlm(f1, "--end of generation dispatch")
end
end
