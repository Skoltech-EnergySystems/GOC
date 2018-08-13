using GLPKMathProgInterface;
#Grid Optimization Competition
using DataFrames, DataArrays, CSV, Ipopt, JuMP;
# Preventive Security Constrained Optimal Power Flow
PSCOPF = Model(solver=IpoptSolver())

#---------------------------------- 1. IMPORTING DATA -----------------------------------#
cd("C:/Users/Ильгиз/Documents/Документы_Ильгиз/Skoltech_2018/Grid Competition/Phase_0_IEEE14_1Scenario/scenario_1")
all_data=readdlm("pscopf.m",',');
#Assume that in all files all_busdata will be written on the third row
# Number of buses in the system
bus_number=parse(Int,split(split(all_data[3,:][1],"a")[end]," ")[1]);

# All starting indeces
start_idx_busdata=findin(all_data,["all_busdata$bus_number = [ "]);
start_idx_loaddata=findin(all_data,["all_loaddata$bus_number = [ "]);
start_idx_shuntdata=findin(all_data,["all_fixedshuntdata$bus_number = [ "]);
start_idx_gendata=findin(all_data,["all_generatordata$bus_number = [ "]);
start_idx_linedata=findin(all_data,["all_linedata$bus_number = [ "]);
start_idx_trandata=findin(all_data,["all_transformerdata$bus_number = [ "]);

# All ending indeces
ends_idx_all=findin(all_data,["]; "]);
end_idx_busdata=ends_idx_all[1];
end_idx_loaddata=ends_idx_all[2];
end_idx_shuntdata=ends_idx_all[3];
end_idx_gendata=ends_idx_all[4];
end_idx_linedata=ends_idx_all[5];
end_idx_trandata=ends_idx_all[6];

# In bus table 12 non-empty columns
all_busdata=all_data[start_idx_busdata[1]+1:end_idx_busdata-1,1:12];
# In load table 13 non-empty columns
all_loaddata=all_data[start_idx_loaddata[1]+1:end_idx_loaddata-1,1:13];
# In shunt table 4 non-empty columns
all_shuntdata=all_data[start_idx_shuntdata[1]+1:end_idx_shuntdata-1,1:4]
# In generator table 19 non-empty columns
all_gendata=all_data[start_idx_gendata[1]+1:end_idx_gendata-1,1:19];
# In line table 17 non-empty columns
all_linedata=all_data[start_idx_linedata[1]+1:end_idx_linedata-1,1:17];
# In transformer table 40 non-empty columns
all_trandata=all_data[start_idx_trandata[1]+1:end_idx_trandata-1,1:40];

# Generator cost function
gen_cost=CSV.read("generator.csv");
gen_cost=convert(Array,gen_cost);

# Contingency
contingency=CSV.read("contingency.csv");

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
# 12: Pg, pu
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

BLGS=zeros(14,25);
# Bus ID
BLGS[:,1]=all_busdata[:,1];
# base kVA
BLGS[:,2]=all_busdata[:,2];
# bus Type
BLGS[:,3]=all_busdata[:,3];
# Initial voltage magnitude
BLGS[:,4]=all_busdata[:,7];
# Initial voltage Angle
BLGS[:,5]=deg2rad.(all_busdata[:,8]);
# normal Vmax
BLGS[:,6]=all_busdata[:,9];
# normal Vmin
BLGS[:,7]=all_busdata[:,10];
# emergency Vmax
BLGS[:,8]=all_busdata[:,11];
# emergency Vmin
BLGS[:,9]=all_busdata[:,12];
# base MVA
bMVA=hcat(all_gendata[:,1],all_gendata[:,15]);
for i=1:size(bMVA,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==bMVA[i,1]
            BLGS[j,17]=bMVA[i,2]
        end
    end
end
baseMVA = maximum(BLGS[:,17]);
# P demand
Pd=hcat(all_loaddata[:,1],all_loaddata[:,5]/baseMVA);
for i=1:size(Pd,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Pd[i,1]
            BLGS[j,10]=Pd[i,2]
        end
    end
end

# Q demand
Qd=hcat(all_loaddata[:,1],all_loaddata[:,6]/baseMVA);
for i=1:size(Qd,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Qd[i,1]
            BLGS[j,11]=Qd[i,2]
        end
    end
end
# P generation
Pg=hcat(all_gendata[:,1],all_gendata[:,2]/baseMVA);
for i=1:size(Pg,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Pg[i,1]
            BLGS[j,12]=Pg[i,2]
        end
    end
end
# Q generation
Qg=hcat(all_gendata[:,1],all_gendata[:,3]/baseMVA);
for i=1:size(Qg,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Qg[i,1]
            BLGS[j,13]=Qg[i,2]
        end
    end
end
# Q max
Qmax=hcat(all_gendata[:,1],all_gendata[:,4]/baseMVA);
for i=1:size(Qmax,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Qmax[i,1]
            BLGS[j,14]=Qmax[i,2]
        end
    end
end
# Q min
Qmin=hcat(all_gendata[:,1],all_gendata[:,5]/baseMVA);
for i=1:size(Qmin,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Qmin[i,1]
            BLGS[j,15]=Qmin[i,2]
        end
    end
end
# Voltage magnitude on generator
Vg=hcat(all_gendata[:,1],all_gendata[:,6]);
for i=1:size(Vg,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Vg[i,1]
            BLGS[j,16]=Vg[i,2]
        end
    end
end
# P max
Pmax=hcat(all_gendata[:,1],all_gendata[:,16]/baseMVA);
for i=1:size(Pmax,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Pmax[i,1]
            BLGS[j,18]=Pmax[i,2]
        end
    end
end
# P Min
Pmin=hcat(all_gendata[:,1],all_gendata[:,17]/baseMVA);
for i=1:size(Pmin,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Pmin[i,1]
            BLGS[j,19]=Pmin[i,2]
        end
    end
end
# Conductance Gs
Gs=hcat(all_shuntdata[:,1],all_shuntdata[:,3]/baseMVA);
for i=1:size(Gs,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Gs[i,1]
            BLGS[j,20]=Gs[i,2]
        end
    end
end
# Susceptance Bs
Bs=hcat(all_shuntdata[:,1],all_shuntdata[:,4]/baseMVA);
for i=1:size(Bs,1)
    for j=1:size(BLGS,1)
        if BLGS[j,1]==Bs[i,1]
            BLGS[j,21]=Bs[i,2]
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

Br=zeros(size(all_linedata,1)+size(all_trandata,1),7);
Br[:,7]=ones(size(all_linedata,1)+size(all_trandata,1));
Br[1:size(all_linedata,1),1:6]=all_linedata[:,1:6];
Br[size(all_linedata,1)+1:size(all_linedata,1)+size(all_trandata,1),1:2]=all_trandata[1:3,1:2];
Br[size(all_linedata,1)+1:size(all_linedata,1)+size(all_trandata,1),3:4]=all_trandata[1:3,19:20];
Br[size(all_linedata,1)+1:size(all_linedata,1)+size(all_trandata,1),6]=all_trandata[1:3,23];
Br[size(all_linedata,1)+1:size(all_linedata,1)+size(all_trandata,1),7]=all_trandata[1:3,22];
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
Ys = Br[:, 3] + im * Br[:, 4];
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
i = [1:NBr; 1:NBr];
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
RATE_A = Br[:,6]/baseMVA;

# Voltage limits
Vmax = BLGS[1,6]; # All values are the same
Vmin = BLGS[1,7]; # All values are the same

# Demand
Pd = BLGS[:,10];
Qd = BLGS[:,11];

# Starting value of voltage magnitude
#V0 = BLGS[:,4];
# Starting value of voltage phase angle
#δ0 = BLGS[:,5];

####################### Incidence Matrix [I x L] #############################
A = zeros(size(BLGS,1),size(Br,1));
for i=1:size(BLGS,1)
    for j=1:size(Br,1)
        if Br[j,1]==BLGS[i,1]
            A[i,j]=1
        elseif Br[j,2]==BLGS[i,1]
            A[i,j]=-1
        else A[i,j]=0
        end
    end
end

##################### Contingencies ##########################################
# All lines and transformers are online
aL = transpose(ones(size(Br,1),1));

# Contingency FROM buses
Cont_FR=float(contingency[:,3])[1,1];
# Contingency TO buses
Cont_TO=float(contingency[:,4])[1,1];
# Number of contingencies
Cont_size=size(contingency,1);
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
links
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
@variable(PSCOPF, V[i=1:NBus, k=1:NK])
#,Voltage phase angle
@variable(PSCOPF, δ[i=1:NBus, k=1:NK])
# Active power generation cost
@variable(PSCOPF, Cost, lowerbound=0);

#-------------------------- STARTING POINT -----------------------------------
# Warm start from PowerModels optimization: all lines are online
pk0 = transpose([0.37965 0.846428 0.0580179 0 0 1.10503 0 0.00319829 0 0 0 0 0 0]);
# Warm start from PowerModels optimization for the case with line 6-12 offline
pk1 = transpose([0.37965 0.848737 0.0580179 0 0 1.10503 0 0.00319829 0 0 0 0 0 0]);
p0 = hcat(pk0,pk1);
setvalue(p[1:NBus, 1:NK],p0)

# Warm start from PowerModels optimization: all lines are online
qk0 = transpose([0.0184817 0.17529 0.268335 0 0 -0.130702 0 0.0701339 0 0 0 0 0 0]);
# Warm start from PowerModels optimization for the case with line 6-12 offline
qk1 = transpose([0.0199731 0.178544 0.269567 0 0 -0.132471 0 0.0715449 0 0 0 0 0 0]);
q0 = hcat(qk0,qk1);
setvalue(q[1:NBus, 1:NK],q0)

# Warm start from PowerModels optimization: all lines are online
Vk0 = transpose([1.09559 1.07407 1.08165 1.08655 1.0801 1.06287 1.0914 1.06537 1.06967 1.07325 1.1 1.08877 1.1 1.07917]);
# Warm start from PowerModels optimization for the case with line 6-12 offline
Vk1 = transpose([1.09614 1.07352 1.08128 1.05878 1.07387 1.05993 1.09194 1.06584 1.06986 1.07354 1.1 1.08854 1.1 1.07858]);
V0 = hcat(Vk0,Vk1);
setvalue(V[1:NBus,1:NK], V0)

# Warm start from PowerModels optimization: all lines are online
δk0 = transpose([7.83E-34 -0.0460963 -0.00374554 0.0258115 0.0182263 -0.0401794 -0.0107868 -0.0934034 -0.0490856 -0.0268997 0.0447899 -0.0557694 -0.055299 -0.0596323]);
# Warm start from PowerModels optimization for the case with line 6-12 offline
δk1 = transpose([-4.89E-32 -0.047066 -0.00363043 0.00427936 0.010817 -0.0444117 -0.0107745 -0.0934403 -0.0491833 -0.0268099	0.0459912 -0.0567405 -0.05627 -0.0610723]);
δ0 = hcat(δk0,δk1);
setvalue(δ[1:NBus,1:NK], δ0)

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
@constraint(PSCOPF, ObjectiveFunction, (Cost == sum(a[i]*(p[i])^2 + b[i]*p[i] + c[i]  for i=1:NBus)));
@objective(PSCOPF, Min, Cost)

# Generation participation factor
#@constraint(PSCOPF, AGC[i=1:NBus, k=1:NK], p[i,k=2] == p[i,k=1]+part_factor[i].*sum(Pd[i]-p[i,k=1] for i=1:NBus));
#@constraint(PSCOPF, AGC[i=1:NBus, k=1:NK], p[i,k] == p[i]+part_factor[i].*sum(Pd[i]-p[i] for i=1:NBus));

# Print the model
#print(PSCOPF)

## resolution
solve(PSCOPF)

println("Objective value: ", getobjectivevalue(PSCOPF))
println("V: ", getvalue(V))
println("δ: ", getvalue(δ))
println("p: ", getvalue(p))
println("q: ", getvalue(q))

#sum(Pd)
#sum(BLGS[:,19])
