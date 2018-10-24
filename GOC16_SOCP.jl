#Grid Optimization Competition
#cd("C:/Users/Ильгиз/Documents/Документы_Ильгиз/Skoltech_2018/Grid Competition/Code/GOC_v15_dictionaries")
using CSV, JuMP;
using Ipopt;
using Gurobi, Mosek
println("Data reading")
    tic()
# Input files
#
rawFile = "IEEE14-1_powersystem.raw"

contFile = "IEEE14-1_contingency.csv"
genFile = "IEEE14-1_generator.csv"
#
#=
rawFile = "RTS96-1_powersystem.raw"
contFile = "RTS96-1_contingency.csv"
genFile = "RTS96-1_generator.csv"
=#
# Preventive Security Constrained Optimal Power Flow
#PSCOPF = Model(solver=IpoptSolver(print_level=0))
#PSCOPF = Model(solver=IpoptSolver())
PSCOPF = Model(solver=GurobiSolver())
# Need to define environment for mosek first: check mosek.jl package in julia. some command containing myhome
#env = Mosek.Env();
#PSCOPF = Model(solver=MosekSolver())
    toc()
####################################################################################
println("Defining External functions")
    tic()
### def.jl
type busData
  ID_ext :: Int64
  ID_int :: Int64
  Name :: String
  gen :: Array{Any,1}

  Vmax :: Float64
  Vmin :: Float64
  Pd :: Float64
  Qd :: Float64

  gsh :: Float64
  bsh :: Float64
end

type genData
  ID_ext :: Any
  Name :: Any
  Loc :: Int64
  ID_int :: Any

  cn :: Array{Int64,1}
  cParams :: Dict{Any,Any}

  Pmax :: Float64
  Pmin :: Float64
  Qmax :: Float64
  Qmin :: Float64

  alpha :: Float64
end

type branchData
  From :: Int64
  To :: Int64
  CKT :: Any
  ID_ext :: Any
  #revID :: Any

  r :: Float64
  x :: Float64
  bc :: Float64

  t :: Float64
  tau :: Float64
end

type fixedData
  baseMVA :: Float64

  busList :: Array{Any,1}
  busDList :: Dict{Any,Any}

  genList :: Array{Any,1}
  genDList :: Dict{Any,Any}

  brList :: Array{Any,1}
  brListSingle :: Array{Any,1}
  brDList :: Dict{Any,Any}
end

type uncertainData
  contList :: Array{Int64,1}
  contDList :: Dict{Any,Any}
end

type contingencyData
  ID :: Int64
  Type :: String
  Loc :: Array{Any,1}
end
###############################################################################
### preProc
function preProc(rawFile)
  rawData = readdlm(rawFile,',', skipblanks=false);

  # separate the data into different segments
  titleLine1 = rawData[1,:];
#  println("titleLine1: ",titleLine1);
  titleLine2 = rawData[2,:];
#  println("titleLine2: ",titleLine2);
  titleLine3 = rawData[3,:];
#  println("titleLine3: ",titleLine3);
  n,m = size(rawData);
  busStartL = 4; # first three lines are headers
  busEndL = 0
  loadStartL = 0
  loadEndL = 0
  shuntStartL = 0
  shuntEndL = 0
  genStartL = 0
  genEndL = 0
  branchStartL = 0
  branchEndL = 0
  transformerStartL = 0
  transformerEndL = 0
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

  return baseMVA,busSeg,loadSeg,shuntSeg,genSeg,branchSeg,transformerSeg;
end

########################################################################################
### busgenProc
function busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile)
bn,bm = size(busSeg);
gn,gm = size(genSeg);
sn,sm = size(shuntSeg);
ln,lm = size(loadSeg);

busSeg_int = collect(1:bn);
### BUS DATA
busList = [];
busDList = Dict();
busNo = 0;
for i in 1:bn
  # build the bus item
  busID_ext = busSeg[i,1];
  busID_int = busSeg_int[i,1];
  push!(busList,busID_ext);
  busName = string(busSeg[i,2]);
  busVmax = busSeg[i,10];
  busVmin = busSeg[i,11];
  busItem = busData(busID_ext,busID_int,busName,[],busVmax,busVmin,0,0,0,0);
  busDList[busID_ext] = busItem;
end

# process the shuntSeg
for i in 1:sn
  busID_ext = shuntSeg[i,1];
  busDList[busID_ext].gsh += shuntSeg[i,4]/baseMVA;
  busDList[busID_ext].bsh += shuntSeg[i,5]/baseMVA;
end

# process the loadSeg
for i in 1:ln
  busID_ext = loadSeg[i,1];
  busDList[busID_ext].Pd += loadSeg[i,6]/baseMVA;
  busDList[busID_ext].Qd += loadSeg[i,7]/baseMVA;
end

### GENERATOR DATA
# process the cost data from genFile
genCostData,genCostTitle = readdlm(genFile,',',header = true);
gcn,gcm = size(genCostData);
genCost = Dict();
for i in 1:gcn
  if !((genCostData[i,1],genCostData[i,2]) in keys(genCost))
    genCost[(genCostData[i,1],genCostData[i,2])] = Dict();
    genCost[(genCostData[i,1],genCostData[i,2])][genCostData[i,3]] = genCostData[i,4];
  else
    genCost[(genCostData[i,1],genCostData[i,2])][genCostData[i,3]] = genCostData[i,4];
  end
end

# process the genSeg
genList = [];
genDList = Dict();
for i in 1:gn
  busID_ext = genSeg[i,1];
  NGen = size(genSeg,1);
  genID_int = (1:NGen)[i];
  genLoc = busID_ext;
  genName = genSeg[i,2];
  genID = (genLoc,genName,genID_int);
  push!(genList,genID);
  push!(busDList[busID_ext].gen,genID);

  genQmax = genSeg[i,5]/baseMVA;
  genQmin = genSeg[i,6]/baseMVA;
  genPmax = genSeg[i,17]/baseMVA;
  genPmin = genSeg[i,18]/baseMVA;

  genCn = [];
  gencParams = Dict();
  genalpha = 1;
  for j in keys(genCost[(genLoc,genName)])
    if j != 9
      push!(genCn,j);
      gencParams[j] = genCost[(genLoc,genName)][j];
    else
      genalpha = genCost[(genLoc,genName)][j];
    end
  end

  genItem = genData(busID_ext,genName,genLoc,genID_int,genCn,gencParams,genPmax,genPmin,genQmax,genQmin,genalpha);
  genDList[genID] = genItem;
end
return busList,busDList,genList,genDList;
end

################################################################################
### branchProc
function branchProc(baseMVA,branchSeg,transformerSeg)
brn,brm = size(branchSeg);
brList = [];
brListSingle = [];
brDList = Dict();
for i in 1:brn
  brFrom = branchSeg[i,1];
  brTo = branchSeg[i,2];
  brCKT = branchSeg[i,3];
  branchID_ext = (brFrom,brTo,brCKT);

  brr = branchSeg[i,4];
  brx = branchSeg[i,5];
  brbc = branchSeg[i,6];
  brt = branchSeg[i,7]/baseMVA;

  brtau = 1;
  brItem = branchData(brFrom,brTo,brCKT,branchID_ext,brr,brx,brbc,brt,brtau);

  push!(brList,branchID_ext);
  push!(brListSingle,branchID_ext);
  brDList[branchID_ext] = brItem;
end

# process the transformerSeg
trn,trm = size(transformerSeg);
lineNo = 1;
while lineNo <= trn
  trFrom = transformerSeg[lineNo,1];
  trTo = transformerSeg[lineNo,2];
  if transformerSeg[lineNo,3] == 0
    lineNoNew = lineNo + 4;
  end
  trName = transformerSeg[lineNo,4];
  trID_ext = (trFrom,trTo,trName);

  trr = transformerSeg[lineNo+1,1];
  trx = transformerSeg[lineNo+1,2];
  trbc = 0;
  trt = transformerSeg[lineNo+2,4]/baseMVA;

  trtau = transformerSeg[lineNo+2,1]/transformerSeg[lineNo+3,1];
  trItem = branchData(trFrom,trTo,trName,trID_ext,trr,trx,trbc,trt,trtau);

  push!(brList,trID_ext);
  push!(brListSingle,trID_ext);
  brDList[trID_ext] = trItem;
  lineNo = lineNoNew;
end
return brList,brListSingle,brDList;
end
##############################################################################
### ContProc
function contProc(contFile)
contingency = CSV.read(contFile);
contData,contTitle = readdlm(contFile,',',header = true);
contn,contm = size(contData);
contList = [];
contDList = Dict();
for i in 1:contn
  contID = contData[i,1];
  push!(contList,contID);
  contType = contData[i,2];
  if (contType == "B")|(contType == "T")
    contLoc = [(contData[i,3],contData[i,4],contData[i,5])];
  else
    contLoc = [contData[i,3],contData[i,4]];
  end
  contItem = contingencyData(contID,contType,contLoc);
  contDList[contID] = contItem;
end
return contList,contDList, contingency;
end
    toc()
###############################################################################
### MyJulia1
println("Creating all dictionaties and structures")
    tic()
    rawFile
baseMVA,busSeg,loadSeg,shuntSeg,genSeg,branchSeg,transformerSeg = preProc(rawFile);
brList,brListSingle,brDList = branchProc(baseMVA,branchSeg,transformerSeg);
busList,busDList,genList,genDList = busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile);
contList,contDList, contingency = contProc(contFile);
fData = fixedData(baseMVA,busList,busDList,genList,genDList,brList,brListSingle,brDList);
uData = uncertainData(contList,contDList);
#############################################################################
### BuildMod
baseMVA = fData.baseMVA;
bList = fData.busList;
bData = fData.busDList;
gList = fData.genList;
gData = fData.genDList;
brList = fData.brList;
brData = fData.brDList;
S = uData.contList;
contData = uData.contDList;
    toc()

println("Defining rest of variables for model")
    tic()

####################################################################################### ALL UP GOES TO dProc function
# Constant term c in cost function
Dictcost = collect(gData[k].cParams for k in sort(collect(keys(gData))))
con = collect(Dictcost[k][:0] for k in sort(collect(keys(Dictcost))))
# Linear term b in cost function
b = collect(Dictcost[k][:1] for k in sort(collect(keys(Dictcost))))
# Quadratic term a in cost function
a = collect(Dictcost[k][:2] for k in sort(collect(keys(Dictcost))))
# Participation factor of generator
part_factor = collect(gData[k].alpha for k in sort(collect(keys(gData))))
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
Ys = 1 ./ (collect(brData[k].r for k in sort(collect(keys(brData)))) + im * collect(brData[k].x for k in sort(collect(keys(brData)))));
# line charging susceptance
Bc = collect(brData[k].bc for k in sort(collect(keys(brData))));
# transformer phase shirters
shift = zeros(size(brList,1),1);
# tap ratios
tap = collect(brData[k].tau for k in sort(collect(keys(brData)))).*exp.(im*pi/180*shift);
# admittance co-matrices
Ytt = Ys + im*Bc/2;

Yff = Ytt ./ (tap .* conj(tap));

Yft = - Ys ./ conj(tap);

Ytf = - Ys ./ tap;
# Shunt admittance
Ysh = (collect(bData[k].gsh for k in sort(collect(keys(bData)))) + im * collect(bData[k].bsh for k in sort(collect(keys(bData)))));

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
    toc()
println("Building the model")
    tic()
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
# Auxiliary variable
@variable(PSCOPF, c[i=1:NBus, j=1:NBus, k=1:NK])
# Auxiliary variable
@variable(PSCOPF, s[i=1:NBus, j=1:NBus, k=1:NK])
# Active power generation cost
#@variable(PSCOPF, Cost, lowerbound=0);
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
sum(Ω)
#--------------------------------- CONSTRAINTS -------------------------------
# Nodal active power balance
@constraint(PSCOPF, EqConstrP[i=1:NBus,k=1:NK], sum(Ω[i,g]*p[g,k] for g=1:NGen) - Pd[i] == sum(GL[i,j]*c[i,j] - BL[i,j]*s[i,j] for j=1:NBus));
# Nodal reactive power balance
@constraint(PSCOPF, EqConstrQ[i=1:NBus,k=1:NK], sum(Ω[i,g]*q[g,k] for g=1:NGen) - Qd[i] == sum(-BL[i,j]*c[i,j] - GL[i,j]*s[i,j] for j=1:NBus));

# Constraints on auxiliary variable c
@constraint(PSCOPF, SymmetryC[i=1:NBus, j=1:NBus; (GL[i,j] + BL[i,j]≠0.)], c[i,j] == c[j,i]);
# Constraints on auxiliary variable s
@constraint(PSCOPF, SymmetryS[i=1:NBus, j=1:NBus; ((GL[i,j] + BL[i,j]≠0.)&&i≠j)], s[i,j] == -s[j,i]);
# SOCP relaxation
@constraint(PSCOPF, SOCP[i=1:NBus, j=1:NBus; (GL[i,j] + BL[i,j]≠0.)], c[i,j]^2 + s[i,j]^2 + (0.5*(c[i,i]-c[j,j]))^2 ≤ (0.5*(c[i,i]+c[j,j]))^2);
# Vm => 0.9
@constraint(PSCOPF, LimitLCii[i=1:NBus,k=1:NK], 0.81 ≤ c[i,i,k]);
# Vm <= 1.1
@constraint(PSCOPF, LimitUCii[i=1:NBus,k=1:NK], c[i,i,k] ≤ 1.21);

# Active power flow
@constraint(PSCOPF,DefFlowP_ij[l=1:NBr,k=1:NK], pl[br_int[l,1],br_int[l,2],k] ==
(c[br_int[l,1],br_int[l,2],k]*GL[br_int[l,1],br_int[l,2]] - s[br_int[l,1],br_int[l,2],k]*BL[br_int[l,1],br_int[l,2]] - c[br_int[l,1],br_int[l,1],k]*GL[br_int[l,1],br_int[l,1]])*aL[k,l]);

@constraint(PSCOPF,DefFlowP_ji[l=1:NBr,k=1:NK], pl[br_int[l,2],br_int[l,1],k] ==
(c[br_int[l,2],br_int[l,1],k]*GL[br_int[l,2],br_int[l,1]] - s[br_int[l,2],br_int[l,1],k]*BL[br_int[l,2],br_int[l,1]] - c[br_int[l,2],br_int[l,2],k]*GL[br_int[l,2],br_int[l,2]])*aL[k,l]);

# Reactive power flow
@constraint(PSCOPF,DefFlowQ_ij[l=1:NBr,k=1:NK], ql[br_int[l,1],br_int[l,2],k] ==
(-s[br_int[l,1],br_int[l,2],k]*GL[br_int[l,1],br_int[l,2]] - c[br_int[l,1],br_int[l,2],k]*BL[br_int[l,1],br_int[l,2]] + c[br_int[l,1],br_int[l,1],k]*BL[br_int[l,1],br_int[l,1]])*aL[k,l]);

@constraint(PSCOPF,DefFlowQ_ji[l=1:NBr,k=1:NK], ql[br_int[l,2],br_int[l,1],k] ==
(-s[br_int[l,2],br_int[l,1],k]*GL[br_int[l,2],br_int[l,1]] - c[br_int[l,2],br_int[l,1],k]*BL[br_int[l,2],br_int[l,1]] + c[br_int[l,2],br_int[l,2],k]*BL[br_int[l,2],br_int[l,2]])*aL[k,l]);

# Line limits on Apparent Power
@NLconstraint(PSCOPF, FlowLimits_ij[l=1:NBr,k=1:NK], (pl[br_int[l,1],br_int[l,2],k]^2 + ql[br_int[l,1],br_int[l,2],k]^2) <= RATE_A[l]^2);
@NLconstraint(PSCOPF, FlowLimits_ji[l=1:NBr,k=1:NK], (pl[br_int[l,2],br_int[l,1],k]^2 + ql[br_int[l,2],br_int[l,1],k]^2) <= RATE_A[l]^2);

# Generator active power output upper bound
@constraint(PSCOPF, GenULimP[g=1:NGen,k=1:NK], p[g,k] <= Pmax[g]);
# Generator active power output lower bound
@constraint(PSCOPF, GenLLimP[g=1:NGen,k=1:NK], Pmin[g] <= p[g,k]);

# Generator reactive power output upper bound
@constraint(PSCOPF, GenULimQ[g=1:NGen,k=1:NK], q[g,k] <= Qmax[g]);
# Generator reactive power output lower bound
@constraint(PSCOPF, GenLLimQ[g=1:NGen,k=1:NK], Qmin[g] <= q[g,k]);

# Delta - powerfall constraint
@constraint(PSCOPF, Delta == sum(p[g,2] for g = 1:NGen) - sum(p[g,1] for g = 1:NGen));

# Post-contingency real power shortfall
total_part_factor = sum(part_factor);
@constraint(PSCOPF, PostCont[g=1:NGen], p[g,2] == p[g,1] + (part_factor[g]/total_part_factor).*Delta);

## Objective function
#@NLconstraint(PSCOPF, ObjectiveFunction, (Cost == sum(a[g]*(p[g,1])^2 + b[g]*p[g,1] + con[g]  for g=1:NGen)));
@objective(PSCOPF, Min, sum(a[g]*(p[g,1])^2 + b[g]*p[g,1] + con[g]  for g=1:NGen))

################################# AGC #########################################
# Voltage control only on generator buses (leads to unknown in Julia)
#@NLconstraint(PSCOPF, AGCLower[g=1:NGen, i=1:NBus], (q[g,2] - Qmin[g])*(sqrt(c[i,i,2]) - sqrt(c[i,i,1])) <= 0);
#@NLconstraint(PSCOPF, AGCUpper[g=1:NGen, i=1:NBus], (q[g,2] - Qmax[g])*(sqrt(c[i,i,2]) - sqrt(c[i,i,1])) <= 0);

## resolution
println("Model solution")
tic()
solve(PSCOPF)
toc()

# Deriving voltage phase values
δ_ij = atan2.(getvalue(s),getvalue(c))

# Deriving voltage magnitude values in base case
Vmb_list = [];
for i=1:NBus
    Vmb = sqrt(getvalue(c[i,i,1]))
    Vmb_list = [Vmb_list;Vmb]
end
Vmb_list

# Deriving voltage magnitude values in contingency case
Vmc_list = [];
for i=1:NBus
    Vmc = sqrt(getvalue(c[i,i,2]))
    Vmc_list = [Vmc_list;Vmc]
end
Vmc_list
################################ Deriving FROM and TO injections ##############
#=
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
Sf_k1 = Diagonal(Cf*Vk1c)*conj(Yf)*conj(Vk1c);

# Apparent power TO for pre-contingency case
St_k0 = Diagonal(Ct*Vk0c)*conj(Yt)*conj(Vk0c);
# Apparent power TO for post-contingency case
St_k1 = Diagonal(Ct*Vk1c)*conj(Yt)*conj(Vk1c);

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

# List of origin buses
orig_bus = collect(brDList[k].From for k in sort(collect(keys(brDList))));

# List of destination buses
dest_bus = collect(brDList[k].To for k in sort(collect(keys(brDList))));

# origin bus ID
pr5[1:NBr,3] = orig_bus; # 0. base case
pr5[NBr+1:2*NBr,3] = orig_bus; # 1. contingency case

# destination bus ID
pr5[1:NBr,4] = dest_bus; # 0. base case
pr5[NBr+1:2*NBr,4] = dest_bus; # 1. contingency case

# circuit ID
pr5[1:2*NBr,5] = repeat(["'BL'"],outer=[2*NBr]);

# real power in megawatts at origin
pr5[1:NBr,6] = Pf_k0*100;
pr5[NBr+1:2*NBr,6] = Pf_k1*100;

# reactive power in MVar at origin
pr5[1:NBr,7] = Qf_k0*100;
pr5[NBr+1:2*NBr,7] = Qf_k1*100;

# real power in megawatts at destination
pr5[1:NBr,8] = Pt_k0*100;
pr5[NBr+1:2*NBr,8] = Pt_k1*100;

# reactive power in MVar at destination
pr5[1:NBr,9] = Qt_k0*100;
pr5[NBr+1:2*NBr,9] = Qt_k1*100;
# Deleting row of contingency case
for i in size(pr5)[1]
    if pr5[i,1] == "1"
        for j in 1:NBr
            if aL[2,j] == 0
                pr5 = pr5[1:size(pr5,1) .!= j+NBr,: ]
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
    toc()
=#
