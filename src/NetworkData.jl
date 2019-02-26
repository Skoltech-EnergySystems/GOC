#=
This file containts a structure of The Network with all parameters of the network needed to create a JuMP optimization model

"""
General rules for new variables

    1. \undersocre{x}_i^j translates to xJ_i_min
    2. superscipt indeces become capital letters right afret a variable name
    3. subscript indeces become lower case letters after underscore
    4. all additions like overscore and underscore translate to _min and _max at the end of a      variable name
    5. all caligraphic variables go to cali*, i.e. caliA
"""
=#

"""
function to initialize switchedShunt data since it is easier to do line by line
"""
function temp_switchedShunt_init(Data, s)
    BL_len = 8 # the number of (Ni, Bi) pairs
    rows = size(Data)[1]
    bCS0 = Dict()
    bCS_max = Dict()
    bCS_min = Dict()
    for i = 1:rows
        bCS0[ Data[:I][i] ] = Data[:BINIT][i] ./ s
        BL = []
        for j = 1:BL_len
            n = Symbol("N$j")
            b = Symbol("B$j")
            push!(BL, Data[n][i] .* Data[b][i] ./ s)
        end
        NBL = findfirst(BL, 0) - 1
        bCS_max[ Data[:I][i] ] = sum([max(0,x) for x in BL])
        bCS_min[ Data[:I][i] ] = sum([min(0,x) for x in BL])
    end
    return bCS0, bCS_max, bCS_min
end

##############################################
# BUS
mutable struct Bus
    # caliI::Int
    # caliA::Int
    i::Int # i
    a::Int # a
    v0::Float64 # VM
    theta0::Float64 # VA
    v_max::Float64 # NVHI
    v_min::Float64 # NVLO
    vK_max::Float64 # EVHI
    vK_min::Float64 # EVLO

    # Bus() = new()
    Bus(Data) = bus_constr(new(), Data)
end
function bus_constr(B, Data)
    # B = new()
    # B.a = Data.bus[:AREA]
    # B.v0 = Data.bus[:VM]
    # B.theta0 = Data.bus[:VA]
    # B.v_max =Data.bus[:NVHI]
    # B.v_min = Data.bus[:NVLO]
    # B.vK_max = Data.bus[:EVHI]
    # B.vK_min = Data.bus[:EVLO]
    B.i = Data[1]
    B.a = Data[5]
    B.v0 = Data[8]
    B.theta0 = Data[9]
    B.v_max = Data[10]
    B.v_min = Data[11]
    B.vK_max = Data[12]
    B.vK_min = Data[13]
    return B
end

function Bus_init!(Data, PN)
  n = size(Data)[1]
  PN.BusList = []
  for i = 1:n
    B = Bus(Data[i,:])
    push!(PN.BusList, B)
    push!(PN.caliI, B.i)
    push!(PN.caliA, B.a)
  end
end

# LOAD
mutable struct Load
    pL::Float64
    qL::Float64

    Load(Data, s) = load_constr(new(), Data, s)
end
function load_constr(Ld, Data, s)
    Ld.pL = 0.0
    Ld.qL = 0.0
    if Data[3] == true
        Ld.pL += Data[6] ./ s
        Ld.qL += Data[7] ./ s
    end
    return Ld
end

function Load_init!(Data, PN)
  n = size(Data)[1]
  PN.LoadList = []
  for i = 1:n
    Ld = Load(Data[i,:], PN.s)
    push!(PN.LoadList, Ld)
  end
end


# FIXED SHUNT
mutable struct FixedShunt
    gFS::Float64
    bFS::Float64

    FixedShunt(Data) = fShunt_constr(new(), Data, s)
end
function fShunt_constr(fS, Data, s)
    fS.gFS = 0.0
    fS.bFS = 0.0
    if Data[3] == true
        fS.gFS += Data[4] ./ s
        fS.bFS += Data[5] ./ s
    end
    return fS
end

"""
in challenge 1 it is empty so I need to create an empty list otherwise i get an error "access to undefined reference"
"""
function fShunt_init!(Data, PN)
  n = size(Data)[1]
  PN.fShuntList = []
  for i = 1:n
    B = FixedShunt(Data[i,:], PN.s)
    push!(PN.fShuntList, B)
  end
end

# GENERATOR
mutable struct Generator
    # caliG::Int
    # G::Int
    i::Int
    id::Int
    g::Tuple
    p0::Float64
    q0::Float64
    q_max::Float64
    q_min::Float64
    p_max::Float64
    p_min::Float64
    stat::Bool
    # costs
    p::Array
    c::Array
    N::Int
    # participation factor
    alpha::Float64

    Generator(Data, s) = generator_constr(new(), Data, s)
end
function generator_constr(G, Data, s)
    G.i = Data[1]
    G.id = parse(Int,replace(strip(Data[2]), "'" => ""))
    G.g = (G.i, G.id)
    G.p0 = Data[3] ./ s
    G.q0 = Data[4] ./ s
    G.q_max = Data[5] ./ s
    G.q_min = Data[6] ./ s
    G.p_max = Data[17] ./ s
    G.p_min = Data[18] ./ s
    G.stat = Data[15]
    return G
end

function Generator_init!(Data, PN)
  n = size(Data)[1]
  PN.GeneratorList = []
  PN.gen_ind_I = Dict{Int, Int}()
  for i = 1:n
    G = Generator(Data[i,:], PN.s)
    push!(PN.GeneratorList, G)
    push!(PN.caliG, G.g)
    if G.stat == true
        push!(PN.G, G.g)
    end
    PN.gen_ind_I[G.g] = length(PN.GeneratorList)
  end
end

# LINE
mutable struct Line
    iO::Int
    iD::Int
    id::Int
    e::Tuple
    g::Float64
    b::Float64
    bCH::Float64
    R_max::Float64
    RK_max::Float64
    st::Bool
    Line(Data, s) = line_constr(new(), Data, s)
end
function line_constr(L, Data, s)
    L.iO = Data[1]
    L.iD = Data[2]
    L.id = parse(Int,replace(strip(Data[3]), "'" => ""))
    L.e = (L.iO, L.iD, L.id)
    L.g = Data[4] / (Data[4] ^ 2 + Data[5] ^ 2)
    L.b = Data[5] / (Data[4] ^ 2 + Data[5] ^ 2)
    L.bCH = Data[6]
    L.R_max = Data[7] / s
    L.RK_max = Data[9] / s
    L.st = Data[14]

    return L
end

function Line_init!(Data, PN)
  n = size(Data)[1]
  PN.LineList = []
  for i = 1:n
    L = Line(Data[i,:], PN.s)
    push!(PN.LineList, L)
    push!(PN.Eps, L.e)
    if L.st == true
        push!(PN.E, L.e)
    end
  end
end

# TRANSFORMER
mutable struct Transformer
    iO::Int
    iD::Int
    id::Int
    f::Tuple
    gM::Float64
    bM::Float64
    g::Float64
    b::Float64
    tau::Float64
    theta::Float64
    s_max::Float64
    sK_max::Float64
    st::Bool

    Transformer(Data, s) = transformer_constr(new(), Data, s)
end
function transformer_constr(T, Data, s)
    T.iO = Data[1]
    T.iD = Data[2]
    T.id = parse(Int,replace(strip(Data[4]), "'" => ""))
    T.f = (T.iO, T.iD, T.id)
    T.gM = Data[8]
    T.bM = Data[9]
    T.g = Data[22] / (Data[22] ^ 2 + Data[23] ^ 2)
    T.b = Data[23] / (Data[22] ^ 2 + Data[23] ^ 2)
    T.tau = Data[25] ./ Data[42]
    T.theta = Data[27] .* pi ./ 180
    T.s_max = Data[28] ./ s
    T.sK_max = Data[30] ./ s
    T.st = Data[12]

    return T
end

function Transformer_init!(Data, PN)
  n = size(Data)[1]
  PN.TransformerList = []
  step = 4 # these 4 lines should be contacatinated as one to put into the frame
  lines_len = [21, 3, 17, 2] # length of each of 4 lines. 43 attributes of transformer in total
  for i = 1:step:n
    full_line = []
    for j = 0:step-1
      full_line = [full_line; Data[i+j, 1:lines_len[j+1]]];
    end
    T = Transformer(full_line, PN.s)
    push!(PN.TransformerList, T)
    push!(PN.caliF, T.f)
    if T.st == true
        push!(PN.F, T.f)
    end
  end
end

# SWITCHED SHUNT
mutable struct SwitchedShunt
    i::Int
    bCS0::Float64
    bCS_max::Float64
    bCS_min::Float64
    stat::Bool

    SwitchedShunt(Data, s) = sShunt_constr(new(), Data, s)
end
function sShunt_constr(sS, Data, s)
    sS.i = Data[1]
    BL_len = 8
    BL = []
    st_ind = 11 # index of N1, starting the sequence of Ni, Bi
    if Data[4] == true
        sS.stat = true
        sS.bCS0 = Data[10]
        for j = st_ind:2:2*BL_len-1
            push!(BL, Data[j] * Data[j+1] / s)
        end
        NBL = findfirst(isequal(0), BL) - 1
        sS.bCS_max = sum([max(0,x) for x in BL])
        sS.bCS_min = sum([min(0,x) for x in BL])
    end
    return sS
end

function sShunt_init!(Data, PN)
  n = size(Data)[1]
  PN.sShuntList = []
  PN.get_shunt_index = Dict{Int, Int}()
  for i = 1:n
    S = SwitchedShunt(Data[i,:], PN.s)
    push!(PN.sShuntList, S)
    PN.get_shunt_index[S.i] = length(PN.sShuntList)
  end
end

#########################################
# main data descrbing network parameters
mutable struct PNetwork
    # base MVA
    s::Float16

    #bus data
    caliI::Set
    caliA::Set
    BusList::Array{Bus, 1}

    # load data
    LoadList::Array{Load, 1}

    # fixed shunt data
    fShuntList::Array{FixedShunt}

    # generator data
    caliG::Set
    G::Set
    GeneratorList::Array{Generator}
    gen_ind_I::Dict{Tuple, Int} # change to get_gener_index

    # line data BRANCH
    Eps::Set
    E::Set
    LineList::Array{Line}

    # transformer data
    caliF::Set
    F::Set
    TransformerList::Array{Transformer}

    # swtiched shunt
    sShuntList::Array{SwitchedShunt}
    get_shunt_index::Dict{Int, Int}

    PNetwork() = constr_network(new())
end
function constr_network(PN)
    PN.caliI = Set(Array{Int64, 1}())
    PN.caliA = Set(Array{Int64, 1}())
    PN.caliG = Set(Array{Tuple{Int64, Int64}, 1}())
    PN.G = Set(Array{Tuple{Int64, Int64}, 1}())
    PN.Eps = Set(Array{Tuple{Int64, Int64, Int64}, 1}())
    PN.E = Set(Array{Tuple{Int64, Int64, Int64}, 1}())
    PN.caliF = Set(Array{Tuple{Int64, Int64, Int64}, 1}())
    PN.F = Set(Array{Tuple{Int64, Int64, Int64}, 1}())

    return PN
end


# =====================================
# Costs data

function costs_init!(Data, PN)
    D0 = join(Data.genDispatch, Data.activeDispatch, on=[:TBL], kind=:left)
    D = join(D0, Data.linearCurve, on=[:CTBL], kind=:left)
    # println(names(D))
    for i = 1:size(D)[1]
        L = D[i,:]
        g = (L[:BUS][1], L[:GENID][1])
        j = PN.gen_ind_I[g]

        PN.GeneratorList[j].N = L[:NPAIRS][1]
        # return L[:XYmatrix][1][:,1]
        # println(names(L))
        # PN.GeneratorList[j].p = L[:XYmatrix][1][:,1] / PN.s
        PN.GeneratorList[j].p = L[:XYmatrix][:,1] / PN.s
        PN.GeneratorList[j].c = L[:XYmatrix][:,2] / PN.s
    end
    # return D
end


# OLD STRUCTURE TO DELETE
mutable struct PCosts
    caliH::Dict
    p_gh::Dict
    c_gh::Dict
    PCosts() = new()
end
